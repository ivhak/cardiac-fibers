// Copyright (C) 2022 Iver Håkonsen
//
// cardiac-fibers is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// cardiac-fibers is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// cardiac-fibers.  If not, see <https://www.gnu.org/licenses/>.
//
// Authors: Iver Håkonsen <hakonseniver@yahoo.no

#include <iostream>
#include <fstream>

#include "mfem.hpp"
#include "mfem/general/forall.hpp"

#include "cardiac-fibers.hpp"
#include "fem.hpp"
#include "util.hpp"
#include "calculus.hpp"

using namespace mfem;
using namespace util;

// Solve the Laplace equation
//
//      delta u = 0
//
// using the solver `solver`, storing the solution in `x`.
//
// This routine expects that
//   1. ParGridFunction `x` and its associated ParFiniteElementSpace has been
//      set up before hand.
//   2. The wanted solver, `solver`, and possible preconditioner, has been configured.
//
//
// Notes:
//
// - The boundary surfaces that are essential for a given solution are set in
//   `essential_boundaries`. For this use case these surfaces are the
//   epicardium, the left ventricle endocardium, the right ventricle and the
//   base. Each of these surfaces are expected to have an id/attribute
//   associated to it, which is dependent on the input mesh.
//
// - The essential boundary surfaces with boundary condition 0 are set
//   `zero_essential_boundaries`. For the boundary condition to be properly
//   enforced on a given surface, it also has to be set as an essential
//   boundary in `essential_boundaries`.
//
// - The boundary surfaces with boundary condition 1 are set in
//   `nonzero_essential_boundaries`. For the boundary condition to be properly
//   enforced on a given surface, it also has to be set as an essential boundary
//   in `essential_boundaries`.
//
// - The `apex` variable handles a corner case when solving the Laplace
//   equation from the base to the epicardium (psi_ab). In this special case, a
//   boundary condition of 0 needs to be enforced in a single vertex, namely
//   the apex. This is handled by setting `apex` to be the index of the apex
//   vertex in the mesh.
void laplace(
    ParGridFunction* x,
    Solver* solver,
    Array<int> &essential_boundaries,
    Array<int> &nonzero_essential_boundaries,
    Array<int> &zero_essential_boundaries,
    int apex,
    int verbose,
    std::ostream& tout,
    std::string tag,
    int rank)
{
    struct timespec t0, t1;

    std::string log_prefix = "[" + tag + "]:";

    // Determine the list of true (i.e. parallel conforming) essential boundary
    // dofs, defined by the boundary attributes marked as essential (Dirichlet)
    // and converting to a list of true dofs..
    Array<int> ess_tdof_list;

    ParFiniteElementSpace *fespace = x->ParFESpace();

    fespace->GetEssentialTrueDofs(essential_boundaries, ess_tdof_list);

    // Set up boundary conditions
    tracing::roctx_range_push("Setup boundary conditions");
    timing::tick(&t0);
    {
        // Initialize x with initial guess of zero, which satisfies the boundary
        // conditions.
        *x = 0.0;

        // Project the constant value 1.0 to all the essential boundaries marked as nonzero.
        ConstantCoefficient nonzero_bdr(1.0);
        x->ProjectBdrCoefficient(nonzero_bdr, nonzero_essential_boundaries);

        // Project the constant value 0.0 to all the essential boundaries marked as zero.
        ConstantCoefficient zero_bdr(0.0);
        x->ProjectBdrCoefficient(zero_bdr, zero_essential_boundaries);

        // For the laplacian involving the apex we need to treat the boundary
        // conditions a little different, as it is to be enforced on a single node
        // rather than a whole boundary surface. To do so, we make sure that the
        // apex is en essential true dof, and then we project the wanted value, in
        // this case 0.0, to only that node.
        if (apex >= 0) {
            // Initialize the internal data needed in the finite element space
            fespace->BuildDofToArrays();

            // Find the dof at the local apex vertex
            int apex_dof = fespace->GetLocalTDofNumber(apex);

            // This assertion should never fail, as we should have made sure that
            // the apex is in the form of a local vertex number.
            MFEM_ASSERT(apex_dof != -1, "No local DoF for the given local apex vertex");

            // Make sure the apex is in the list of essential true Dofs
            ess_tdof_list.Append(apex_dof);

            Array<int> node_disp(1);
            node_disp[0] = apex_dof;
            Vector node_disp_value(1);
            node_disp_value[0] = 0.0;

            VectorConstantCoefficient node_disp_value_coeff(node_disp_value);

            x->ProjectCoefficient(node_disp_value_coeff, node_disp);
        }

    }
    timing::tick(&t1);
    if (verbose >= 2 && rank == 0) {
        logging::timestamp(tout, log_prefix + " Setup boundary cond.", timing::duration(t0, t1), 3);
    }
    tracing::roctx_range_pop(); // Setup boundary conditions

    // Set up the parallel linear form b(.) which corresponds to the right-hand
    // side of the FEM linear system, which in this case is (1,phi_i) where
    // phi_i are the basis functions in fespace.
    tracing::roctx_range_push("Assemble RHS");
    timing::tick(&t0);

    ParLinearForm b(fespace);
    ConstantCoefficient zero(0.0);
    b.AddDomainIntegrator(new DomainLFIntegrator(zero));

    if (b.SupportsDevice()) {
        // This does not work on a simplex mesh, probably for the same reason that
        // AssemblyLevel::FULL does not work.
        if (verbose >= 3 && rank == 0)
            logging::info(std::cout, "MFEM_VERSION >= 4.5, using FastAssembly on RHS b.");
        b.UseFastAssembly(true);
    }

    b.Assemble();
    timing::tick(&t1);
    if (verbose >= 2 && rank == 0) {
        logging::timestamp(tout, log_prefix + " Assemble RHS", timing::duration(t0, t1), 3);
    }
    tracing::roctx_range_pop(); // Assemble RHS

    // Set up the parallel bilinear form a(.,.) on the finite element space
    // corresponding to the Laplacian operator -Delta, by adding the
    // Diffusion domain integrator.
    tracing::roctx_range_push("Assemble LHS");
    timing::tick(&t0);
    ParBilinearForm a(fespace);

    // Assemble the parallel bilinear form and the corresponding linear system.
    ConstantCoefficient one(1.0);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));

#if 0
    // AssemblyLevel::FULL makes Assemble() call some partial assembly
    // functions, and since partial assembly is not supported for simplex
    // elements, e.g. tetrahedrons, full assembly does not work.
    a.SetAssemblyLevel(AssemblyLevel::FULL);
#endif

    a.Assemble();
    timing::tick(&t1);
    if (verbose >= 2 && rank == 0) {
        logging::timestamp(tout, log_prefix + " Assemble LHS", timing::duration(t0, t1), 3);
    }
    tracing::roctx_range_pop(); // Assemble LHS

    // Form the linear system
    tracing::roctx_range_push("Form linear system");
    timing::tick(&t0);

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, *x, b, A, X, B);
    timing::tick(&t1);
    if (verbose >= 2 && rank == 0) {
        logging::timestamp(tout, log_prefix + " Form linear system", timing::duration(t0, t1), 3);
    }
    tracing::roctx_range_pop(); // Form linear system

    // Solve the linear system A X = B.
    tracing::roctx_range_push("Solve");
    timing::tick(&t0);

    solver->SetOperator(*A);
    solver->Mult(B, X);
    timing::tick(&t1, /* barrier */ true, /* device barrier */ true);
    if (verbose >= 2 && rank == 0) {
        logging::timestamp(tout, log_prefix + " Solve", timing::duration(t0, t1), 3);
    }
    tracing::roctx_range_pop(); // Solve

    // Recover the parallel grid function corresponding to X. This is the local
    // finite element solution on each processor.
    a.RecoverFEMSolution(X, b, *x);
}

int main(int argc, char *argv[])
{
    // Initialize MPI and HYPRE
    Mpi::Init();
    int rank = Mpi::WorldRank();
    int nranks = Mpi::WorldSize();
    Hypre::Init();

    Options opts;

    // Set program defaults
    opts.mesh_file = NULL;
    opts.par_mesh = false;
    opts.output_dir = ".";
    opts.time_to_file = false;

#if defined(MFEM_USE_HIP)
    opts.device_config = "hip";
#elif defined(MFEM_USE_CUDA)
    opts.device_config = "cuda";
#else
    opts.device_config = "cpu";
#endif

    opts.verbose = 0;
    opts.prescribed_apex = Vector(3);

    opts.save_paraview   = true;
    opts.save_mfem       = false;

    opts.save_laplacians = true;
    opts.save_gradients  = false;

    opts.save_partitioning = false;

    opts.alpha_endo =  60.0;
    opts.alpha_epi  = -60.0;
    opts.beta_endo  = -60.0;
    opts.beta_epi   =  60.0;

    opts.base_id = 1;
    opts.epi_id  = 2;
    opts.lv_id   = 3;
    opts.rv_id   = 4;

    opts.uniform_refinement = 0;

#if defined(MFEM_USE_HIP) || defined(MFEM_USE_CUDA)
    opts.gpu_ids = Array<int>(nranks);
    opts.gpu_ids = 0;

    opts.gpu_tuned_amg = false;
#endif

    opts.fibers_per_element = true;

    // Parse command-line options
    OptionsParser args(argc, argv);
    args.AddOption(&opts.verbose,
            "-v", "--verbose",
            "Set verbosity level:\n"
            "\t    1: Print timestamps for each major computation step.\n"
            "\t    2: Print more granular timestamps.\n"
            "\t    3: Print additional info.\n"
            "\t    4: Print BoomerAmg and solver output.\n"
            "\tIf the -t (--time-to-file) flag is passed, the output from level 1 and 2\n"
            "\twill be output to <OUT>/time.txt, where <OUT> is set with the -o (--out) flag.");

    args.AddOption(&opts.mesh_file,
            "-m", "--mesh",
            "Mesh file to use. "
            "See https://mfem.org/mesh-formats/ for a list of suppored formats.", true);

    args.AddOption(&opts.par_mesh,
            "-pm", "--par-mesh",
            "-npm", "--no-par-mesh",
            "Mesh is parallel, i.e, one file for each rank.\n"
            "\tEach rank will read the mesh <MESH_FILE>.<rank>, e.g., mesh.000000.");

    args.AddOption(&opts.prescribed_apex,
            "-a", "--apex",
            "Coordinate of apex, space separated list: 'x y z'.", true);

    args.AddOption(&opts.output_dir,
            "-o", "--out",
            "Directory for output files.");

    args.AddOption(&opts.time_to_file,
            "-t", "--time-to-file",
            "-nt","--no-time-to-file",
            "Output time log to <OUT>/time.txt rather than stdout, "
            "where <OUT> is set with the -o (--out) flag.");

    args.AddOption(&opts.device_config,
            "-d", "--device",
            "Device configuration string, see Device::Configure().");

    args.AddOption(&opts.save_paraview,
            "-p",  "--save-paraview",
            "-np", "--no-save-paraview",
            "Save data files for ParaView (paraview.org) visualization.");

    args.AddOption(&opts.save_mfem,
            "-s",  "--save-mfem",
            "-ns", "--no-save-mfem",
            "Save data files in the native MFEM format.");

    args.AddOption(&opts.save_laplacians,
            "-sl",  "--save-laplacians",
            "-nsl", "--no-save-laplacians",
            "Save the Laplacians phi_epi, phi_lv and phi_rv");

    args.AddOption(&opts.save_gradients,
            "-sg",  "--save-gradients",
            "-nsg", "--no-save-gradients",
            "Save the gradients of phi_epi, phi_lv, phi_rv and psi_ab");

    args.AddOption(&opts.save_partitioning,
            "-sp",  "--save-partitioning",
            "-nso", "--no-save-partitioning",
            "Save the mesh partitioning.");

    args.AddOption(&opts.alpha_endo,
            "-ao", "--alpha-endo",
            "Alpha angle in endocardium.");

    args.AddOption(&opts.alpha_epi,
            "-ai", "--alpha-epi",
            "Alpha angle in epicardium.");

    args.AddOption(&opts.beta_endo,
            "-bo", "--beta-endo",
            "Beta angle in endocardium.");

    args.AddOption(&opts.beta_epi,
            "-bi", "--beta-epi",
            "Beta angle in epicardium.");

    args.AddOption(&opts.base_id,
            "-base", "--base-id",
            "Id of the base surface.");

    args.AddOption(&opts.epi_id,
            "-epi", "--epi-id",
            "Id of the epicardium surface.");

    args.AddOption(&opts.lv_id,
            "-lv", "--lv-id",
            "Id of the left ventricle endocardium surface.");

    args.AddOption(&opts.rv_id,
            "-rv", "--rv-id",
            "Id of the right ventricle endocardium surface.\n"
            "\tSet to -1 if there is no right ventricle in the geometry, "
            "e.g. for a single ventricle geometry.");

    args.AddOption(&opts.uniform_refinement,
            "-u", "--uniform-refinement",
            "Perform n levels of uniform refinement on the mesh.");

    args.AddOption(&opts.fibers_per_element,
            "-fpe",  "--fibers-per-element",
            "-fpv",  "--fibers-per-vertex",
            "Calculate fibers in the L2 space (one fiber per element) "
            "or H1 (one fiber per vertex).");

#if defined (MFEM_USE_HIP) || defined (MFEM_USE_CUDA)
    args.AddOption(&opts.gpu_ids,
            "-gid",  "--gpu-ids",
            "Set the id of the GPU each rank should use");

    args.AddOption(&opts.gpu_tuned_amg,
            "-gamg",  "--gpu-tuned-amg",
            "-ngamg", "--no-gpu-tuned-amg",
            "Tune the BoomerAmg preconditioner for (hopefully) better GPU performance.");
#endif

    args.Parse();

    if (!args.Good()) {
        if (rank == 0) {
            args.PrintUsage(std::cout);
        }
        exit(1);
    }

    if (opts.verbose >= 3 && rank == 0) {
        args.PrintOptions(std::cout);
    }

    bool mesh_has_right_ventricle = opts.rv_id != -1;

    // Make sure the output directory exists if we are saving anything
    if (opts.save_mfem || opts.save_paraview || opts.time_to_file) {
        fs::mksubdir(opts.output_dir);
    }

    // If a filename has been passed with the -t flag, we output everything
    // that uses logging::timestamp and logging::marker to the filename instead
    // of stdout.
    std::ofstream f;
    std::string tout_file(opts.output_dir);
    tout_file += "/time.txt";
    std::ostream& tout = opts.time_to_file ? (f.open(tout_file.c_str()), f) : std::cout;
    if (opts.verbose >= 3 && rank == 0 && opts.time_to_file) {
        std::string msg = "Outputing timing log to " + tout_file;
        logging::info(std::cout, msg);
    }

    struct timespec t0, t1;     // For all intermediate timings
    struct timespec start, end; // Start to finish timing
    timing::tick(&start);

    // Enable hardware devices such as GPUs, and programming models such as
    // HIP, CUDA, OCCA, RAJA and OpenMP based on command line options.

#if defined (MFEM_USE_HIP) || defined (MFEM_USE_CUDA)
    Device device(opts.device_config, opts.gpu_ids[rank]);
#else
    Device device(opts.device_config);
#endif
    if (opts.verbose >= 3 && rank == 0) {
        device.Print();
    }

    int *partitioning = NULL;
    int global_num_elements = 0;

    // Load the mesh
    ParMesh *pmesh = NULL;
    if (!opts.par_mesh) {
        timing::tick(&t0);
        Mesh mesh(opts.mesh_file, 1, 1);
        timing::tick(&t1);

        global_num_elements = mesh.GetNE();

        if (opts.verbose >= 3 && rank == 0) {
            std::string msg = "Loaded meshfile '" + std::string(opts.mesh_file) + "' "
                 + "consisting of " + std::to_string(mesh.GetNV()) + " vertices "
                 + "and " + std::to_string(mesh.GetNE()) + " elements";
            logging::info(std::cout, msg);
        }
        if (opts.verbose && rank == 0)
            logging::timestamp(tout, "Mesh load", timing::duration(t0, t1));

        MFEM_VERIFY((mesh.MeshGenerator() & 0b1110) == 0,
                    "The mesh has non-simplex elements. Only tetrahedral meshes are supported");

        // Define a parallel mesh by a partitioning of the serial mesh.
        timing::tick(&t0);
        if (opts.save_partitioning) {
            partitioning = mesh.GeneratePartitioning(nranks, 1);
            pmesh = new ParMesh(MPI_COMM_WORLD, mesh, partitioning);
        } else {
            pmesh = new ParMesh(MPI_COMM_WORLD, mesh);
        }
        mesh.Clear();
        timing::tick(&t1);
        if (opts.verbose && rank == 0)
            logging::timestamp(tout, "Partition mesh", timing::duration(t0, t1));
    } else {
        timing::tick(&t0);
        std::stringstream fname;
        fname << opts.mesh_file  << "." << std::setw(6) << std::setfill('0') << rank;
        std::string filename = fname.str();

        if (opts.verbose >= 3 && rank == 0) {
            std::stringstream msg;
            msg << "The -pm flag was passed, loading mesh in parallel from "
                << "'" << opts.mesh_file << ".{000000 ... "
                << std::setw(6) << std::setfill('0') << nranks-1 << "}'";
            logging::info(std::cout, msg.str());

        }
        mfem::ifgzstream mesh_ifs(filename);
        pmesh = new ParMesh(MPI_COMM_WORLD, mesh_ifs, 1);
        timing::tick(&t1);
        if (opts.verbose && rank == 0) {
            logging::timestamp(tout, "Parallel mesh load", timing::duration(t0, t1));
        }
    }

    // Uniform refinement of the mesh
    if (opts.uniform_refinement > 0) {
        struct timespec uf0, uf1;
        timing::tick(&uf0);
        for (int i = 0; i < opts.uniform_refinement; i++)
            pmesh->UniformRefinement();
        timing::tick(&uf1);
        if (opts.verbose && rank == 0) {
            std::string msg = "Uniform refinement ("
                             + std::to_string(opts.uniform_refinement)
                             + ")";
            logging::timestamp(tout, msg, timing::duration(uf0, uf1));
        }

        if (opts.verbose >= 3) {
            int num_verts = pmesh->GetNV();
            std::vector<int> verts(nranks, 0);
            MPI_Allgather(&num_verts,    1, MPI_INT,
                           verts.data(), 1, MPI_INT, MPI_COMM_WORLD);

            int num_elems = pmesh->GetNE();
            std::vector<int> elems(nranks, 0);
            MPI_Allgather(&num_elems,   1, MPI_INT,
                          elems.data(), 1, MPI_INT, MPI_COMM_WORLD);

            if (rank == 0) {
                int sum_verts = std::accumulate(verts.begin(), verts.end(), 0);
                int sum_elems = std::accumulate(elems.begin(), elems.end(), 0);

                std::string msg =
                    "Performed " + std::to_string(opts.uniform_refinement)
                  + " rounds of uniform refinement on the mesh. Mesh now consists"
                  + " of " + std::to_string(sum_verts) + " vertices"
                  + " and " + std::to_string(sum_elems) + " elements.";

                logging::info(std::cout, msg);
            }
        }
    }



    // Each rank finds the vertex in its submesh that is closest to the
    // prescribed apex. The rank with the closest one sets up the boundary
    // conditions in that vertex.
    int apex = -1;
    {
        struct timespec begin_apex, end_apex;
        timing::tick(&begin_apex);
        timing::tick(&t0, /* barrier */ false);
        tracing::roctx_range_push("Find apex");

        if (opts.verbose >= 2 && rank == 0) {
            logging::marker(tout, "Find apex", 0);
        }

        Vector local_vertices;
        local_vertices.UseDevice(true);
        pmesh->GetVertices(local_vertices);
        const double *vert = local_vertices.Read();
        int n = pmesh->GetNV();

        Vector local_distances(n);
        double *dist = local_distances.Write();

        const double apex_x = opts.prescribed_apex[0];
        const double apex_y = opts.prescribed_apex[1];
        const double apex_z = opts.prescribed_apex[2];

        // The vertices we query from the parallel mesh are in SOA format:
        //     { x0,x1,x2,...,xn,y0,y1,y0,...,yn,z0,z1,z2,...,zn }
        // When running on device we also have to make sure that we synchronize
        // the device after each element is computed.

        timing::tick(&t1);
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "Setup", timing::duration(t0, t1), 1);
        }

        timing::tick(&t0);
        MFEM_FORALL(i, n, {
            dist[i]  = (vert[0*n+i]-apex_x) * (vert[0*n+i]-apex_x);
            MFEM_SYNC_THREAD;
            dist[i] += (vert[1*n+i]-apex_y) * (vert[1*n+i]-apex_y);
            MFEM_SYNC_THREAD;
            dist[i] += (vert[2*n+i]-apex_z) * (vert[2*n+i]-apex_z);
        });
        MFEM_DEVICE_SYNC;
        timing::tick(&t1);
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "Calculate distances", timing::duration(t0, t1), 1);
        }

        timing::tick(&t0);
        double min_distance = std::numeric_limits<double>::max();

        const double *host_distance = local_distances.HostRead();

        // TODO (ivhak): This could (and should) probably be implemented as a
        // parallel reduction and be run on the GPU as well, but this work for
        // now...
        for (int i = 0; i < pmesh->GetNV(); i++) {
            if (host_distance[i] < min_distance) {
                apex = i;
                min_distance = host_distance[i];
            }
        }
        timing::tick(&t1);
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "Find local minimum", timing::duration(t0, t1), 1);
        }

        timing::tick(&t0);
        struct double_int local_min = {.val=min_distance, .rank=rank};
        struct double_int global_min;

        MPI_Allreduce(&local_min, &global_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
        if (rank != global_min.rank) {
            apex = -1;
        } else if (opts.verbose >= 3) {
            const double *found_apex = pmesh->GetVertex(apex);
            std::string msg = "Rank " + std::to_string(rank)
                            + " found the vertex closest to the prescribed apex ("
                            + std::to_string(opts.prescribed_apex[0]) +", "
                            + std::to_string(opts.prescribed_apex[1]) +", "
                            + std::to_string(opts.prescribed_apex[2]) +") at ("
                            + std::to_string(found_apex[0]) + ", "
                            + std::to_string(found_apex[1]) + ", "
                            + std::to_string(found_apex[2]) + ")";
            logging::info(std::cout, msg);
        }
        timing::tick(&t1);
        timing::tick(&end_apex, /* barrier*/ false);
        tracing::roctx_range_pop();
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "Find global minimum", timing::duration(t0, t1), 1);
        }

        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "Total", timing::duration(begin_apex, end_apex), 1, '=');
        } else if (opts.verbose && rank == 0) {
            logging::timestamp(tout, "Find apex",
                               timing::duration(begin_apex, end_apex), 0);
        }

    }


    // Setup the linear solver and preconditioner.
    timing::tick(&t0);
    HypreBoomerAMG *prec = NULL;
    {
        prec = new HypreBoomerAMG;
        prec->SetPrintLevel(opts.verbose >= 4 ? 1 : 0);

#ifdef MFEM_USE_CUDA_OR_HIP
        // FIXME (ivhak): Currently, these setting end up giving the wrong
        // solution to grad_psi_ab, probably caused by enforcing a boundary
        // condition in a single vertex in psi_ab.
        if (opts.gpu_tuned_amg) {
            prec->SetCoarsening(8);           // PMIS
            prec->SetCycleNumSweeps(1,1);     // 1 sweep on the up and down cycle
            prec->SetInterpolation(17);       // extended+i, matrix-matrix
            prec->SetAggressiveCoarsening(1); // Number of levels of aggressive coarsening
            prec->SetStrengthThresh(0.5);
            prec->SetRelaxType(7);            // weighted Jacobi
        }
#endif
    }

    Solver *solver;
    HyprePCG *cg = new HyprePCG(MPI_COMM_WORLD);
    cg->SetTol(1e-20);
    cg->SetMaxIter(2000);
    cg->SetPrintLevel(opts.verbose >= 4 ? 2 : 0);
    if (prec) cg->SetPreconditioner(*prec);
    solver = cg;

    timing::tick(&t1);
    if (opts.verbose && rank == 0) {
        logging::timestamp(tout, "Setup preconditioner and solver", timing::duration(t0,t1));
    }


    if (opts.verbose && rank == 0) {
        logging::marker(tout, "Run LDRB algorithm");
    }

    struct timespec begin_laplace, end_laplace;
    timing::tick(&begin_laplace);
    if (opts.verbose >= 2 && rank == 0) {
        logging::marker(tout, "[LDRB] Compute Laplacians", 1);
    }
    // Time the fiber orientation calculations
    struct timespec begin_fiber, end_fiber;
    timing::tick(&begin_fiber);
    timing::tick(&t0, /* barrier */ false);


    // Set up the finite element collection. We use  first order H1-conforming finite elements.
    H1_FECollection h1_fec(1, pmesh->Dimension());

    // Set up two finite element spaces: one for scalar values (1D) , and one for 3d vectors
    ParFiniteElementSpace fespace_scalar_h1(pmesh, &h1_fec);

    ParGridFunction *x_phi_epi = new ParGridFunction(&fespace_scalar_h1);
    ParGridFunction *x_phi_lv  = new ParGridFunction(&fespace_scalar_h1);
    ParGridFunction *x_phi_rv  = new ParGridFunction(&fespace_scalar_h1);
    ParGridFunction *x_psi_ab  = new ParGridFunction(&fespace_scalar_h1);

    x_phi_epi->UseDevice(true);
    x_phi_lv->UseDevice(true);
    x_phi_rv->UseDevice(true);
    x_psi_ab->UseDevice(true);

    int nattr = pmesh->bdr_attributes.Max();
    Array<int> ess_bdr(nattr);          // Essential boundaries
    Array<int> nonzero_ess_bdr(nattr);  // Essential boundaries with value 1.0
    Array<int> zero_ess_bdr(nattr);     // Essential boundaries with value 0.0

    timing::tick(&t1);
    if (opts.verbose >= 2 && rank == 0) {
        logging::timestamp(tout, "[Laplacians]: Setup", timing::duration(t0, t1), 2);
    }
    // Solve the Laplace equation from EPI (1.0) to (LV_ENDO union RV_ENDO) (0.0)
    {
        if (opts.verbose >= 2 && rank == 0) {
            logging::marker(tout, "[Laplacians]: phi_epi", 2);
        }
        timing::tick(&t0);

        ess_bdr = 0;
        ess_bdr[opts.epi_id-1] = 1;
        ess_bdr[opts.lv_id -1] = 1;
        if (mesh_has_right_ventricle)
            ess_bdr[opts.rv_id-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.epi_id-1] = 1;

        zero_ess_bdr = 0;
        zero_ess_bdr[opts.lv_id-1] = 1;
        if (mesh_has_right_ventricle)
            zero_ess_bdr[opts.rv_id-1] = 1;

        tracing::roctx_range_push("laplace x_phi_epi");

        laplace(x_phi_epi, solver,
                ess_bdr, nonzero_ess_bdr, zero_ess_bdr,
                -1, opts.verbose, tout, "phi_epi", rank);

        tracing::roctx_range_pop();

        timing::tick(&t1);
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "[phi_epi]: Total", timing::duration(t0, t1), 3, '=');
        }
    }


    // Solve the Laplace equation from LV_ENDO (1.0) to (RV_ENDO union EPI) (0.0)
    {
        if (opts.verbose >= 2 && rank == 0) {
            logging::marker(tout, "[Laplacians]: phi_lv", 2);
        }
        timing::tick(&t0);

        ess_bdr = 0;
        ess_bdr[opts.epi_id-1] = 1;
        ess_bdr[opts.lv_id -1] = 1;
        if (mesh_has_right_ventricle)
            ess_bdr[opts.rv_id-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.lv_id-1] = 1;

        zero_ess_bdr = 0;
        zero_ess_bdr[opts.epi_id-1] = 1;
        if (mesh_has_right_ventricle)
            zero_ess_bdr[opts.rv_id-1] = 1;

        tracing::roctx_range_push("laplace x_phi_lv");

        laplace(x_phi_lv, solver,
                ess_bdr, nonzero_ess_bdr, zero_ess_bdr,
                -1, opts.verbose, tout, "phi_lv", rank);

        tracing::roctx_range_pop();

        timing::tick(&t1);
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "[phi_lv]: Total", timing::duration(t0, t1), 3, '=');
        }
    }


    // Solve the Laplace equation from RV_ENDO (1.0) to (LV_ENDO union EPI) (0.0)
    if (mesh_has_right_ventricle) {
        if (opts.verbose >= 2 && rank == 0) {
            logging::marker(tout, "[Laplacians]: phi_rv", 2);
        }
        timing::tick(&t0);

        ess_bdr = 0;
        ess_bdr[opts.epi_id    -1] = 1;
        ess_bdr[opts.lv_id-1] = 1;
        ess_bdr[opts.rv_id-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.rv_id-1] = 1;

        zero_ess_bdr = 0;
        zero_ess_bdr[opts.epi_id-1] = 1;
        zero_ess_bdr[opts.lv_id -1] = 1;

        tracing::roctx_range_push("laplace x_phi_rv");

        laplace(x_phi_rv, solver,
                ess_bdr, nonzero_ess_bdr, zero_ess_bdr,
                -1, opts.verbose, tout, "phi_rv", rank);

        tracing::roctx_range_pop();
        timing::tick(&t1);
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "[phi_rv]: Total", timing::duration(t0, t1), 3, '=');
        }
    } else {
        *x_phi_rv = 0.0;
    }


    // Solve the Laplace equation from BASE (1.0) to APEX (0.0)
    {
        if (opts.verbose >= 2 && rank == 0) {
            logging::marker(tout, "[Laplacians]: psi_ab", 2);
        }
        timing::tick(&t0);

        ess_bdr = 0;
        ess_bdr[opts.base_id-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.base_id-1] = 1;

        zero_ess_bdr = 0;

        tracing::roctx_range_push("laplace x_psi_ab");

        laplace(x_psi_ab, solver,
                ess_bdr, nonzero_ess_bdr, zero_ess_bdr,
                apex, opts.verbose, tout, "psi_ab", rank);

        tracing::roctx_range_pop();

        timing::tick(&t1);
        if (opts.verbose >= 2 && rank == 0 ) {
            logging::timestamp(tout, "[psi_ab]: Total", timing::duration(t0, t1), 3, '=');
        }
    }

    // We are done using the preconditioner and solver
    delete solver;
    delete prec;

    timing::tick(&end_laplace);
    if (opts.verbose >= 2 && rank == 0) {
        logging::timestamp(tout, "[Laplacians]: Total", timing::duration(begin_laplace, end_laplace), 2, '=');
    } else if (opts.verbose && rank == 0) {
        logging::timestamp(tout, "Compute Laplacians",
                           timing::duration(begin_laplace, end_laplace), 1);
    }

    // If we want the fibers to be in L2, we first have to project the
    // Laplacians from H1 to L2. Then we can compute the gradients in L2.
    if (opts.verbose >= 2 && rank == 0) {
        logging::marker(tout, "[LDRB]: Compute gradients", 1);
    }
    struct timespec begin_grad, end_grad;
    timing::tick(&begin_grad);
    timing::tick(&t0, /* barrier */ false);

    L2_FECollection l2_fec(0, pmesh->Dimension());
    ParFiniteElementSpace fespace_vector_l2(pmesh, &l2_fec, 3, Ordering::byVDIM);

    // Setup the things needed for gradient computation. The tables will also
    // be used for the H1 to L2 projection when calculating the fibers in L2,
    // and for the gradient interpolation when calculating fibers in H1.
    //
    // NOTE: These device side arrays (vertices, h1_I, h1_J, l2_I, l2_J) are
    // not freed by MFEM's memory manager until the end of this scope, i.e.,
    // after the fiber computations. However, they are only needed until we are
    // done with the gradients/projections. Make sure they are manually freed
    // after that.
    Vector local_vertices;
    pmesh->GetVertices(local_vertices);
    const double *vert = local_vertices.Read();

    const Table& h1_element_to_dof = fespace_scalar_h1.GetElementToDofTable();
    const Table& l2_element_to_dof = fespace_vector_l2.GetElementToDofTable();

    const int *h1_I = h1_element_to_dof.ReadI();
    const int *h1_J = h1_element_to_dof.ReadJ();

    const int *l2_I = l2_element_to_dof.ReadI();
    const int *l2_J = l2_element_to_dof.ReadJ();

    const int num_elements = pmesh->GetNE();
    const int num_vertices = pmesh->GetNV();

    ParGridFunction *grad_phi_epi, *grad_phi_lv, *grad_phi_rv, *grad_psi_ab;

    grad_phi_epi = new ParGridFunction(&fespace_vector_l2);
    grad_phi_lv  = new ParGridFunction(&fespace_vector_l2);
    grad_phi_rv  = new ParGridFunction(&fespace_vector_l2);
    grad_psi_ab  = new ParGridFunction(&fespace_vector_l2);

    grad_phi_epi->UseDevice(true);
    grad_phi_lv->UseDevice(true);
    grad_phi_rv->UseDevice(true);
    grad_psi_ab->UseDevice(true);

    timing::tick(&t1);
    if (opts.verbose >= 2 && rank == 0) {
        logging::timestamp(tout, "[Gradients]: Setup", timing::duration(t0, t1), 2);
    }

    // Gradient of phi_epi
    {
        timing::tick(&t0);
        tracing::roctx_range_push("grad_phi_epi");

        const double *epi_vals = x_phi_epi->Read();
        double *epi_grads = grad_phi_epi->Write();

        compute_gradient(epi_grads, epi_vals, vert,
                         num_elements, num_vertices,
                         h1_I, h1_J, l2_I, l2_J);
        tracing::roctx_range_pop();
        timing::tick(&t1, /* barrier */ true, /* device barrier*/ true);
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "[Gradients]: grad_phi_epi", timing::duration(t0, t1), 2);
        }
    }

    // Gradient of phi_rv
    {
        timing::tick(&t0);
        tracing::roctx_range_push("grad_phi_lv");

        const double *lv_vals = x_phi_lv->Read();
        double *lv_grads = grad_phi_lv->Write();

        compute_gradient(lv_grads, lv_vals, vert,
                         num_elements, num_vertices,
                         h1_I, h1_J, l2_I, l2_J);
        tracing::roctx_range_pop();
        timing::tick(&t1, /* barrier */ true, /* device barrier*/ true);
        if (opts.verbose >= 2&& rank == 0) {
            logging::timestamp(tout, "[Gradients]: grad_phi_lv", timing::duration(t0, t1), 2);
        }
    }

    // Gradient of phi_lv
    if (mesh_has_right_ventricle) {
        timing::tick(&t0);
        tracing::roctx_range_push("grad_phi_rv");

        const double *rv_vals = x_phi_rv->Read();
        double *rv_grads = grad_phi_rv->Write();

        compute_gradient(rv_grads, rv_vals, vert,
                         num_elements, num_vertices,
                         h1_I, h1_J, l2_I, l2_J);
        tracing::roctx_range_pop();
        timing::tick(&t1, /* barrier */ true, /* device barrier*/ true);
        if (opts.verbose >= 2&& rank == 0) {
            logging::timestamp(tout, "[Gradients]: grad_phi_rv", timing::duration(t0, t1), 2);
        }
    } else {
        *grad_phi_rv = 0.0;
    }

    // Gradient of psi_ab
    {
        timing::tick(&t0);
        tracing::roctx_range_push("grad_psi_ab");

        const double *ab_vals = x_psi_ab->Read();
        double *ab_grads = grad_psi_ab->Write();

        compute_gradient(ab_grads, ab_vals, vert,
                         num_elements, num_vertices,
                         h1_I, h1_J, l2_I, l2_J);
        tracing::roctx_range_pop();
        timing::tick(&t1, /* barrier */ true, /* device barrier*/ true);
        timing::tick(&t1);
        if (opts.verbose >= 2&& rank == 0) {
            logging::timestamp(tout, "[Gradients]: grad_psi_ab", timing::duration(t0, t1), 2);
        }
    }
    timing::tick(&end_grad);
    if (opts.verbose >= 2 && rank == 0) {
        logging::timestamp(tout, "[Gradients]: Total", timing::duration(begin_grad, end_grad), 2, '=');
    } else if (opts.verbose && rank == 0) {
        logging::timestamp(tout, "Compute gradients", timing::duration(begin_grad, end_grad), 1);
    }

    // If we want one fiber per element, i.e. in L2, then we have to project
    // the Laplacians from H1 to L2. The gradients are already in the correct space.
    //
    // If we want one fiber per vertex, i.e. in H1, then we have to interpolate
    // the gradients from H1 to L2. the Laplacians are already in the correct
    // space
    if (opts.fibers_per_element) {
        if (opts.verbose >= 2 && rank == 0) {
            logging::marker(tout, "[LDRB]: Project Laplacians H1 -> L2", 1);
        }
        struct timespec begin_proj, end_proj;
        timing::tick(&begin_proj);
        tracing::roctx_range_push("Project Laplacians H1 -> L2");
        tracing::roctx_range_push("Setup");

        timing::tick(&t0, /* barrier */ false);
        ParFiniteElementSpace *fespace_scalar_l2 = new ParFiniteElementSpace(pmesh, &l2_fec);
        ParGridFunction *x_phi_epi_l2 = new ParGridFunction(fespace_scalar_l2);
        ParGridFunction *x_phi_lv_l2 = new ParGridFunction(fespace_scalar_l2);
        ParGridFunction *x_phi_rv_l2 = new ParGridFunction(fespace_scalar_l2);

        x_phi_epi_l2->UseDevice(true);
        x_phi_lv_l2->UseDevice(true);
        x_phi_rv_l2->UseDevice(true);

        timing::tick(&t1);
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "[Projection]: Setup", timing::duration(t0, t1), 2);
        }
        tracing::roctx_range_pop();

        // Project x_phi_epi, x_phi_lv and x_phi_rv into L2, then free the old
        // solutions and set the solution pointers to the projected ones.

        // Epicardium
        {
            timing::tick(&t0);
            const double *x_phi_epi_h1_vals = x_phi_epi->Read();
            double *x_phi_epi_l2_vals = x_phi_epi_l2->Write();
            project_h1_to_l2(x_phi_epi_l2_vals, x_phi_epi_h1_vals, num_elements,
                              h1_I, h1_J, l2_I, l2_J);
            timing::tick(&t1, /* barrier */ true, /* device barrier*/ true);
            if (opts.verbose >= 2 && rank == 0) {
                logging::timestamp(tout, "[Projection]: x_phi_epi", timing::duration(t0, t1), 2);
            }
        }

        // Left ventricle endocardium
        {
            timing::tick(&t0);
            const double *x_phi_lv_h1_vals = x_phi_lv->Read();
            double *x_phi_lv_l2_vals = x_phi_lv_l2->Write();
            project_h1_to_l2(x_phi_lv_l2_vals, x_phi_lv_h1_vals, num_elements,
                              h1_I, h1_J, l2_I, l2_J);
            timing::tick(&t1, /* barrier */ true, /* device barrier*/ true);
        }
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "[Projection]: x_phi_lv", timing::duration(t0, t1), 2);
        }

        // Right ventricle endocardium
        if (mesh_has_right_ventricle) {
            x_phi_rv_l2->UseDevice(true);
            timing::tick(&t0);
            {
                const double *x_phi_rv_h1_vals = x_phi_rv->Read();
                double *x_phi_rv_l2_vals = x_phi_rv_l2->Write();
                project_h1_to_l2(x_phi_rv_l2_vals, x_phi_rv_h1_vals, num_elements,
                                  h1_I, h1_J, l2_I, l2_J);
                timing::tick(&t1, /* barrier */ true, /* device barrier*/ true);
            }
            timing::tick(&t1);
            if (opts.verbose >= 2 && rank == 0) {
                logging::timestamp(tout, "[Projection]: x_phi_rv", timing::duration(t0, t1), 2);
            }
        } else {
            *x_phi_rv_l2 = 0.0;
        }

        delete x_phi_epi;
        delete x_phi_lv;
        delete x_phi_rv;

        x_phi_epi = x_phi_epi_l2;
        x_phi_lv  = x_phi_lv_l2;
        x_phi_rv  = x_phi_rv_l2;

        tracing::roctx_range_pop();
        timing::tick(&end_proj);
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "[Projection]: Total", timing::duration(begin_proj, end_proj), 2, '=');
        } else if (opts.verbose && rank == 0) {
            logging::timestamp(tout, "Project Laplacians from H1 -> L2",
                               timing::duration(begin_proj, end_proj), 1);
        }
    } else {
        // If we want the fibers on a per vertex basis (i.e. in H1), we need to
        // interpolate the gradients into H1 rather than projecting the
        // Laplacians into L2.


        if (opts.verbose >= 2 && rank == 0) {
            logging::marker(tout, "[LDRB]: Interpolate gradients L2 -> H1", 1);
        }
        struct timespec begin_interp, end_interp;
        timing::tick(&begin_interp);
        tracing::roctx_range_push("Interpolate gradients L2 -> H1");

        timing::tick(&t0);
        ParFiniteElementSpace *fespace_vector_h1 = new ParFiniteElementSpace(pmesh, &h1_fec, 3, Ordering::byVDIM);
        Table *vertex_to_element = pmesh->GetVertexToElementTable();

        const int *v2e_I = vertex_to_element->ReadI();
        const int *v2e_J = vertex_to_element->ReadJ();
        ParGridFunction *grad_phi_epi_h1 = new ParGridFunction(fespace_vector_h1);
        ParGridFunction *grad_phi_lv_h1 = new ParGridFunction(fespace_vector_h1);
        ParGridFunction *grad_phi_rv_h1 = new ParGridFunction(fespace_vector_h1);
        ParGridFunction *grad_psi_ab_h1 = new ParGridFunction(fespace_vector_h1);

        grad_phi_epi_h1->UseDevice(true);
        grad_phi_lv_h1->UseDevice(true);
        grad_phi_rv_h1->UseDevice(true);
        grad_psi_ab->UseDevice(true);

        timing::tick(&t1);
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "[Interpolate]: Setup", timing::duration(t0, t1), 2);
        }

        {
            timing::tick(&t0);
            double *grads_epi_h1 = grad_phi_epi_h1->Write();
            const double *grads_epi_l2 = grad_phi_epi->Read();

            interpolate_gradient_to_h1(grads_epi_h1, grads_epi_l2, num_vertices,
                                       v2e_I, v2e_J, l2_I, l2_J);
            timing::tick(&t1, /* barrier */ true, /* device barrier*/ true);
            if (opts.verbose >= 2 && rank == 0) {
                logging::timestamp(tout, "[Interpolate]: grad_phi_epi", timing::duration(t0, t1), 2);
            }
        }

        {
            timing::tick(&t0);
            double *grads_lv_h1 = grad_phi_lv_h1->Write();
            const double *grads_lv_l2 = grad_phi_lv->Read();

            interpolate_gradient_to_h1(grads_lv_h1, grads_lv_l2, num_vertices,
                                       v2e_I, v2e_J, l2_I, l2_J);
            timing::tick(&t1, /* barrier */ true, /* device barrier*/ true);
            if (opts.verbose >= 2 && rank == 0) {
                logging::timestamp(tout, "[Interpolate]: grad_phi_lv", timing::duration(t0, t1), 2);
            }
        }

        if (mesh_has_right_ventricle) {
            timing::tick(&t0);
            double *grads_rv_h1 = grad_phi_rv_h1->Write();
            const double *grads_rv_l2 = grad_phi_rv->Read();

            interpolate_gradient_to_h1(grads_rv_h1, grads_rv_l2, num_vertices,
                                       v2e_I, v2e_J, l2_I, l2_J);
            timing::tick(&t1, /* barrier */ true, /* device barrier*/ true);
            if (opts.verbose >= 2 && rank == 0) {
                logging::timestamp(tout, "[Interpolate]: grad_phi_rv", timing::duration(t0, t1), 2);
            }
        } else {
            *grad_phi_rv_h1 = 0.0;
        }

        {
            timing::tick(&t0);
            double *grads_ab_h1 = grad_psi_ab_h1->Write();
            const double *grads_ab_l2 = grad_psi_ab->Read();

            interpolate_gradient_to_h1(grads_ab_h1, grads_ab_l2, num_vertices,
                                       v2e_I, v2e_J, l2_I, l2_J);
            timing::tick(&t1, /* barrier */ true, /* device barrier*/ true);
            if (opts.verbose >= 2 && rank == 0) {
                logging::timestamp(tout, "[Interpolate]: grad_psi_ab", timing::duration(t0, t1), 2);
            }
        }

        delete grad_phi_epi;
        delete grad_phi_lv;
        delete grad_phi_rv;
        delete grad_psi_ab;

        grad_phi_epi = grad_phi_epi_h1;
        grad_phi_lv = grad_phi_lv_h1;
        grad_phi_rv = grad_phi_rv_h1;
        grad_psi_ab = grad_psi_ab_h1;

        delete vertex_to_element;
        timing::tick(&end_interp);
        if (opts.verbose >= 2 && rank == 0) {
            logging::timestamp(tout, "[Interpolate]: Total", timing::duration(begin_interp, end_interp), 2, '=');
        } else if (opts.verbose && rank == 0) {
            logging::timestamp(tout, "Interpolate gradients from L2 -> H1",
                               timing::duration(begin_interp, end_interp), 1);
        }
    }

    // We are done with the gradients/projections, so we can free the device
    // side vertices, as well as the element-to-DoF tables for H1 and L2.
    {
        Memory<int> h1_element_to_dof_I_mem = h1_element_to_dof.GetIMemory();
        h1_element_to_dof_I_mem.DeleteDevice(/* copy_back */ false);

        Memory<int> h1_element_to_dof_J_mem = h1_element_to_dof.GetJMemory();
        h1_element_to_dof_J_mem.DeleteDevice(/* copy_back */ false);

        Memory<int> l2_element_to_dof_I_mem = l2_element_to_dof.GetIMemory();
        l2_element_to_dof_I_mem.DeleteDevice(/* copy_back */ false);

        Memory<int> l2_element_to_dof_J_mem = l2_element_to_dof.GetJMemory();
        l2_element_to_dof_J_mem.DeleteDevice(/* copy_back */ false);

        local_vertices.DeleteDevice(/* copy_back */ false);
    }

    if (opts.verbose >= 2 && rank == 0) {
        logging::marker(tout, "[LDRB]: Compute fibers", 1);
    }
    struct timespec begin_ldrb, end_ldrb;
    timing::tick(&begin_ldrb);
    timing::tick(&t0);

    // Get read-only pointers to the internal arrays of the laplacian solutions
    const double *phi_epi = x_phi_epi->Read();
    const double *phi_lv  = x_phi_lv->Read();
    const double *phi_rv  = x_phi_rv->Read();

    // Get read-only pointers to the internal arrays of the gradients
    const double *grad_phi_epi_vals = grad_phi_epi->Read();
    const double *grad_phi_lv_vals  = grad_phi_lv->Read();
    const double *grad_phi_rv_vals  = grad_phi_rv->Read();
    const double *grad_psi_ab_vals  = grad_psi_ab->Read();

    // Setup GridFunctions to store the fibre directions in

    ParFiniteElementSpace *fiber_space = grad_phi_epi->ParFESpace();
    ParGridFunction F(fiber_space);
    ParGridFunction S(fiber_space);
    ParGridFunction T(fiber_space);

    // Get write-only pointers to the internal arrays
    double *F_vals = F.Write();
    double *S_vals = S.Write();
    double *T_vals = T.Write();

    timing::tick(&t1);
    if (opts.verbose >= 2 && rank == 0) {
        logging::timestamp(tout, "[Fiber]: Setup", timing::duration(t0, t1), 2);
    }

    // Calculate the fiber orientation
    timing::tick(&t0);
    tracing::roctx_range_push("define_fibers");
    define_fibers(x_phi_epi->Size(),
                  phi_epi, phi_lv, phi_rv,
                  grad_phi_epi_vals, grad_phi_lv_vals, grad_phi_rv_vals, grad_psi_ab_vals,
                  opts.alpha_endo, opts.alpha_epi, opts.beta_endo, opts.beta_epi,
                  F_vals, S_vals, T_vals
    );

    tracing::roctx_range_pop();
    timing::tick(&t1, /* barrier */ true, /* device barrier */ true);
    timing::tick(&end_ldrb, /* barrier */ false);
    if (opts.verbose >= 2 && rank == 0) {
        logging::timestamp(tout, "[Fiber]: Run define_fibers", timing::duration(t0, t1), 2);
        logging::timestamp(tout, "[Fiber]: Total", timing::duration(begin_ldrb, end_ldrb), 2, '=');
    } else if (opts.verbose && rank == 0) {
        logging::timestamp(tout, "Run LDRB", timing::duration(begin_ldrb, end_ldrb), 1, '-');
    }

    timing::tick(&end_fiber);
    if (opts.verbose >= 2 && rank == 0) {
        logging::timestamp(tout, "[LDRB]: Total", timing::duration(begin_fiber, end_fiber), 1, '=');
    } else if (opts.verbose && rank == 0) {
        logging::timestamp(tout, "Total", timing::duration(begin_fiber, end_fiber), 1, '=');
    }


    // Calculate the angle of the fibers. Useful for visualization.
    ParGridFunction angle(x_phi_epi->ParFESpace());
    double *angle_vals = angle.Write();
    const double *F_read_vals = F.Read();
    MFEM_FORALL(i, angle.Size(), { angle_vals[i] = asin(-F_read_vals[3*i])*180.0 / PI; });
    angle.HostRead();

    // Make sure the fiber directions are read back to the host before saving
    // the solutions.
    F.HostRead();
    S.HostRead();
    T.HostRead();

    if (opts.save_paraview || opts.save_mfem) {
        timing::tick(&t0);
        // Save the mesh and solutions
        // Set the basename of the mesh
        std::string mesh_basename = fs::remove_extension(
                fs::basename(std::string(opts.mesh_file)));

        if (opts.save_mfem) {
            // Output the normal solutions in the mfem subdirectory
            std::string mfem_output_dir(opts.output_dir);
            mfem_output_dir += "/mfem";
            fs::mksubdir(mfem_output_dir);

            // Save the MFEM mesh
            save::save_mesh(pmesh, mfem_output_dir, mesh_basename, rank);

            if (opts.save_laplacians) {
                save::save_solution(x_phi_epi, mfem_output_dir, mesh_basename, "_phi_epi.gf", rank);
                save::save_solution(x_phi_lv,  mfem_output_dir, mesh_basename, "_phi_lv.gf", rank);
                save::save_solution(x_phi_rv,  mfem_output_dir, mesh_basename, "_phi_rv.gf", rank);
            }

            if (opts.save_gradients) {
                save::save_solution(grad_phi_epi, mfem_output_dir, mesh_basename, "_grad_phi_epi.gf", rank);
                save::save_solution(grad_phi_lv,  mfem_output_dir, mesh_basename, "_grad_phi_lv.gf", rank);
                save::save_solution(grad_phi_rv,  mfem_output_dir, mesh_basename, "_grad_phi_rv.gf", rank);
                save::save_solution(grad_psi_ab,  mfem_output_dir, mesh_basename, "_grad_psi_ab.gf", rank);
            }
            if (opts.save_partitioning && rank == 0) {
                std::ofstream partitioning_out;
                std::ostringstream partitioning_out_fname;
                partitioning_out_fname << mfem_output_dir << "/partitioning.txt";
                partitioning_out.open(partitioning_out_fname.str().c_str());
                if (partitioning_out.fail()) {
                    std::cerr << "Could not open output file '" << partitioning_out_fname.str().c_str() << "'" << std::endl;
                    exit(1);
                }
                for (int i = 0; i < global_num_elements; i++) {
                    partitioning_out << partitioning[i] << std::endl;
                }
                partitioning_out.close();
            }
            // Save the solutions
            save::save_solution(&F,  mfem_output_dir, mesh_basename, "_F.gf", rank);
            save::save_solution(&S,  mfem_output_dir, mesh_basename, "_S.gf", rank);
            save::save_solution(&T,  mfem_output_dir, mesh_basename, "_T.gf", rank);
            save::save_solution(&angle,  mfem_output_dir, mesh_basename, "_angle.gf", rank);
        }

        // Save in paraview as well
        ParaViewDataCollection *pd = NULL;
        if (opts.save_paraview) {
            std::string paraview_path(opts.output_dir);
            paraview_path += "/paraview";
            fs::mksubdir(paraview_path);

            pd = new ParaViewDataCollection(mesh_basename, pmesh);
            pd->SetPrefixPath(paraview_path);
            if (opts.save_gradients) {
                pd->RegisterField("grad phi epi", grad_phi_epi);
                pd->RegisterField("grad phi lv",  grad_phi_lv);
                pd->RegisterField("grad phi rv",  grad_phi_rv);
                pd->RegisterField("grad psi ab",  grad_psi_ab);
            }
            if (opts.save_laplacians) {
                pd->RegisterField("phi epi", x_phi_epi);
                pd->RegisterField("phi lv",  x_phi_lv);
                pd->RegisterField("phi rv",  x_phi_rv);
                pd->RegisterField("psi ab",  x_psi_ab);
            }
            pd->RegisterField("F", &F);
            pd->RegisterField("S", &S);
            pd->RegisterField("T", &T);
            pd->RegisterField("angle", &angle);
            pd->SetLevelsOfDetail(1);
            pd->SetDataFormat(VTKFormat::BINARY);
            pd->SetHighOrderOutput(false);
            pd->Save();
        }
        timing::tick(&t1);
        if (opts.verbose && rank == 0)
            logging::timestamp(tout, "Save", timing::duration(t0, t1));
    } else if (opts.verbose >= 3 && rank == 0) {
        logging::info(std::cout,
                "Both mfem and paraview output is disabled "
                "(--no-save-mfem and --no-save-paraview)");
    }

    timing::tick(&end);
    if (opts.verbose && rank == 0) {
        logging::timestamp(tout, "Total time", timing::duration(start, end), 0, '=');
    }
}

