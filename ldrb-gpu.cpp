// Copyright (C) 2022 Iver Håkonsen
//
// ldrb-gpu is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ldrb-gpu is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with ldrb-gpu.  If not, see <https://www.gnu.org/licenses/>.
//
// Authors: Iver Håkonsen <hakonseniver@yahoo.no

#include <iostream>
#include <fstream>

#include "mfem.hpp"
#include "mfem/general/forall.hpp"
#include "util.hpp"
#include "ldrb-gpu.hpp"
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
    int rank)
{
    struct timespec t0, t1;


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
        logging::timestamp(tout, "Setup boundary cond.", timing::duration(t0, t1), 2);
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
        logging::timestamp(tout, "Assemble RHS", timing::duration(t0, t1), 2);
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
        logging::timestamp(tout, "Assemble LHS", timing::duration(t0, t1), 2);
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
        logging::timestamp(tout, "Form linear system", timing::duration(t0, t1), 2);
    }
    tracing::roctx_range_pop(); // Form linear system

    // Solve the linear system A X = B.
    tracing::roctx_range_push("Solve");
    timing::tick(&t0);

    solver->SetOperator(*A);
    solver->Mult(B, X);
    timing::tick(&t1);
    if (verbose >= 2 && rank == 0) {
        logging::timestamp(tout, "Solve", timing::duration(t0, t1), 2);
    }
    tracing::roctx_range_pop(); // Solve

    // Recover the parallel grid function corresponding to X. This is the local
    // finite element solution on each processor.
    a.RecoverFEMSolution(X, b, *x);
}

int main(int argc, char *argv[])
{
    // 1. Initialize MPI and HYPRE
    Mpi::Init();
    int rank = Mpi::WorldRank();
    int nranks = Mpi::WorldSize();
    Hypre::Init();

    Options opts;

    // Set program defaults
    opts.mesh_file = NULL;
    opts.output_dir = "./out";
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

    opts.save_mfem = true;
    opts.save_paraview = false;

    opts.itol = 1e-1;

    opts.alpha_endo =  60.0;
    opts.alpha_epi  = -60.0;
    opts.beta_endo  = -60.0;
    opts.beta_epi   =  60.0;

    opts.base_id = 1;
    opts.epi_id  = 2;
    opts.lv_id   = 3;
    opts.rv_id   = 4;

    opts.geom_has_rv = true;

    opts.solver = 0;

#if defined(MFEM_USE_HIP) || defined(MFEM_USE_CUDA)
    opts.gpu_tuned_amg = false;
#endif

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
            "Mesh file to use.\n"
            "\tSee https://mfem.org/mesh-formats/ for a list of suppored formats.", true);

    args.AddOption(&opts.prescribed_apex,
            "-a", "--apex",
            "Coordinate of apex, space separated list: 'x y z'.", true);

    args.AddOption(&opts.output_dir,
            "-o", "--out",
            "Directory for output files.");

    args.AddOption(&opts.time_to_file,
            "-t", "--time-to-file",
            "-nt","--no-time-to-file",
            "Output time log to <OUT>/time.txt rather than stdout,\n"
            "\twhere <OUT> is set with the -o (--out) flag.");

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

    args.AddOption(&opts.itol,
            "-it", "--interpolation-tolerance",
            "Tolerance for LDRB interpolations.");

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
            "\tSet to -1 if there is no right ventricle in the geometry,\n "
            "\te.g. for a single ventricle geometry.");

    args.AddOption(&opts.solver,
            "-s", "--solver",
            "Solver to use. Options are:\n"
            "\t    0: HyprePCG (See mfem::HyprePCG).\n"
            "\t    1: CGSolver (See mfem::CGSolver).");

#if defined (MFEM_USE_HIP) || defined (MFEM_USE_CUDA)
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

    if (opts.rv_id == -1) {
        opts.geom_has_rv = false;
    }

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
    Device device(opts.device_config);
    if (opts.verbose >= 3 && rank == 0) {
        device.Print();
    }

    // Load the mesh
    timing::tick(&t0);
    Mesh mesh(opts.mesh_file, 1, 1);
    timing::tick(&t1);

    if (opts.verbose >= 3 && rank == 0) {
        std::string msg = "Loaded meshfile '" + std::string(opts.mesh_file) + "' "
             + "consisting of " + std::to_string(mesh.GetNV()) + " vertices "
             + "and " + std::to_string(mesh.GetNE()) + " elements";
        logging::info(std::cout, msg);

    }
    if (opts.verbose >= 2 && rank == 0)
        logging::timestamp(tout, "Mesh load", timing::duration(t0, t1));

    // Define a parallel mesh by a partitioning of the serial mesh.
    timing::tick(&t0);
    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();
    timing::tick(&t1);
    if (opts.verbose && rank == 0)
        logging::timestamp(tout, "Partition mesh", timing::duration(t0, t1));


    // Each rank finds the vertex in its submesh that is closest to the
    // prescribed apex. The rank with the closest one sets up the boundary
    // conditions in that vertex.
    //
    // XXX (ivhak): What if the apex vertex is in multiple submeshes?
    int apex = -1;
    {
        timing::tick(&t0);
        std::vector<double> min_buffer(nranks, -1);
        double min_distance = std::numeric_limits<double>::max();
        for (int i = 0; i < pmesh.GetNV(); i++) {
            double *vertices = pmesh.GetVertex(i);
            double this_distance =
                (vertices[0]-opts.prescribed_apex[0]) * (vertices[0]-opts.prescribed_apex[0])
              + (vertices[1]-opts.prescribed_apex[1]) * (vertices[1]-opts.prescribed_apex[1])
              + (vertices[2]-opts.prescribed_apex[2]) * (vertices[2]-opts.prescribed_apex[2]);
            if (this_distance < min_distance) {
                apex = i;
                min_distance = this_distance;
            }
        }

        MPI_Allgather(&min_distance,     1, MPI_DOUBLE,
                      min_buffer.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD);
        auto result = std::min_element(min_buffer.begin(), min_buffer.end());
        auto target_rank = result - min_buffer.begin();
        if (target_rank != rank) {
            apex = -1;
        }

        timing::tick(&t1);
        if (opts.verbose && rank == 0)
            logging::timestamp(tout, "Find apex", timing::duration(t0, t1));
    }


    // TODO (ivhak): Test different solvers + preconditioners.
    //
    // Setup the linear solver.

    timing::tick(&t0);
    HypreBoomerAMG *prec = NULL;
    {
        prec = new HypreBoomerAMG;
        prec->SetPrintLevel(opts.verbose >= 4 ? 1 : 0);

#if defined(MFEM_USE_HIP) || defined(MFEM_USE_CUDA)
        // FIXME (ivhak): Gives wrong solution!
        if (opts.gpu_tuned_amg) {
            prec->SetCoarsening(8);           // PMIS
            prec->SetCycleNumSweeps(1,1);     // 1 sweep on the up and down cycle
            prec->SetInterpolation(17);       // extended+i, matrix-matrix
            prec->SetAggressiveCoarsening(1); // Number of levels of aggressive coarsening
            // prec->SetStrengthThresh(0.5);
            prec->SetRelaxType(7);            // weighted Jacobi
        }
#endif
    }

    Solver *solver;
    if (opts.solver == 0) {
        HyprePCG *cg = new HyprePCG(MPI_COMM_WORLD);
        cg->SetTol(1e-12);
        cg->SetMaxIter(2000);
        cg->SetPrintLevel(opts.verbose >= 4 ? 1 : 0);
        if (prec) cg->SetPreconditioner(*prec);
        solver = cg;
    } else if (opts.solver == 1) {
        CGSolver *cg = new CGSolver(MPI_COMM_WORLD);
        cg->SetAbsTol(1e-12);
        cg->SetMaxIter(2000);
        cg->SetPrintLevel(opts.verbose >= 4 ? 1 : 0);
        if (prec) cg->SetPreconditioner(*prec);
        solver = cg;
    } else {
        std::cerr << "Invalid solver: " << opts.solver << std::endl;
        exit(1);
    }


    timing::tick(&t1);
    if (opts.verbose && rank == 0) {
        logging::timestamp(tout, "Setup precond. and solver", timing::duration(t0,t1));
    }

    // Time the fiber orientation calculations
    struct timespec start_fiber, end_fiber;
    timing::tick(&start_fiber);

    if (opts.verbose && rank == 0) {
        logging::marker(tout, "Compute fiber orientation");
    }


    // Set up the finite element collection. We use  first order H1-conforming finite elements.
    H1_FECollection fec(1, pmesh.Dimension());

    // Set up two finite element spaces: one for scalar values (1D) , and one for 3d vectors
    ParFiniteElementSpace fespace_scalar_h1(&pmesh, &fec);
    ParFiniteElementSpace fespace_vector_h1(&pmesh, &fec, 3, mfem::Ordering::byVDIM);

    int nattr = pmesh.bdr_attributes.Max();
    Array<int> ess_bdr(nattr);          // Essential boundaries
    Array<int> nonzero_ess_bdr(nattr);  // Essential boundaries with value 1.0
    Array<int> zero_ess_bdr(nattr);     // Essential boundaries with value 0.0

    // Solve the Laplace equation from EPI (1.0) to (LV_ENDO union RV_ENDO) (0.0)
    ParGridFunction x_phi_epi(&fespace_scalar_h1);
    x_phi_epi.UseDevice(true);
    {
        if (opts.verbose > 1 && rank == 0) {
            logging::marker(tout, "Compute phi_epi", 1);
        }
        timing::tick(&t0);

        ess_bdr = 0;
        ess_bdr[opts.epi_id-1] = 1;
        ess_bdr[opts.lv_id -1] = 1;
        if (opts.geom_has_rv)
            ess_bdr[opts.rv_id-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.epi_id-1] = 1;

        zero_ess_bdr = 0;
        zero_ess_bdr[opts.lv_id-1] = 1;
        if (opts.geom_has_rv)
            zero_ess_bdr[opts.rv_id-1] = 1;

        tracing::roctx_range_push("laplace x_phi_epi");

        laplace(&x_phi_epi, solver,
                ess_bdr, nonzero_ess_bdr, zero_ess_bdr,
                -1, opts.verbose, tout, rank);

        tracing::roctx_range_pop();

        timing::tick(&t1);
        if (opts.verbose && rank == 0) {
            if (opts.verbose > 1) {
                logging::timestamp(tout, "Total", timing::duration(t0, t1), 2, '=');
            } else {
                logging::timestamp(tout, "Compute phi_epi", timing::duration(t0, t1), 1);
            }
        }
    }


    // Solve the Laplace equation from LV_ENDO (1.0) to (RV_ENDO union EPI) (0.0)
    ParGridFunction x_phi_lv(&fespace_scalar_h1);
    x_phi_lv.UseDevice(true);
    {
        if (opts.verbose > 1 && rank == 0) {
            logging::marker(tout, "Compute phi_lv", 1);
        }
        timing::tick(&t0);

        ess_bdr = 0;
        ess_bdr[opts.epi_id-1] = 1;
        ess_bdr[opts.lv_id -1] = 1;
        if (opts.geom_has_rv)
            ess_bdr[opts.rv_id-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.lv_id-1] = 1;

        zero_ess_bdr = 0;
        zero_ess_bdr[opts.epi_id-1] = 1;
        if (opts.geom_has_rv)
            zero_ess_bdr[opts.rv_id-1] = 1;

        tracing::roctx_range_push("laplace x_phi_lv");

        laplace(&x_phi_lv, solver,
                ess_bdr, nonzero_ess_bdr, zero_ess_bdr,
                -1, opts.verbose, tout, rank);

        tracing::roctx_range_pop();

        timing::tick(&t1);
        if (opts.verbose && rank == 0) {
            if (opts.verbose > 1) {
                logging::timestamp(tout, "Total", timing::duration(t0, t1), 2, '=');
            } else {
                logging::timestamp(tout, "Compute phi_lv", timing::duration(t0, t1), 1);
            }
        }
    }


    // Solve the Laplace equation from RV_ENDO (1.0) to (LV_ENDO union EPI) (0.0)
    ParGridFunction x_phi_rv(&fespace_scalar_h1);
    x_phi_rv.UseDevice(true);
    if (opts.geom_has_rv) {
        if (opts.verbose > 1 && rank == 0) {
            logging::marker(tout, "Compute phi_rv", 1);
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

        laplace(&x_phi_rv, solver,
                ess_bdr, nonzero_ess_bdr, zero_ess_bdr,
                -1, opts.verbose, tout, rank);

        tracing::roctx_range_pop();
        timing::tick(&t1);
        if (opts.verbose && rank == 0) {
            if (opts.verbose > 1) {
                logging::timestamp(tout, "Total", timing::duration(t0, t1), 2, '=');
            } else {
                logging::timestamp(tout, "Compute phi_rv", timing::duration(t0, t1), 1);
            }
        }
    } else {
        x_phi_rv = 0.0;
    }


    // Solve the Laplace equation from BASE (1.0) to APEX (0.0)
    ParGridFunction x_psi_ab(&fespace_scalar_h1);
    x_psi_ab.UseDevice(true);
    {
        if (opts.verbose > 1 && rank == 0) {
            logging::marker(tout, "Compute psi_ab", 1);
        }
        timing::tick(&t0);

        ess_bdr = 0;
        ess_bdr[opts.base_id-1] = 1;

        nonzero_ess_bdr = 0;
        nonzero_ess_bdr[opts.base_id-1] = 1;

        zero_ess_bdr = 0;

        tracing::roctx_range_push("laplace x_psi_ab");

        laplace(&x_psi_ab, solver,
                ess_bdr, nonzero_ess_bdr, zero_ess_bdr,
                apex, opts.verbose, tout, rank);

        tracing::roctx_range_pop();

        timing::tick(&t1);
        if (opts.verbose && rank == 0 ) {
            if (opts.verbose > 1) {
                logging::timestamp(tout, "Total", timing::duration(t0, t1), 2, '=');
            } else {
                logging::timestamp(tout, "Compute psi_epi", timing::duration(t0, t1), 1);
            }
        }
    }


    // Compute gradients for phi_epi, phi_lv, phi_rv and psi_ab
    ParGridFunction grad_phi_epi(&fespace_vector_h1);
    grad_phi_epi.UseDevice(true);
    {
        timing::tick(&t0);
        tracing::roctx_range_push("grad_phi_epi");

        GradientGridFunctionCoefficient ggfc(&x_phi_epi);
        grad_phi_epi.ProjectCoefficient(ggfc);

        tracing::roctx_range_pop();
        timing::tick(&t1);
        if (opts.verbose && rank == 0) {
            logging::timestamp(tout, "Compute grad_phi_epi", timing::duration(t0, t1), 1);
        }
    }

    ParGridFunction grad_phi_lv(&fespace_vector_h1);
    grad_phi_lv.UseDevice(true);
    {
        timing::tick(&t0);
        tracing::roctx_range_push("grad_phi_lv");

        GradientGridFunctionCoefficient ggfc(&x_phi_lv);
        grad_phi_lv.ProjectCoefficient(ggfc);

        tracing::roctx_range_pop();
        timing::tick(&t1);
        if (opts.verbose && rank == 0) {
            logging::timestamp(tout, "Compute grad_phi_lv", timing::duration(t0, t1), 1);
        }
    }

    ParGridFunction grad_phi_rv(&fespace_vector_h1);
    grad_phi_rv.UseDevice(true);
    if (opts.geom_has_rv) {
        timing::tick(&t0);
        tracing::roctx_range_push("grad_phi_rv");

        GradientGridFunctionCoefficient ggfc(&x_phi_rv);
        grad_phi_rv.ProjectCoefficient(ggfc);

        tracing::roctx_range_pop();
        timing::tick(&t1);
        if (opts.verbose && rank == 0) {
            logging::timestamp(tout, "Compute grad_phi_rv", timing::duration(t0, t1), 1);
        }
    } else {
        grad_phi_rv = 0.0;
    }

    ParGridFunction grad_psi_ab(&fespace_vector_h1);
    grad_psi_ab.UseDevice(true);
    {
        timing::tick(&t0);
        tracing::roctx_range_push("grad_psi_ab");

        GradientGridFunctionCoefficient ggfc(&x_psi_ab);
        grad_psi_ab.ProjectCoefficient(ggfc);

        tracing::roctx_range_pop();
        timing::tick(&t1);
        if (opts.verbose && rank == 0) {
            logging::timestamp(tout, "Compute grad_psi_ab", timing::duration(t0, t1), 1);
        }
    }

    // Get read-only pointers to the internal arrays of the laplacian solutions
    const double *phi_epi = x_phi_epi.Read();
    const double *phi_lv  = x_phi_lv.Read();
    const double *phi_rv  = x_phi_rv.Read();

    // Get read-only pointers to the internal arrays of the gradients
    const double *grad_phi_epi_vals = grad_phi_epi.Read();
    const double *grad_phi_lv_vals  = grad_phi_lv.Read();
    const double *grad_phi_rv_vals  = grad_phi_rv.Read();
    const double *grad_psi_ab_vals  = grad_psi_ab.Read();

    // Setup GridFunctions to store the fibre directions in
    ParGridFunction F(&fespace_vector_h1);
    ParGridFunction S(&fespace_vector_h1);
    ParGridFunction T(&fespace_vector_h1);

    // Get write-only pointers to the internal arrays
    double *F_vals = F.Write();
    double *S_vals = S.Write();
    double *T_vals = T.Write();


    // Calculate the fiber orientation
    timing::tick(&t0);
    util::tracing::roctx_range_push("define_fibers");
    define_fibers(
            pmesh.GetNV(),
            phi_epi, phi_lv, phi_rv,
            grad_phi_epi_vals, grad_phi_lv_vals, grad_phi_rv_vals, grad_psi_ab_vals,
            opts.alpha_endo, opts.alpha_epi, opts.beta_endo, opts.beta_epi,
            F_vals, S_vals, T_vals, opts.itol
    );

    util::tracing::roctx_range_pop();
    timing::tick(&t1);
    if (opts.verbose && rank == 0)
        logging::timestamp(tout, "Define fiber orientation", timing::duration(t0, t1), 1);

    timing::tick(&end_fiber);
    if (opts.verbose && rank == 0) {
        logging::timestamp(tout, "Total (fiber)", timing::duration(start_fiber, end_fiber), 1, '=');
    }

    // Make sure the fiber directions are read back to the host before saving
    // the solutions.
    F.HostRead();
    S.HostRead();
    T.HostRead();

    if (opts.save_paraview || opts.save_mfem) {
        timing::tick(&t0);
        // Save the mesh and solutions
        // Set the basename of the mesh
        opts.mesh_basename = fs::remove_extension(
                fs::basename(std::string(opts.mesh_file)));

        if (opts.save_mfem) {
            // Output the normal solutions in the mfem subdirectory
            std::string mfem_output_dir(opts.output_dir);
            mfem_output_dir += "/mfem";
            fs::mksubdir(mfem_output_dir);

            // Save the MFEM mesh
            save::save_mesh(&pmesh, mfem_output_dir, opts.mesh_basename, rank);

#ifdef DEBUG
            std::string debug_dir = mfem_output_dir + "/debug";
            fs::mksubdir(debug_dir);
            save::save_solution(&x_phi_epi, debug_dir, opts.mesh_basename, "_phi_epi.gf", rank);
            save::save_solution(&x_phi_lv,  debug_dir, opts.mesh_basename, "_phi_lv.gf", rank);
            save::save_solution(&x_phi_rv,  debug_dir, opts.mesh_basename, "_phi_rv.gf", rank);
            save::save_solution(&x_psi_ab,  debug_dir, opts.mesh_basename, "_psi_ab.gf", rank);

            save::save_solution(&grad_phi_epi, debug_dir, opts.mesh_basename, "_grad_phi_epi.gf", rank);
            save::save_solution(&grad_phi_lv,  debug_dir, opts.mesh_basename, "_grad_phi_lv.gf", rank);
            save::save_solution(&grad_phi_rv,  debug_dir, opts.mesh_basename, "_grad_phi_rv.gf", rank);
            save::save_solution(&grad_psi_ab,  debug_dir, opts.mesh_basename, "_grad_psi_ab.gf", rank);
#endif
            // Save the solutions
            save::save_solution(&F,  mfem_output_dir, opts.mesh_basename, "_F.gf", rank);
            save::save_solution(&S,  mfem_output_dir, opts.mesh_basename, "_S.gf", rank);
            save::save_solution(&T,  mfem_output_dir, opts.mesh_basename, "_T.gf", rank);
        }

        // Save in paraview as well
        ParaViewDataCollection *pd = NULL;
        if (opts.save_paraview) {
            std::string paraview_path(opts.output_dir);
            paraview_path += "/paraview";
            fs::mksubdir(paraview_path);

            pd = new ParaViewDataCollection(opts.mesh_basename, &pmesh);
            pd->SetPrefixPath(paraview_path);
#ifdef DEBUG
            pd->RegisterField("grad phi epi", &grad_phi_epi);
            pd->RegisterField("grad phi lv",  &grad_phi_lv);
            pd->RegisterField("grad phi rv",  &grad_phi_rv);
            pd->RegisterField("grad psi ab",  &grad_psi_ab);
            pd->RegisterField("phi epi", &x_phi_epi);
            pd->RegisterField("phi lv",  &x_phi_lv);
            pd->RegisterField("phi rv",  &x_phi_rv);
            pd->RegisterField("psi ab",  &x_psi_ab);
#endif
            pd->RegisterField("F", &F);
            pd->RegisterField("S", &S);
            pd->RegisterField("T", &T);
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

