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

#include <fstream>

#include "mfem.hpp"
#include "calculus.hpp"
#include "util.hpp"
#include "ldrb-gpu.hpp"

using namespace std;
using namespace mfem;

// Save a solution (in form of a ParGridFunction) to a file named "<dir>/<prefix><suffix>".
void save_solution(
    ParGridFunction *x,
    string const& dir,
    string const& base_name,
    string const& suffix)
{
    string filename(dir);
    filename += "/";
    filename += base_name;
    filename += suffix;
    ofstream x_ofs(filename.c_str());
    x_ofs.precision(8);
    x->Save(x_ofs);
}


void laplace(
    ParGridFunction *x,
    ParMesh *pmesh,
    Array<int> &essential_boundaries,
    Array<int> &nonzero_essential_boundaries,
    Array<int> &zero_essential_boundaries,
    int apex,
    int dim,
    Options *opts,
    int rank)
{
    // Define a finite element space on the mesh.
    FiniteElementCollection *fec;

    fec = new H1_FECollection(1, dim);
    ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
    HYPRE_BigInt size = fespace->GlobalTrueVSize();
    if (opts->verbose > 1 && rank == 0)
        cout << "Number of finite element unknowns: " << size << endl;

    // Determine the list of true (i.e. parallel conforming) essential boundary
    // dofs, defined by the boundary attributes marked as essential (Dirichlet)
    // and converting to a list of true dofs..
    Array<int> ess_tdof_list;
    MFEM_ASSERT(pmesh->bdr_attributes.Size() != 0, "Boundary size cannot be zero.");

    fespace->GetEssentialTrueDofs(essential_boundaries, ess_tdof_list);

    // Set up the parallel linear form b(.) which corresponds to the right-hand
    // side of the FEM linear system, which in this case is (1,phi_i) where
    // phi_i are the basis functions in fespace.
    ParLinearForm b(fespace);
    ConstantCoefficient zero(0.0);
    b.AddDomainIntegrator(new DomainLFIntegrator(zero));
    b.Assemble();

    // Define the solution vector x as a parallel finite element grid function
    // corresponding to fespace. Initialize x with initial guess of zero, which
    // satisfies the boundary conditions.
    x->SetSpace(fespace);
    *x = 0.0;

    // Project the constant value 1.0 to all the essential boundaries marked as nonzero.
    ConstantCoefficient nonzero_bdr(1.0);
    x->ProjectBdrCoefficient(nonzero_bdr, nonzero_essential_boundaries);

    // Project the constant value 0.0 to all the essential boundaries marked as zero.
    ConstantCoefficient zero_bdr(0.0);
    x->ProjectBdrCoefficient(zero_bdr, zero_essential_boundaries);

    // For the laplacians involving the apex we need to treat the boundary
    // conditions a little different, as it is to be enforced on a single node
    // rather than a whole boundary surface. To do so, we make sure that the
    // apex is en essential true dof, and then we project the wanted value, in
    // this case 1.0, to only that node.
    if (apex >= 0) {
        // Initialize the internal data needed in the finite element space
        x->FESpace()->BuildDofToArrays();

        // Make sure the apex is in the list of essential true Dofs
        ess_tdof_list.Append(apex);

        Array<int> node_disp(1);
        node_disp[0] = apex;
        Vector node_disp_value(1);
        node_disp_value[0] = 1.0;

        VectorConstantCoefficient node_disp_value_coeff(node_disp_value);

        x->ProjectCoefficient(node_disp_value_coeff, node_disp);
    }

    // Set up the parallel bilinear form a(.,.) on the finite element space
    // corresponding to the Laplacian operator -Delta, by adding the
    // Diffusion domain integrator.
    ParBilinearForm a(fespace);
    ConstantCoefficient one(1.0);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));

    // Assemble the parallel bilinear form and the corresponding linear system.
    a.Assemble();

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, *x, b, A, X, B);

    // Solve the linear system A X = B.
    // Use the BoomerAMG preconditioner from hypre.
    HypreBoomerAMG *prec = NULL;
    prec = new HypreBoomerAMG;
    prec->SetPrintLevel(opts->verbose > 1 ? 1 : 0);
    CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(1e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(opts->verbose > 1 ? 1 : 0);
    if (prec) { cg.SetPreconditioner(*prec); }
    cg.SetOperator(*A);
    cg.Mult(B, X);
    delete prec;

    // Recover the parallel grid function corresponding to X. This is the local
    // finite element solution on each processor.
    a.RecoverFEMSolution(X, b, *x);
}


void laplace_phi_epi(
    ParGridFunction *x,
    ParMesh *mesh,
    int dim,
    Options *opts,
    int rank)
{
    // Define the following three arrays to determine (1) which boundary surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said boundary surfaces
    // should be set to zero.
    int nattr = mesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

    // Solve the Laplace equation from EPI (1.0) to (LV_ENDO union RV_ENDO) (0.0)
    essential_boundaries = 0;
    essential_boundaries[EPI    -1] = 1;
    essential_boundaries[LV_ENDO-1] = 1;
    essential_boundaries[RV_ENDO-1] = 1;

    nonzero_essential_boundaries = 0;
    nonzero_essential_boundaries[EPI-1] = 1;

    zero_essential_boundaries = 0;
    zero_essential_boundaries[LV_ENDO-1] = 1;
    zero_essential_boundaries[RV_ENDO-1] = 1;

    laplace(x, mesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, -1, dim, opts, rank);
}

void laplace_phi_lv(
    ParGridFunction *x,
    ParMesh *pmesh,
    int dim,
    Options *opts,
    int rank)
{
    // Define the following three arrays to determine (1) which border surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said border surfaces
    // should be set to zero.
    int nattr = pmesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

    // Solve the Laplace equation from LV_ENDO (1.0) to (RV_ENDO union EPI) (0.0)
    essential_boundaries = 0;
    essential_boundaries[EPI    -1] = 1;
    essential_boundaries[LV_ENDO-1] = 1;
    essential_boundaries[RV_ENDO-1] = 1;

    nonzero_essential_boundaries = 0;
    nonzero_essential_boundaries[LV_ENDO-1] = 1;

    zero_essential_boundaries = 0;
    zero_essential_boundaries[EPI    -1] = 1;
    zero_essential_boundaries[RV_ENDO-1] = 1;

    laplace(x, pmesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, -1, dim, opts, rank);
}

void laplace_phi_rv(
    ParGridFunction *x,
    ParMesh *pmesh,
    int dim,
    Options *opts,
    int rank)
{
    // Define the following three arrays to determine (1) which border surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said border surfaces
    // should be set to zero.
    int nattr = pmesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

    // Solve the Laplace equation from RV_ENDO (1.0) to (LV_ENDO union EPI) (0.0)
    essential_boundaries = 0;
    essential_boundaries[EPI    -1] = 1;
    essential_boundaries[LV_ENDO-1] = 1;
    essential_boundaries[RV_ENDO-1] = 1;

    nonzero_essential_boundaries = 0;
    nonzero_essential_boundaries[RV_ENDO-1] = 1;

    zero_essential_boundaries = 0;
    zero_essential_boundaries[EPI    -1] = 1;
    zero_essential_boundaries[LV_ENDO-1] = 1;

    laplace(x, pmesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, -1, dim, opts, rank);
}

void laplace_psi_ab(
    ParGridFunction *x,
    ParMesh *pmesh,
    int apex,
    int dim,
    Options *opts,
    int rank)
{
    // Define the following three arrays to determine (1) which border surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said border surfaces
    // should be set to zero.
    int nattr = pmesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

    // Find the apex by solving a laplacian with base solution = 0. The apex
    // will be at he maximum of the solution.
    essential_boundaries = 0;
    essential_boundaries[BASE-1] = 1;
    nonzero_essential_boundaries = 0;
    zero_essential_boundaries = 0;
    zero_essential_boundaries[BASE-1] = 1;
    laplace(x, pmesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, apex, dim, opts, rank);
}

// Find the vertex closest to the prescribed apex, Euclidean distance.
int find_apex_vertex(Mesh *mesh, Vector& apex)
{
    int apex_vertex = 0;
    double distance = numeric_limits<double>::max();
    for (int i = 0; i < mesh->GetNV(); i++) {
        double *vertices = mesh->GetVertex(i);
        double this_distance = (vertices[0]-apex[0]) * (vertices[0]-apex[0])
                             + (vertices[1]-apex[1]) * (vertices[1]-apex[1])
                             + (vertices[2]-apex[2]) * (vertices[2]-apex[2]);
        if (this_distance < distance) {
            apex_vertex = i;
            distance = this_distance;
        }
    }
    return apex_vertex;
}

int main(int argc, char *argv[])
{
    // 1. Initialize MPI and HYPRE
    Mpi::Init();
    int rank = Mpi::WorldRank();
    Hypre::Init();

    Options opts;

    // Set program defaults
    opts.mesh_file = NULL;
    opts.output_dir = "./out";
    opts.device_config = "hip";
    opts.verbose = 0;
    opts.prescribed_apex = Vector(3);
    opts.paraview = false;

    // Parse command-line options
    OptionsParser args(argc, argv);
    args.AddOption(&opts.verbose,         "-v", "--verbose", "Be verbose");
    args.AddOption(&opts.mesh_file,       "-m", "--mesh",    "Mesh file to use", true);
    args.AddOption(&opts.prescribed_apex, "-a", "--apex",    "Coordinate of apex, space separated list: 'x y z'.", true);
    args.AddOption(&opts.output_dir,      "-o", "--out",     "Directory for output files.");
    args.AddOption(&opts.device_config,   "-d", "--device",  "Device configuration string, see Device::Configure().");
    args.AddOption(&opts.paraview,        "-p", "--paraview", "-np", "--no-paraview", "Save data files for ParaView (paraview.org) visualization.");

    args.Parse();

    if (!args.Good()) {
        if (rank == 0)
            args.PrintUsage(cout);
        exit(1);
    }

    if (rank == 0 && opts.verbose > 1)
        args.PrintOptions(cout);

    // Set the basename of the mesh
    opts.mesh_basename = remove_extension(basename(std::string(opts.mesh_file)));

    // Make sure the output direcory exists
    mksubdir(opts.output_dir);

    struct timespec t0, t1;

    // 3. Enable hardware devices such as GPUs, and programming models such as
    //    HIP, CUDA, OCCA, RAJA and OpenMP based on command line options.
    Device device(opts.device_config);
    if (opts.verbose > 1 && rank == 0)
        device.Print();

    // 4. Load the mesh
    clock_gettime(CLOCK_MONOTONIC, &t0);
    Mesh mesh(opts.mesh_file, 1, 1);
    clock_gettime(CLOCK_MONOTONIC, &t1);

    if (opts.verbose > 1 && rank == 0) {
        cout << "Loaded meshfile '" << opts.mesh_file << "' "
             << "consisting of " << mesh.GetNV() << " vertices "
             << "and " << mesh.GetNE() << " elements" << endl;

    }
    if (opts.verbose && rank == 0)
        log_timing(cout, "Mesh load", timespec_duration(t0, t1));

    int dim = mesh.Dimension();

    // Set the apex node based on the prescribed apex coordinate.
    int apex = 0;
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);
        apex = find_apex_vertex(&mesh, opts.prescribed_apex);
        clock_gettime(CLOCK_MONOTONIC, &t1);

        if (opts.verbose > 1 && rank == 0) {
            double *closest_vertex = mesh.GetVertex(apex);
            cout << setprecision(2)
                 << "Found closest vertex to prescribed apex "
                 << "("  << opts.prescribed_apex[0]
                 << ", " << opts.prescribed_apex[1]
                 << ", " << opts.prescribed_apex[2]
                 << ") at "
                 << "("  << closest_vertex[0]
                 << ", " << closest_vertex[1]
                 << ", " << closest_vertex[2]
                 << ")." << endl;

        }
        if (opts.verbose && rank == 0)
            log_timing(cout, "Find apex", timespec_duration(t0, t1));
    }

    // 6. Define a parallel mesh by a partitioning of the serial mesh.
    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();

    // Solve the Laplace equation from EPI (1.0) to (LV_ENDO union RV_ENDO) (0.0)
    clock_gettime(CLOCK_MONOTONIC, &t0);
    ParGridFunction x_phi_epi;
    laplace_phi_epi(&x_phi_epi, &pmesh, dim, &opts, rank);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (opts.verbose && rank == 0)
        log_timing(cout, "phi_epi", timespec_duration(t0, t1));


    // Solve the Laplace equation from LV_ENDO (1.0) to (RV_ENDO union EPI) (0.0)
    clock_gettime(CLOCK_MONOTONIC, &t0);
    ParGridFunction x_phi_lv;
    laplace_phi_lv(&x_phi_lv, &pmesh, dim, &opts, rank);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (opts.verbose && rank == 0)
        log_timing(cout, "phi_lv", timespec_duration(t0, t1));


    // Solve the Laplace equation from RV_ENDO (1.0) to (LV_ENDO union EPI) (0.0)
    clock_gettime(CLOCK_MONOTONIC, &t0);
    ParGridFunction x_phi_rv;
    laplace_phi_rv(&x_phi_rv, &pmesh, dim, &opts, rank);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (opts.verbose && rank == 0)
        log_timing(cout, "phi_rv", timespec_duration(t0, t1));


    // Solve the Laplace equation from BASE (1.0) to APEX (0.0)
    clock_gettime(CLOCK_MONOTONIC, &t0);
    ParGridFunction x_psi_ab;
    laplace_psi_ab(&x_psi_ab, &pmesh, apex, dim, &opts, rank);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (opts.verbose && rank == 0)
        log_timing(cout, "psi_ab", timespec_duration(t0, t1));

    // Save the mesh
    {

        // Output the normal solutions in the mfem subdirectory
        string mfem_output_dir(opts.output_dir);
        mfem_output_dir += "/mfem";
        mksubdir(mfem_output_dir);

        // Save the MFEM mesh
        string mesh_out(mfem_output_dir);
        mesh_out += "/" + opts.mesh_basename + ".mesh";
        ofstream mesh_ofs(mesh_out.c_str());
        mesh_ofs.precision(8);
        mesh.Print(mesh_ofs);

        // Save the solutions
        save_solution(&x_phi_epi, mfem_output_dir, opts.mesh_basename, "_phi_epi.gf");
        save_solution(&x_phi_lv,  mfem_output_dir, opts.mesh_basename, "_phi_lv.gf");
        save_solution(&x_phi_rv,  mfem_output_dir, opts.mesh_basename, "_phi_rv.gf");
        save_solution(&x_psi_ab,  mfem_output_dir, opts.mesh_basename, "_psi_ab.gf");

    }
}

