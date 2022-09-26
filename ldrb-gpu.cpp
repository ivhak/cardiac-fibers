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
// along with hyprep.  If not, see <https://www.gnu.org/licenses/>.
//
// Authors: Iver Håkonsen <hakonseniver@yahoo.no

#include "mfem.hpp"
#include <fstream>

using namespace std;
using namespace mfem;


typedef struct {
    const char *mesh_file;
    const char *out_mesh;
    char *out_sol;
    const char *device_config;
    int order;
} options_t;

typedef enum {
    BASE   = 1,
    EPI    = 2,
    LV_ENDO = 3,
    RV_ENDO = 4
} mesh_attributes_t;

void laplace(
        ParMesh *pmesh,
        Array<int> &essential_borders,
        Array<int> &nonzero_essential_borders,
        Array<int> &zero_essential_borders,
        int dim,
        const char *filename,
        options_t *opts,
        int rank)
{
    // Define a finite element space on the mesh.
    FiniteElementCollection *fec;

    fec = new H1_FECollection(opts->order, dim);
    ParFiniteElementSpace fespace(pmesh, fec);
    HYPRE_BigInt size = fespace.GlobalTrueVSize();
    if (rank == 0)
        cout << "Number of finite element unknowns: " << size << endl;

    // Determine the list of true (i.e. parallel conforming) essential boundary
    // dofs, defined by the boundary attributes marked as essential (Dirichlet)
    // and converting to a list of true dofs..
    Array<int> ess_tdof_list;
    MFEM_ASSERT(pmesh->bdr_attributes.Size() != 0, "Boundary size cannot be zero.");

    fespace.GetEssentialTrueDofs(essential_borders, ess_tdof_list);

    // Set up the parallel linear form b(.) which corresponds to the right-hand
    // side of the FEM linear system, which in this case is (1,phi_i) where
    // phi_i are the basis functions in fespace.
    ParLinearForm b(&fespace);
    ConstantCoefficient one(1.0);
    b.AddDomainIntegrator(new DomainLFIntegrator(one));
    b.Assemble();

    // Define the solution vector x as a parallel finite element grid function
    // corresponding to fespace. Initialize x with initial guess of zero, which
    // satisfies the boundary conditions.
    ParGridFunction x(&fespace);
    x = 0.0;

    // Project the constant value 1.0 to all the essential borders marked as nonzero.
    ConstantCoefficient nonzero_bdr(1.0);
    x.ProjectBdrCoefficient(nonzero_bdr, nonzero_essential_borders);

    // Project the constant value 0.0 to all the essential borders marked as zero.
    ConstantCoefficient zero_bdr(0.0);
    x.ProjectBdrCoefficient(zero_bdr, zero_essential_borders);

    // Set up the parallel bilinear form a(.,.) on the finite element space
    // corresponding to the Laplacian operator -Delta, by adding the
    // Diffusion domain integrator.
    ParBilinearForm a(&fespace);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));

    // Assemble the parallel bilinear form and the corresponding linear system.
    a.Assemble();

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

    // Solve the linear system A X = B.
    // * With full assembly, use the BoomerAMG preconditioner from hypre.
    Solver *prec = NULL;
    prec = new HypreBoomerAMG;
    CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(1e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(1);
    if (prec) { cg.SetPreconditioner(*prec); }
    cg.SetOperator(*A);
    cg.Mult(B, X);
    delete prec;

    // Recover the parallel grid function corresponding to X. This is the local
    // finite element solution on each processor.
    a.RecoverFEMSolution(X, b, x);

    ofstream sol(filename);
    sol.precision(8);
    x.Save(sol);
}

void laplace_psi_epi(ParMesh *pmesh, int dim, options_t *opts, int rank)
{
    // Define the following three arrays to determine (1) which border surfaces
    // to include in the Laplace equation, (2) which of said borders should be
    // set to a nonzero value (1.0) and (3) which of said border surfaces
    // should be set to zero.
    int nattr = pmesh->bdr_attributes.Max();
    Array<int> essential_borders(nattr);
    Array<int> nonzero_essential_borders(nattr);
    Array<int> zero_essential_borders(nattr);

    // Laplace PSI_EPI:
    // Solve the Laplace equation from EPI (1.0) to (LV_ENDO union RV_ENDO) (0.0)
    essential_borders[BASE   -1] = 0;
    essential_borders[EPI    -1] = 1;
    essential_borders[LV_ENDO-1] = 1;
    essential_borders[RV_ENDO-1] = 1;

    nonzero_essential_borders[EPI    -1] = 1;
    nonzero_essential_borders[LV_ENDO-1] = 0;
    nonzero_essential_borders[RV_ENDO-1] = 0;

    zero_essential_borders[EPI    -1] = 0;
    zero_essential_borders[LV_ENDO-1] = 1;
    zero_essential_borders[RV_ENDO-1] = 1;

    laplace(pmesh, essential_borders, nonzero_essential_borders, zero_essential_borders, dim, "out/heart01_psi_epi.gf", opts, rank);
}

void laplace_psi_lv(ParMesh *pmesh, int dim, options_t *opts, int rank)
{
    // Define the following three arrays to determine (1) which border surfaces
    // to include in the Laplace equation, (2) which of said borders should be
    // set to a nonzero value (1.0) and (3) which of said border surfaces
    // should be set to zero.
    int nattr = pmesh->bdr_attributes.Max();
    Array<int> essential_borders(nattr);
    Array<int> nonzero_essential_borders(nattr);
    Array<int> zero_essential_borders(nattr);

    // Laplace PSI_LV;
    // Solve the Laplace equation from LV_ENDO (1.0) to (RV_ENDO union EPI) (0.0)
    essential_borders[BASE   -1] = 0;
    essential_borders[EPI    -1] = 1;
    essential_borders[LV_ENDO-1] = 1;
    essential_borders[RV_ENDO-1] = 1;

    nonzero_essential_borders[EPI    -1] = 0;
    nonzero_essential_borders[LV_ENDO-1] = 1;
    nonzero_essential_borders[RV_ENDO-1] = 0;

    zero_essential_borders[EPI    -1] = 1;
    zero_essential_borders[LV_ENDO-1] = 0;
    zero_essential_borders[RV_ENDO-1] = 1;

    laplace(pmesh, essential_borders, nonzero_essential_borders, zero_essential_borders, dim,  "out/heart01_psi_lv.gf",opts, rank);
}

void laplace_psi_rv(ParMesh *pmesh, int dim, options_t *opts, int rank)
{
    // Define the following three arrays to determine (1) which border surfaces
    // to include in the Laplace equation, (2) which of said borders should be
    // set to a nonzero value (1.0) and (3) which of said border surfaces
    // should be set to zero.
    int nattr = pmesh->bdr_attributes.Max();
    Array<int> essential_borders(nattr);
    Array<int> nonzero_essential_borders(nattr);
    Array<int> zero_essential_borders(nattr);

    // Laplace PSI_RV;
    // Solve the Laplace equation from RV_ENDO (1.0) to (LV_ENDO union EPI) (0.0)
    essential_borders[BASE   -1] = 0;
    essential_borders[EPI    -1] = 1;
    essential_borders[LV_ENDO-1] = 1;
    essential_borders[RV_ENDO-1] = 1;

    nonzero_essential_borders[EPI    -1] = 0;
    nonzero_essential_borders[LV_ENDO-1] = 0;
    nonzero_essential_borders[RV_ENDO-1] = 1;

    zero_essential_borders[EPI    -1] = 1;
    zero_essential_borders[LV_ENDO-1] = 1;
    zero_essential_borders[RV_ENDO-1] = 0;

    laplace(pmesh, essential_borders, nonzero_essential_borders, zero_essential_borders, dim, "out/heart01_psi_rv.gf", opts, rank);
}

int main(int argc, char *argv[])
{
    // 1. Initialize MPI and HYPRE
    Mpi::Init();
    int rank = Mpi::WorldRank();
    Hypre::Init();

    options_t opts;
    // 2. Parse command-line options

    // Set program defaults
    opts.mesh_file = NULL;
    opts.out_mesh = "./out.mesh";
    opts.out_sol = NULL;
    opts.device_config = "cpu";
    opts.order = 1;

    OptionsParser args(argc, argv);
    args.AddOption(&opts.mesh_file, "-m", "--mesh",
            "Mesh file to use (required)");
    args.AddOption(&opts.out_mesh, "-o", "--out",
            "Output file");
    args.AddOption(&opts.device_config, "-d", "--device",
            "Device configuration string, see Device::Configure().");
    args.AddOption(&opts.order, "-o", "--order",
            "Finite element order (polynomial degree) or -1 for "
            "isoparametric space.");
    args.Parse();

    // Set the solution file to be the same as the output mesh but with the
    // ".gf" suffix
    {
        // TODO: Fix this hacky thing
        opts.out_sol = (char *)malloc(strlen(opts.out_mesh) * sizeof(char));
        strncpy(opts.out_sol, opts.out_mesh, strlen(opts.out_mesh));
        char *suff = strstr(opts.out_sol, ".mesh");
        strncpy(suff, ".gf", strlen(".gf"));
        suff[strlen(".gf")] = 0;
    }

    if (!args.Good() || opts.mesh_file == NULL) {
        if (rank == 0)
            args.PrintUsage(cout);
        exit(1);
    }

    if (rank == 0)
        args.PrintOptions(cout);

    // 3. Enable hardware devices such as GPUs, and programming models such as
    //    HIP, CUDA, OCCA, RAJA and OpenMP based on command line options.
    Device device(opts.device_config);
    if (rank == 0)
        device.Print();

    // 4. Load the mesh
    Mesh mesh(opts.mesh_file, 1, 1);
    int dim = mesh.Dimension();

#ifdef DEBUG
    {
        // Check that the border elements have to proper attributes.
        //                   Base Epi  Lv   Rv
        int attributes[4] = {0,   0,   0,   0};

        int num_border_elements = mesh.GetNBE();
        for (int i = 0; i < num_border_elements; i++) {
            Element *ele = mesh.GetBdrElement(i);
            int attr = ele->GetAttribute();
            MFEM_ASSERT(attr >= BASE && attr <=RV_ENDO,
                    "The element attributes should be 1 (base), 2 (epicardium), "
                    "3 (left ventricle endocardium) or 4 (right ventricle endocardium)");
            attributes[attr-1]++;
        }

        if (rank == 0) {
            cout << "Number of border elements: " << num_border_elements << endl
                 << "Number of border elements at each surface:" << endl
                 << "\tBase: " << attributes[0] << " elements" << endl
                 << "\tEpi:  " << attributes[1] << " elements" << endl
                 << "\tLv: "   << attributes[2] << " elements" << endl
                 << "\tRv: "   << attributes[3] << " elements" << endl;
        }
    }
#endif

    // 6. Define a parallel mesh by a partitioning of the serial mesh.
    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();

    // Laplace PSI_EPI:
    // Solve the Laplace equation from EPI (1.0) to (LV_ENDO union RV_ENDO) (0.0)
    laplace_psi_epi(&pmesh, dim,&opts, rank);

    // Laplace PSI_LV;
    // Solve the Laplace equation from LV_ENDO (1.0) to (RV_ENDO union EPI) (0.0)
    laplace_psi_lv(&pmesh, dim, &opts, rank);

    // Laplace PSI_RV
    // Solve the Laplace equation from RV_ENDO (1.0) to (LV_ENDO union EPI) (0.0)
    laplace_psi_rv(&pmesh, dim, &opts, rank);

    // TODO: Laplace GAMMA_AB (apex -> base)
    // TODO: Solve the Laplace equation from BASE (1.0) to APEX (0.0)

    // Save the mesh
    {
        ofstream mesh_ofs(opts.out_mesh);
        mesh_ofs.precision(8);
        pmesh.Print(mesh_ofs);
    }
}

