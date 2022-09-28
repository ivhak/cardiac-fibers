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

#include "mfem.hpp"
#include <fstream>

using namespace std;
using namespace mfem;


typedef struct {
    const char *mesh_file;
    const char *out;
    const char *device_config;
    int order;
} Options;

typedef enum {
    BASE   = 1,
    EPI    = 2,
    LV_ENDO = 3,
    RV_ENDO = 4
} MeshAttributes;

void laplace(
        ParGridFunction *x,
        ParMesh *pmesh,
        Array<int> &essential_boundaries,
        Array<int> &nonzero_essential_boundaries,
        Array<int> &zero_essential_boundaries,
        int dim,
        Options *opts,
        int rank)
{
    // Define a finite element space on the mesh.
    FiniteElementCollection *fec;

    fec = new H1_FECollection(opts->order, dim);
    ParFiniteElementSpace *fespace = new ParFiniteElementSpace(pmesh, fec);
    HYPRE_BigInt size = fespace->GlobalTrueVSize();
    if (rank == 0)
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
    ConstantCoefficient one(1.0);
    b.AddDomainIntegrator(new DomainLFIntegrator(one));
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

    // Set up the parallel bilinear form a(.,.) on the finite element space
    // corresponding to the Laplacian operator -Delta, by adding the
    // Diffusion domain integrator.
    ParBilinearForm a(fespace);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));

    // Assemble the parallel bilinear form and the corresponding linear system.
    a.Assemble();

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, *x, b, A, X, B);

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
    a.RecoverFEMSolution(X, b, *x);
}

void laplace_psi_epi(ParGridFunction *x, ParMesh *pmesh, int dim, Options *opts, int rank)
{
    // Define the following three arrays to determine (1) which border surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said border surfaces
    // should be set to zero.
    int nattr = pmesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

    // Laplace PSI_EPI:
    // Solve the Laplace equation from EPI (1.0) to (LV_ENDO union RV_ENDO) (0.0)
    essential_boundaries[BASE   -1] = 0;
    essential_boundaries[EPI    -1] = 1;
    essential_boundaries[LV_ENDO-1] = 1;
    essential_boundaries[RV_ENDO-1] = 1;

    nonzero_essential_boundaries[EPI    -1] = 1;
    nonzero_essential_boundaries[LV_ENDO-1] = 0;
    nonzero_essential_boundaries[RV_ENDO-1] = 0;

    zero_essential_boundaries[EPI    -1] = 0;
    zero_essential_boundaries[LV_ENDO-1] = 1;
    zero_essential_boundaries[RV_ENDO-1] = 1;

    laplace(x, pmesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, dim,opts, rank);
}

void laplace_psi_lv(ParGridFunction *x, ParMesh *pmesh, int dim, Options *opts, int rank)
{
    // Define the following three arrays to determine (1) which border surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said border surfaces
    // should be set to zero.
    int nattr = pmesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

    // Laplace PSI_LV;
    // Solve the Laplace equation from LV_ENDO (1.0) to (RV_ENDO union EPI) (0.0)
    essential_boundaries[BASE   -1] = 0;
    essential_boundaries[EPI    -1] = 1;
    essential_boundaries[LV_ENDO-1] = 1;
    essential_boundaries[RV_ENDO-1] = 1;

    nonzero_essential_boundaries[EPI    -1] = 0;
    nonzero_essential_boundaries[LV_ENDO-1] = 1;
    nonzero_essential_boundaries[RV_ENDO-1] = 0;

    zero_essential_boundaries[EPI    -1] = 1;
    zero_essential_boundaries[LV_ENDO-1] = 0;
    zero_essential_boundaries[RV_ENDO-1] = 1;

    laplace(x, pmesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, dim, opts, rank);
}

void laplace_psi_rv(ParGridFunction *x, ParMesh *pmesh, int dim, Options *opts, int rank)
{
    // Define the following three arrays to determine (1) which border surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said border surfaces
    // should be set to zero.
    int nattr = pmesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

    // Laplace PSI_RV;
    // Solve the Laplace equation from RV_ENDO (1.0) to (LV_ENDO union EPI) (0.0)
    essential_boundaries[BASE   -1] = 0;
    essential_boundaries[EPI    -1] = 1;
    essential_boundaries[LV_ENDO-1] = 1;
    essential_boundaries[RV_ENDO-1] = 1;

    nonzero_essential_boundaries[EPI    -1] = 0;
    nonzero_essential_boundaries[LV_ENDO-1] = 0;
    nonzero_essential_boundaries[RV_ENDO-1] = 1;

    zero_essential_boundaries[EPI    -1] = 1;
    zero_essential_boundaries[LV_ENDO-1] = 1;
    zero_essential_boundaries[RV_ENDO-1] = 0;

    laplace(x, pmesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, dim, opts, rank);
}

int main(int argc, char *argv[])
{
    // 1. Initialize MPI and HYPRE
    Mpi::Init();
    int rank = Mpi::WorldRank();
    Hypre::Init();

    Options opts;
    // 2. Parse command-line options

    // Set program defaults
    opts.mesh_file = NULL;
    opts.out = "./out";
    opts.device_config = "cpu";
    opts.order = 1;

    OptionsParser args(argc, argv);
    args.AddOption(&opts.mesh_file, "-m", "--mesh",
            "Mesh file to use (required)");
    args.AddOption(&opts.out, "-o", "--out",
            "Basename of the ouput files. Outputs will be of the form <basename>_*.gf");
    args.AddOption(&opts.device_config, "-d", "--device",
            "Device configuration string, see Device::Configure().");
    args.AddOption(&opts.order, "-o", "--order",
            "Finite element order (polynomial degree) or -1 for "
            "isoparametric space.");
    args.Parse();

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
    ParGridFunction x_psi_epi;
    laplace_psi_epi(&x_psi_epi, &pmesh, dim, &opts, rank);

    {
        string x_psi_epi_out(opts.out);
        x_psi_epi_out += "_psi_epi.gf";
        ofstream x_psi_epi_ofs(x_psi_epi_out.c_str());
        x_psi_epi_ofs.precision(8);
        x_psi_epi.Save(x_psi_epi_ofs);
    }

    // Laplace PSI_LV;
    // Solve the Laplace equation from LV_ENDO (1.0) to (RV_ENDO union EPI) (0.0)
    ParGridFunction x_psi_lv;
    laplace_psi_lv(&x_psi_lv, &pmesh, dim, &opts, rank);

    {
        string x_psi_lv_out(opts.out);
        x_psi_lv_out += "_psi_lv.gf";
        ofstream x_psi_lv_ofs(x_psi_lv_out.c_str());
        x_psi_lv_ofs.precision(8);
        x_psi_lv.Save(x_psi_lv_ofs);
    }

    // Laplace PSI_RV
    // Solve the Laplace equation from RV_ENDO (1.0) to (LV_ENDO union EPI) (0.0)
    ParGridFunction x_psi_rv;
    laplace_psi_rv(&x_psi_rv, &pmesh, dim, &opts, rank);

    {
        string x_psi_rv_out(opts.out);
        x_psi_rv_out += "_psi_rv.gf";
        ofstream x_psi_rv_ofs(x_psi_rv_out.c_str());
        x_psi_rv_ofs.precision(8);
        x_psi_rv.Save(x_psi_rv_ofs);
    }

    // TODO: Laplace GAMMA_AB (apex -> base)
    // TODO: Solve the Laplace equation from BASE (1.0) to APEX (0.0)

    // Save the mesh
    {

        string mesh_out(opts.out);
        mesh_out += ".mesh";
        ofstream mesh_ofs(mesh_out.c_str());
        mesh_ofs.precision(8);
        pmesh.Print(mesh_ofs);

    }
}

