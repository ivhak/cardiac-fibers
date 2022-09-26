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


int main(int argc, char *argv[])
{
    // 1. Initialize MPI and HYPRE
    Mpi::Init();
    int nprocs = Mpi::WorldSize();
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
    //    CUDA, OCCA, RAJA and OpenMP based on command line options.
    Device device(opts.device_config);
    if (rank == 0)
        device.Print();

    // 4. Load the mesh
    Mesh mesh(opts.mesh_file, 1, 0);
    int dim = mesh.Dimension();

    //                   Base Epi  Lv   Rv
    int attributes[4] = {0,   0,   0,   0};

    int num_border_elements = mesh.GetNBE();
    for (int i = 0; i < num_border_elements; i++) {
        Element *ele = mesh.GetBdrElement(i);
        int attr = ele->GetAttribute();
        MFEM_ASSERT(attr >= 1 && attr <=4,
                "The element attributes should be 1 (base), 2 (epicardium), "
                "3 (left ventricle endocardium) or 4 (right ventricle endocardium)");
        attributes[attr-1]++;
    }

#ifdef DEBUG
    cout << "Number of border elements: " << num_border_elements << endl
         << "Number of border elements at each surface:" << endl
         << "\tBase: " << attributes[0] << " elements" << endl
         << "\tEpi:  " << attributes[1] << " elements" << endl
         << "\tLv: "   << attributes[2] << " elements" << endl
         << "\tRv: "   << attributes[3] << " elements" << endl;
#endif

    // 6. Define a parallel mesh by a partitioning of the serial mesh.
    ParMesh pmesh(MPI_COMM_WORLD, mesh);
    mesh.Clear();

    // 7. Define a finite element space on the mesh.
    FiniteElementCollection *fec;

    fec = new H1_FECollection(opts.order, dim);
    ParFiniteElementSpace fespace(&pmesh, fec);
    HYPRE_BigInt size = fespace.GlobalTrueVSize();
    if (rank == 0)
        cout << "Number of finite element unknowns: " << size << endl;

    // 8. Determine the list of true (i.e. parallel conforming) essential
    //    boundary dofs. In this example, the boundary conditions are defined
    //    by marking all the boundary attributes from the mesh as essential
    //    (Dirichlet) and converting them to a list of true dofs.
    Array<int> ess_tdof_list;
    if (pmesh.bdr_attributes.Size())
    {
        Array<int> ess_bdr(pmesh.bdr_attributes.Max());
        ess_bdr = 1;
        fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    }

    // 9. Set up the parallel linear form b(.) which corresponds to the
    //    right-hand side of the FEM linear system, which in this case is
    //    (1,phi_i) where phi_i are the basis functions in fespace.
    ParLinearForm b(&fespace);
    ConstantCoefficient one(1.0);
    b.AddDomainIntegrator(new DomainLFIntegrator(one));
    b.Assemble();

    // 10. Define the solution vector x as a parallel finite element grid
    //     function corresponding to fespace. Initialize x with initial guess of
    //     zero, which satisfies the boundary conditions.
    ParGridFunction x(&fespace);
    x = 0.0;

    // 11. Set up the parallel bilinear form a(.,.) on the finite element space
    //     corresponding to the Laplacian operator -Delta, by adding the
    //     Diffusion domain integrator.
    ParBilinearForm a(&fespace);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));

    // 12. Assemble the parallel bilinear form and the corresponding linear
    //     system, applying any necessary transformations such as: parallel
    //     assembly, eliminating boundary conditions, applying conforming
    //     constraints for non-conforming AMR, static condensation, etc.
    a.Assemble();

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

    // 13. Solve the linear system A X = B.
    //     * With full assembly, use the BoomerAMG preconditioner from hypre.
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

    // 14. Recover the parallel grid function corresponding to X. This is the
    //     local finite element solution on each processor.
    a.RecoverFEMSolution(X, b, x);

    // 15. Save the refined mesh and the solution in parallel. This output can
    //     be viewed later using GLVis: "glvis -np <np> -m mesh -g sol".
    {
        ofstream mesh_ofs(opts.out_mesh);
        mesh_ofs.precision(8);
        pmesh.Print(mesh_ofs);

        ofstream sol_ofs(opts.out_sol);
        sol_ofs.precision(8);
        x.Save(sol_ofs);
    }
}
