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
        GridFunction *x,
        Mesh *mesh,
        Array<int> &essential_boundaries,
        Array<int> &nonzero_essential_boundaries,
        Array<int> &zero_essential_boundaries,
        int dim,
        Options *opts)
{
    // Define a finite element space on the mesh.
    FiniteElementCollection *fec;

    fec = new H1_FECollection(opts->order, dim);
    FiniteElementSpace * fespace = new FiniteElementSpace(mesh, fec);
    cout << "Number of finite element unknowns: "
         << fespace->GetTrueVSize() << endl;

    // Determine the list of true (i.e. parallel conforming) essential boundary
    // dofs, defined by the boundary attributes marked as essential (Dirichlet)
    // and converting to a list of true dofs..
    Array<int> ess_tdof_list;
    MFEM_ASSERT(mesh->bdr_attributes.Size() != 0, "Boundary size cannot be zero.");

    fespace->GetEssentialTrueDofs(essential_boundaries, ess_tdof_list);

    // Set up the linear form b(.) which corresponds to the right-hand
    // side of the FEM linear system, which in this case is (0,phi_i) where
    // phi_i are the basis functions in fespace.
    LinearForm b(fespace);
    ConstantCoefficient zero(0.0);
    b.AddDomainIntegrator(new DomainLFIntegrator(zero));
    b.Assemble();

    // Define the solution vector x as a finite element grid function
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
    BilinearForm a(fespace);
    ConstantCoefficient one(1.0);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));

    // Assemble the parallel bilinear form and the corresponding linear system.
    a.Assemble();

    OperatorPtr A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, *x, b, A, X, B);

    // Solve the linear system A X = B.
    // Use a simple symmetric Gauss-Seidel preconditioner with PCG.
    GSSmoother M((SparseMatrix&)(*A));
    PCG(*A, M, B, X, 1, 200, 1e-12, 0.0);

    // Recover the grid function corresponding to X.
    a.RecoverFEMSolution(X, b, *x);
}

void laplace_phi_epi(GridFunction *x, Mesh *mesh, int dim, Options *opts)
{
    // Define the following three arrays to determine (1) which border surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said border surfaces
    // should be set to zero.
    int nattr = mesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

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

    laplace(x, mesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, dim, opts);
}

void laplace_phi_lv(GridFunction *x, Mesh *mesh, int dim, Options *opts)
{
    // Define the following three arrays to determine (1) which border surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said border surfaces
    // should be set to zero.
    int nattr = mesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

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

    laplace(x, mesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, dim, opts);
}

void laplace_phi_rv(GridFunction *x, Mesh *mesh, int dim, Options *opts)
{
    // Define the following three arrays to determine (1) which border surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said border surfaces
    // should be set to zero.
    int nattr = mesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

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

    laplace(x, mesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, dim, opts);
}

int main(int argc, char *argv[])
{

    Options opts;
    // 1. Parse command-line options

    // Set program defaults
    opts.mesh_file = NULL;
    opts.out = "./out";
    opts.device_config = "hip";
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
        args.PrintUsage(cout);
        exit(1);
    }

    args.PrintOptions(cout);

    // 2. Enable hardware devices such as GPUs, and programming models such as
    //    HIP, CUDA, OCCA, RAJA and OpenMP based on command line options.
    Device device(opts.device_config);
    device.Print();

    // 4. Load the mesh
    Mesh mesh(opts.mesh_file, 1, 1);
    int dim = mesh.Dimension();

#ifdef DEBUG
    {
        // Check that the boundary elements have to proper attributes.
        //                   Base Epi  Lv   Rv
        int attributes[4] = {0,   0,   0,   0};

        int num_boundary_elements = mesh.GetNBE();
        for (int i = 0; i < num_boundary_elements; i++) {
            Element *ele = mesh.GetBdrElement(i);
            int attr = ele->GetAttribute();
            MFEM_ASSERT(attr >= BASE && attr <=RV_ENDO,
                    "The element attributes should be 1 (base), 2 (epicardium), "
                    "3 (left ventricle endocardium) or 4 (right ventricle endocardium)");
            attributes[attr-1]++;
        }

        cout << "Number of boundary elements: " << num_boundary_elements << endl
             << "Number of boundary elements at each surface:" << endl
             << "\tBase: " << attributes[0] << " elements" << endl
             << "\tEpi:  " << attributes[1] << " elements" << endl
             << "\tLv: "   << attributes[2] << " elements" << endl
             << "\tRv: "   << attributes[3] << " elements" << endl;
    }
#endif

    // Laplace phi_EPI:
    // Solve the Laplace equation from EPI (1.0) to (LV_ENDO union RV_ENDO) (0.0)
    GridFunction x_phi_epi;
    laplace_phi_epi(&x_phi_epi, &mesh, dim, &opts);

    {
        string x_phi_epi_out(opts.out);
        x_phi_epi_out += "_phi_epi.gf";
        ofstream x_phi_epi_ofs(x_phi_epi_out.c_str());
        x_phi_epi_ofs.precision(8);
        x_phi_epi.Save(x_phi_epi_ofs);
    }

    // Laplace phi_LV;
    // Solve the Laplace equation from LV_ENDO (1.0) to (RV_ENDO union EPI) (0.0)
    GridFunction x_phi_lv;
    laplace_phi_lv(&x_phi_lv, &mesh, dim, &opts);

    {
        string x_phi_lv_out(opts.out);
        x_phi_lv_out += "_phi_lv.gf";
        ofstream x_phi_lv_ofs(x_phi_lv_out.c_str());
        x_phi_lv_ofs.precision(8);
        x_phi_lv.Save(x_phi_lv_ofs);
    }

    // Laplace phi_RV
    // Solve the Laplace equation from RV_ENDO (1.0) to (LV_ENDO union EPI) (0.0)
    GridFunction x_phi_rv;
    laplace_phi_rv(&x_phi_rv, &mesh, dim, &opts);

    {
        string x_phi_rv_out(opts.out);
        x_phi_rv_out += "_phi_rv.gf";
        ofstream x_phi_rv_ofs(x_phi_rv_out.c_str());
        x_phi_rv_ofs.precision(8);
        x_phi_rv.Save(x_phi_rv_ofs);
    }

    // TODO: Laplace GAMMA_AB (apex -> base)
    // TODO: Solve the Laplace equation from BASE (1.0) to APEX (0.0)

    // Save the mesh
    {
        string mesh_out(opts.out);
        mesh_out += ".mesh";
        ofstream mesh_ofs(mesh_out.c_str());
        mesh_ofs.precision(8);
        mesh.Print(mesh_ofs);

    }
}

