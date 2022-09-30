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
#include <limits>
#include <iomanip>
#include <time.h>

using namespace std;
using namespace mfem;


typedef struct {
    int verbose;
    const char *mesh_file;
    const char *out;
    const char *device_config;
    Vector apex;
} Options;

typedef enum {
    BASE   = 1,
    EPI    = 2,
    LV_ENDO = 3,
    RV_ENDO = 4,
    APEX = 5
} MeshAttributes;

double timespec_duration(
    struct timespec t0,
    struct timespec t1)
{
    return (t1.tv_sec - t0.tv_sec) +
        (t1.tv_nsec - t0.tv_nsec) * 1e-9;
}

void log_timing(ostream& out, const char *log_string, double seconds) {
    out << "[" << left << setw(12) << log_string << "]: "
        << right << fixed << setw(12) << setprecision(6)<< seconds << " s" << endl;
}

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

    fec = new H1_FECollection(1, dim);
    FiniteElementSpace * fespace = new FiniteElementSpace(mesh, fec);
    if (opts->verbose > 1)
        cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << endl;

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
    PCG(*A, M, B, X, opts->verbose > 1 ? 1 : 0, 200, 1e-12, 0.0);

    // Recover the grid function corresponding to X.
    a.RecoverFEMSolution(X, b, *x);
}

void laplace_phi_epi(GridFunction *x, Mesh *mesh, int dim, Options *opts)
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
    essential_boundaries[APEX   -1] = 1;
    essential_boundaries[LV_ENDO-1] = 1;
    essential_boundaries[RV_ENDO-1] = 1;

    nonzero_essential_boundaries = 0;
    nonzero_essential_boundaries[EPI-1] = 1;
    nonzero_essential_boundaries[APEX-1] = 1;

    zero_essential_boundaries = 0;
    zero_essential_boundaries[LV_ENDO-1] = 1;
    zero_essential_boundaries[RV_ENDO-1] = 1;

    laplace(x, mesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, dim, opts);
}

void laplace_phi_lv(GridFunction *x, Mesh *mesh, int dim, Options *opts)
{
    // Define the following three arrays to determine (1) which boundary surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said boundary surfaces
    // should be set to zero.
    int nattr = mesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

    // Solve the Laplace equation from LV_ENDO (1.0) to (RV_ENDO union EPI) (0.0)
    essential_boundaries = 0;
    essential_boundaries[APEX   -1] = 1;
    essential_boundaries[EPI    -1] = 1;
    essential_boundaries[LV_ENDO-1] = 1;
    essential_boundaries[RV_ENDO-1] = 1;

    nonzero_essential_boundaries = 0;
    nonzero_essential_boundaries[LV_ENDO-1] = 1;

    zero_essential_boundaries = 0;
    zero_essential_boundaries[EPI    -1] = 1;
    zero_essential_boundaries[APEX   -1] = 1;
    zero_essential_boundaries[RV_ENDO-1] = 1;

    laplace(x, mesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, dim, opts);
}

void laplace_phi_rv(GridFunction *x, Mesh *mesh, int dim, Options *opts)
{
    // Define the following three arrays to determine (1) which boundary surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said boundary surfaces
    // should be set to zero.
    int nattr = mesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

    // Solve the Laplace equation from RV_ENDO (1.0) to (LV_ENDO union EPI) (0.0)
    essential_boundaries = 0;
    essential_boundaries[APEX   -1] = 1;
    essential_boundaries[EPI    -1] = 1;
    essential_boundaries[LV_ENDO-1] = 1;
    essential_boundaries[RV_ENDO-1] = 1;

    nonzero_essential_boundaries = 0;
    nonzero_essential_boundaries[RV_ENDO-1] = 1;

    zero_essential_boundaries = 0;
    zero_essential_boundaries[APEX   -1] = 1;
    zero_essential_boundaries[EPI    -1] = 1;
    zero_essential_boundaries[LV_ENDO-1] = 1;

    laplace(x, mesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, dim, opts);
}

void laplace_psi_ab(GridFunction *x, Mesh *mesh, int dim, Options *opts)
{
    // Define the following three arrays to determine (1) which boundary surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said boundary surfaces
    // should be set to zero.
    int nattr = mesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

    // Find the apex by solving the Laplace equation Delta u = 1 with boundary
    // condition 0 on the base. The apex will be at he maximum of the solution.
    essential_boundaries = 0;
    essential_boundaries[BASE-1] = 1;
    essential_boundaries[APEX-1] = 1;

    nonzero_essential_boundaries = 0;
    nonzero_essential_boundaries[APEX-1] = 1;

    zero_essential_boundaries = 0;
    zero_essential_boundaries[BASE-1] = 1;

    laplace(x, mesh, essential_boundaries, nonzero_essential_boundaries, zero_essential_boundaries, dim, opts);

}

int main(int argc, char *argv[])
{

    Options opts;
    // Parse command-line options

    // Set program defaults
    opts.mesh_file = NULL;
    opts.out = "./out";
    opts.device_config = "hip";
    opts.verbose = 0;
    opts.apex = Vector(3);

    OptionsParser args(argc, argv);
    args.AddOption(&opts.verbose,       "-v", "--verbose", "Be verbose");
    args.AddOption(&opts.mesh_file,     "-m", "--mesh",    "Mesh file to use", true);
    args.AddOption(&opts.out,           "-o", "--out",     "Basename of the ouput files. Outputs will be of the form <basename>_*.gf");
    args.AddOption(&opts.device_config, "-d", "--device",  "Device configuration string, see Device::Configure().");
    args.AddOption(&opts.apex,          "-a", "--apex",    "Coordinate of apex, space separated list: 'x y z'.", true);
    args.Parse();

    if (!args.Good() || opts.mesh_file == NULL) {
        args.PrintUsage(cout);
        exit(1);
    }

    if (opts.verbose > 1)
        args.PrintOptions(cout);

    struct timespec t0, t1;

    // Enable hardware devices such as GPUs, and programming models such as
    // HIP, CUDA, OCCA, RAJA and OpenMP based on command line options.
    Device device(opts.device_config);
    if (opts.verbose > 1)
        device.Print();

    // Load the mesh

    clock_gettime(CLOCK_MONOTONIC, &t0);
    Mesh mesh(opts.mesh_file, 1, 1);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    int dim = mesh.Dimension();

    if (opts.verbose > 1) {
        cout << "Loaded meshfile '" << opts.mesh_file << "' "
             << "consisting of " << mesh.GetNV() << " vertices "
             << "and " << mesh.GetNE() << " elements" << endl;

    }

    if (opts.verbose)
        log_timing(cout, "Mesh load", timespec_duration(t0, t1));

    // Set the APEX boundary based on the prescribed apex coordinate.
    {
        // Find the vertex closest to the prescribed apex, Euclidean distance.

        clock_gettime(CLOCK_MONOTONIC, &t0);

        int apex_vertex_index = 0;
        double distance = numeric_limits<double>::max();
        for (int i = 0; i < mesh.GetNV(); i++) {
            double *vertices = mesh.GetVertex(i);
            double this_distance = (vertices[0]-opts.apex[0]) * (vertices[0]-opts.apex[0])
                + (vertices[1]-opts.apex[1]) * (vertices[1]-opts.apex[1])
                + (vertices[2]-opts.apex[2]) * (vertices[2]-opts.apex[2]);
            if (this_distance < distance) {
                apex_vertex_index = i;
                distance = this_distance;
            }
        }
        clock_gettime(CLOCK_MONOTONIC, &t1);

        if (opts.verbose > 1) {
            double *closest_vertex = mesh.GetVertex(apex_vertex_index);
            cout << "Found closest vertex to prescribed apex "
                 << "(" << opts.apex[0] << ", " << opts.apex[1] << ", " << opts.apex[2] << ")"
                 << " at "
                 << "(" << closest_vertex[0] << ", " << closest_vertex[1] << ", " << closest_vertex[2] << ")" << endl;

        }
        if (opts.verbose)
            log_timing(cout, "Find apex", timespec_duration(t0, t1));


        // Set the apex boundary to be the boundary elements that has the apex vertex as one of its vertices
        for (int i = 0; i < mesh.GetNBE(); i++) {
            Element *ele = mesh.GetBdrElement(i);
            const int *vertices = ele->GetVertices();
            for (int j = 0; j < ele->GetNVertices(); j++) {
                if (vertices[j] == apex_vertex_index) {
                    mesh.SetBdrAttribute(i, APEX);
                }
            }
        }
        mesh.SetAttributes(); // Added the APEX attribute, update the mesh
    }

    // Solve the Laplace equation from EPI (1.0) to (LV_ENDO union RV_ENDO) (0.0)
    clock_gettime(CLOCK_MONOTONIC, &t0);
    GridFunction x_phi_epi;
    laplace_phi_epi(&x_phi_epi, &mesh, dim, &opts);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (opts.verbose)
        log_timing(cout, "phi_epi", timespec_duration(t0, t1));


    // Solve the Laplace equation from LV_ENDO (1.0) to (RV_ENDO union EPI) (0.0)
    clock_gettime(CLOCK_MONOTONIC, &t0);
    GridFunction x_phi_lv;
    laplace_phi_lv(&x_phi_lv, &mesh, dim, &opts);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (opts.verbose)
        log_timing(cout, "phi_lv", timespec_duration(t0, t1));


    // Solve the Laplace equation from RV_ENDO (1.0) to (LV_ENDO union EPI) (0.0)
    clock_gettime(CLOCK_MONOTONIC, &t0);
    GridFunction x_phi_rv;
    laplace_phi_rv(&x_phi_rv, &mesh, dim, &opts);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (opts.verbose)
        log_timing(cout, "phi_rv", timespec_duration(t0, t1));


    // Solve the Laplace equation from BASE (1.0) to APEX (0.0)
    clock_gettime(CLOCK_MONOTONIC, &t0);
    GridFunction x_psi_ab;
    laplace_psi_ab(&x_psi_ab, &mesh, dim, &opts);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (opts.verbose)
        log_timing(cout, "psi_ab", timespec_duration(t0, t1));

    // Save the mesh and solutions
    {
        string mesh_out(opts.out);
        mesh_out += ".mesh";
        ofstream mesh_ofs(mesh_out.c_str());
        mesh_ofs.precision(8);
        mesh.Print(mesh_ofs);

    }

    {
        string x_phi_epi_out(opts.out);
        x_phi_epi_out += "_phi_epi.gf";
        ofstream x_phi_epi_ofs(x_phi_epi_out.c_str());
        x_phi_epi_ofs.precision(8);
        x_phi_epi.Save(x_phi_epi_ofs);
    }

    {
        string x_phi_lv_out(opts.out);
        x_phi_lv_out += "_phi_lv.gf";
        ofstream x_phi_lv_ofs(x_phi_lv_out.c_str());
        x_phi_lv_ofs.precision(8);
        x_phi_lv.Save(x_phi_lv_ofs);
    }

    {
        string x_phi_rv_out(opts.out);
        x_phi_rv_out += "_phi_rv.gf";
        ofstream x_phi_rv_ofs(x_phi_rv_out.c_str());
        x_phi_rv_ofs.precision(8);
        x_phi_rv.Save(x_phi_rv_ofs);
    }

    {
        string x_psi_ab_out(opts.out);
        x_psi_ab_out += "_psi_ab.gf";
        ofstream x_psi_ab_ofs(x_psi_ab_out.c_str());
        x_psi_ab_ofs.precision(8);
        x_psi_ab.Save(x_psi_ab_ofs);
    }
}

