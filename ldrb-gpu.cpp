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
#include <limits>

#include "mfem.hpp"
#include "calculus.hpp"
#include "util.hpp"
#include "ldrb-gpu.hpp"

using namespace mfem;

void laplace(
    GridFunction *x,
    Mesh *mesh,
    Array<int> &essential_boundaries,
    Array<int> &nonzero_essential_boundaries,
    Array<int> &zero_essential_boundaries,
    int apex,
    Options *opts)
{
    // Define a finite element space on the mesh.
    FiniteElementCollection *fec;

    fec = new H1_FECollection(1, mesh->Dimension());
    FiniteElementSpace * fespace = new FiniteElementSpace(mesh, fec);
    if (opts->verbose > 1)
        std::cout << "Number of finite element unknowns: " << fespace->GetTrueVSize() << std::endl;

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

    // For the laplacians involving the apex we need to treat the boundary
    // conditions a little different, as it is to be enforced on a single node
    // rather than a whole boundary surface. To do so, we make sure that the
    // apex is en essential true dof, and then we project the wanted value, in
    // this case 0.0, to only that node.
    if (apex >= 0) {
        // Initialize the internal data needed in the finite element space
        x->FESpace()->BuildDofToArrays();

        // Make sure the apex is in the list of essential true Dofs
        ess_tdof_list.Append(apex);

        Array<int> node_disp(1);
        node_disp[0] = apex;
        Vector node_disp_value(1);
        node_disp_value[0] = 0.0;

        VectorConstantCoefficient node_disp_value_coeff(node_disp_value);

        x->ProjectCoefficient(node_disp_value_coeff, node_disp);
    }

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

void laplace_phi_epi(
    GridFunction *x,
    Mesh *mesh,
    Options *opts)
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

    laplace(x, mesh,
            essential_boundaries,
            nonzero_essential_boundaries,
            zero_essential_boundaries,
            -1, opts);
}

void laplace_phi_lv(
    GridFunction *x,
    Mesh *mesh,
    Options *opts)
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
    essential_boundaries[EPI    -1] = 1;
    essential_boundaries[LV_ENDO-1] = 1;
    essential_boundaries[RV_ENDO-1] = 1;

    nonzero_essential_boundaries = 0;
    nonzero_essential_boundaries[LV_ENDO-1] = 1;

    zero_essential_boundaries = 0;
    zero_essential_boundaries[EPI    -1] = 1;
    zero_essential_boundaries[RV_ENDO-1] = 1;

    laplace(x, mesh,
            essential_boundaries,
            nonzero_essential_boundaries,
            zero_essential_boundaries,
            -1, opts);
}

void laplace_phi_rv(
    GridFunction *x,
    Mesh *mesh,
    Options *opts)
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
    essential_boundaries[EPI    -1] = 1;
    essential_boundaries[LV_ENDO-1] = 1;
    essential_boundaries[RV_ENDO-1] = 1;

    nonzero_essential_boundaries = 0;
    nonzero_essential_boundaries[RV_ENDO-1] = 1;

    zero_essential_boundaries = 0;
    zero_essential_boundaries[EPI    -1] = 1;
    zero_essential_boundaries[LV_ENDO-1] = 1;

    laplace(x, mesh,
            essential_boundaries,
            nonzero_essential_boundaries,
            zero_essential_boundaries,
            -1, opts);
}

void laplace_psi_ab(
    GridFunction *x,
    Mesh *mesh,
    int apex,
    Options *opts)
{
    // Define the following three arrays to determine (1) which boundary surfaces
    // to include in the Laplace equation, (2) which of said boundaries should be
    // set to a nonzero value (1.0) and (3) which of said boundary surfaces
    // should be set to zero.
    int nattr = mesh->bdr_attributes.Max();
    Array<int> essential_boundaries(nattr);
    Array<int> nonzero_essential_boundaries(nattr);
    Array<int> zero_essential_boundaries(nattr);

    // Only the base is set as an essential boundary. The boundary condition
    // enforced in the apex is taken care of in `laplace` by setting `apex` >= 0.
    essential_boundaries = 0;
    essential_boundaries[BASE-1] = 1;

    nonzero_essential_boundaries = 0;
    nonzero_essential_boundaries[BASE-1] = 1;

    zero_essential_boundaries = 0;

    laplace(x, mesh,
            essential_boundaries,
            nonzero_essential_boundaries,
            zero_essential_boundaries,
            apex, opts);

}

// Find the vertex closest to the prescribed apex, Euclidean distance.
int find_apex_vertex(Mesh *mesh, Vector& apex)
{
    int apex_vertex = 0;
    double distance = std::numeric_limits<double>::max();
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

// Set the gradient in each vertex to be the average of the gradient in the
// centers of the surrounding elements
void get_vertex_gradients(GridFunction& x, Mesh& mesh, Table* v2e, std::vector<Vector>& grads)
{
    MFEM_ASSERT(grads.size() >= mesh.GetNV(), "grads is to small");

    for (int i = 0; i < mesh.GetNV(); i++) {

        const int num_elements = v2e->RowSize(i);
        const int *elements = v2e->GetRow(i);

        Vector vertex_gradient(3);

        for (int j = 0; j < num_elements; j++) {
            Vector grad(3);
            grad= 0.0;
            int el_id = elements[j];
#if 1
            // Calculate the gradient of an element to be the gradient in its center.

            ElementTransformation *tr = x.FESpace()->GetElementTransformation(el_id);
            Geometry::Type geom = tr->GetGeometryType();
            const IntegrationPoint& center = Geometries.GetCenter(geom);
            tr->SetIntPoint(&center);
            x.GetGradient(*tr, grad);
            vertex_gradient += grad;
#else
            // Calculate the gradient of an element to be the average of the
            // gradients in each of its integration points.

            ElementTransformation *tr = x.FESpace()->GetElementTransformation(el_id);
            const IntegrationRule& ir = x.FESpace()->GetFE(elements[j])->GetNodes();
            for (int k = 0; k < ir.GetNPoints(); k++) {
                Vector grad_point(3);
                grad_point = 0.0;
                const IntegrationPoint& ip = ir.IntPoint(k);
                tr->SetIntPoint(&ip);
                x.GetGradient((*tr), grad_point);
                grad += grad_point;
            }
            grad /= ir.GetNPoints();
            vertex_gradient += grad;

#endif
        }
        vertex_gradient /= num_elements;
        grads[i] = vertex_gradient;
    }
}


void define_fibers(
    Mesh *mesh,
    const double *phi_epi,
    const double *phi_lv,
    const double *phi_rv,
    const double *psi_ab,
    std::vector<Vector>& grad_phi_epi,
    std::vector<Vector>& grad_phi_lv,
    std::vector<Vector>& grad_phi_rv,
    std::vector<Vector>& grad_psi_ab,
    double alpha_endo,
    double alpha_epi,
    double beta_endo,
    double beta_epi,
    std::vector<Vector>& F,
    std::vector<Vector>& S,
    std::vector<Vector>& T)
{

    int num_vertices = mesh->GetNV();

    MFEM_ASSERT(grad_phi_epi.size() == num_vertices, "phi_epi is the wrong size");
    MFEM_ASSERT(grad_phi_lv.size()  == num_vertices, "phi_lv is the wrong size");
    MFEM_ASSERT(grad_phi_rv.size()  == num_vertices, "phi_rv is the wrong size");
    MFEM_ASSERT(grad_psi_ab.size()  == num_vertices, "psi_ab is the wrong size");
    MFEM_ASSERT(F.size()            == num_vertices, "F is the wrong size");
    MFEM_ASSERT(S.size()            == num_vertices, "S is the wrong size");
    MFEM_ASSERT(T.size()            == num_vertices, "T is the wrong size");

#define alpha_s(d) (alpha_endo*(1-(d)) - alpha_endo*(d))
#define alpha_w(d) (alpha_endo*(1-(d)) - alpha_epi *(d))
#define beta_s(d)  (beta_endo *(1-(d)) - beta_endo *(d))
#define beta_w(d)  (beta_endo *(1-(d)) - beta_epi  *(d))

    for (int i = 0; i < num_vertices; i++) {
        const double phi_epi_i = phi_epi[i];
        const double phi_lv_i = phi_lv[i];
        const double phi_rv_i = phi_rv[i];

        Vector grad_phi_epi_i = grad_phi_epi[i];
        Vector grad_phi_lv_i  = grad_phi_lv[i];
        Vector grad_phi_rv_i  = grad_phi_rv[i];
        Vector grad_psi_ab_i  = grad_psi_ab[i];

        const double t = phi_rv_i / (phi_lv_i + phi_rv_i);
        DenseMatrix Q_lv(3,3);
        {
            DenseMatrix T(3,3);
            Vector grad_phi_lv_i_neg = grad_phi_lv_i;
            grad_phi_lv_i_neg.Neg();
            axis(T, grad_psi_ab_i, grad_phi_lv_i_neg);
            orient(Q_lv, T, alpha_s(t), beta_s(t));
        }

        DenseMatrix Q_rv(3,3);
        {
            DenseMatrix T(3,3);
            axis(T, grad_psi_ab_i, grad_phi_rv_i);
            orient(Q_lv, T, alpha_s(t), beta_s(t));
        }

        DenseMatrix Q_endo(3,3);
        bislerp(Q_endo, Q_lv, Q_rv, t);

        DenseMatrix Q_epi(3,3);
        {
            DenseMatrix T(3,3);
            axis(T, grad_psi_ab_i, grad_phi_epi_i);
            orient(Q_epi, T, alpha_w(phi_epi_i), beta_w(phi_epi_i));
        }

        DenseMatrix FST(3,3);
        bislerp(FST, Q_endo, Q_epi, phi_epi_i);

        FST.GetColumn(0, F[i]);
        FST.GetColumn(1, S[i]);
        FST.GetColumn(2, T[i]);
    }
}

int main(int argc, char *argv[])
{

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
    args.AddOption(&opts.verbose,
            "-v", "--verbose",
            "Be verbose");
    args.AddOption(&opts.mesh_file,
            "-m", "--mesh",
            "Mesh file to use", true);
    args.AddOption(&opts.prescribed_apex,
            "-a", "--apex",
            "Coordinate of apex, space separated list: 'x y z'.", true);
    args.AddOption(&opts.output_dir,
            "-o", "--out",
            "Directory for output files.");
    args.AddOption(&opts.device_config,
            "-d", "--device",
            "Device configuration string, see Device::Configure().");
    args.AddOption(&opts.paraview,
            "-p",  "--paraview",
            "-np", "--no-paraview",
            "Save data files for ParaView (paraview.org) visualization.");
    args.Parse();

    if (!args.Good()) {
        args.PrintUsage(std::cout);
        exit(1);
    }

    if (opts.verbose > 1)
        args.PrintOptions(std::cout);

    // Set the basename of the mesh
    opts.mesh_basename = remove_extension(basename(std::string(opts.mesh_file)));

    // Make sure the output direcory exists
    mksubdir(opts.output_dir);

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

    if (opts.verbose > 1) {
        std::cout << "Loaded meshfile '" << opts.mesh_file << "' "
             << "consisting of " << mesh.GetNV() << " vertices "
             << "and " << mesh.GetNE() << " elements" << std::endl;
    }

    if (opts.verbose)
        log_timing(std::cout, "Mesh load", timespec_duration(t0, t1));

    // Set the apex node based on the prescribed apex coordinate.
    int apex = 0;
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);
        apex = find_apex_vertex(&mesh, opts.prescribed_apex);
        clock_gettime(CLOCK_MONOTONIC, &t1);

        if (opts.verbose > 1) {
            double *closest_vertex = mesh.GetVertex(apex);
            std::cout << std::setprecision(2)
                 << "Found closest vertex to prescribed apex "
                 << "("  << opts.prescribed_apex[0]
                 << ", " << opts.prescribed_apex[1]
                 << ", " << opts.prescribed_apex[2]
                 << ") at "
                 << "("  << closest_vertex[0]
                 << ", " << closest_vertex[1]
                 << ", " << closest_vertex[2]
                 << ")." << std::endl;

        }
        if (opts.verbose)
            log_timing(std::cout, "Find apex", timespec_duration(t0, t1));
    }

    // Set up the vertex to element table
    Table *vertex_to_element_table = mesh.GetVertexToElementTable();

    // Solve the Laplace equation from EPI (1.0) to (LV_ENDO union RV_ENDO) (0.0)
    // and calculate the gradients
    GridFunction x_phi_epi;
    std::vector<Vector> grad_phi_epi(mesh.GetNV());
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);
        laplace_phi_epi(&x_phi_epi, &mesh, &opts);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "phi_epi", timespec_duration(t0, t1));

        clock_gettime(CLOCK_MONOTONIC, &t0);
        get_vertex_gradients(x_phi_epi, mesh, vertex_to_element_table, grad_phi_epi);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "grad_phi_epi", timespec_duration(t0, t1));

    }

    // Solve the Laplace equation from LV_ENDO (1.0) to (RV_ENDO union EPI) (0.0)
    // and calculate the gradients
    GridFunction x_phi_lv;
    std::vector<Vector> grad_phi_lv(mesh.GetNV());
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);
        laplace_phi_lv(&x_phi_lv, &mesh, &opts);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "phi_lv", timespec_duration(t0, t1));

        clock_gettime(CLOCK_MONOTONIC, &t0);
        get_vertex_gradients(x_phi_lv, mesh, vertex_to_element_table, grad_phi_lv);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "grad_phi_lv", timespec_duration(t0, t1));
    }


    // Solve the Laplace equation from RV_ENDO (1.0) to (LV_ENDO union EPI) (0.0)
    // and calculate the gradients
    GridFunction x_phi_rv;
    std::vector<Vector> grad_phi_rv(mesh.GetNV());
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);
        laplace_phi_rv(&x_phi_rv, &mesh, &opts);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "phi_rv", timespec_duration(t0, t1));

        clock_gettime(CLOCK_MONOTONIC, &t0);
        get_vertex_gradients(x_phi_rv, mesh, vertex_to_element_table, grad_phi_rv);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "grad_phi_rv", timespec_duration(t0, t1));
    }


    // Solve the Laplace equation from BASE (1.0) to APEX (0.0)
    // and calculate the gradients
    GridFunction x_psi_ab;
    std::vector<Vector> grad_psi_ab(mesh.GetNV());
    {
        clock_gettime(CLOCK_MONOTONIC, &t0);
        laplace_psi_ab(&x_psi_ab, &mesh, apex, &opts);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "psi_ab", timespec_duration(t0, t1));

        clock_gettime(CLOCK_MONOTONIC, &t0);
        get_vertex_gradients(x_psi_ab, mesh, vertex_to_element_table, grad_psi_ab);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        if (opts.verbose)
            log_timing(std::cout, "grad_psi_ab", timespec_duration(t0, t1));
    }

    delete vertex_to_element_table;


    const double *phi_epi = x_phi_epi.Read();
    const double *phi_lv  = x_phi_lv.Read();
    const double *phi_rv  = x_phi_rv.Read();
    const double *psi_ab  = x_psi_ab.Read();

    std::vector<Vector> F(mesh.GetNV()); // Longitudinal
    std::vector<Vector> S(mesh.GetNV()); // Sheet normal
    std::vector<Vector> T(mesh.GetNV()); // Transverse

    double alpha_endo =  40.0;
    double alpha_epi  = -50.0;
    double beta_endo  = -65.0;
    double beta_epi   =  25.0;

    clock_gettime(CLOCK_MONOTONIC, &t0);
    define_fibers(
            &mesh,
            phi_epi,      phi_lv,      phi_rv,      psi_ab,
            grad_phi_epi, grad_phi_lv, grad_phi_rv, grad_psi_ab,
            alpha_endo, alpha_epi, beta_endo, beta_epi,
            F, S, T
    );
    clock_gettime(CLOCK_MONOTONIC, &t1);
    if (opts.verbose)
        log_timing(std::cout, "define_fibers", timespec_duration(t0, t1));

    // Save the mesh and solutions
    {
        // Output the normal solutions in the mfem subdirectory
        std::string mfem_output_dir(opts.output_dir);
        mfem_output_dir += "/mfem";
        mksubdir(mfem_output_dir);

#ifdef DEBUG
        std::string debug_dir = mfem_output_dir + "/debug";
        mksubdir(debug_dir);

        debug_print_to_file(grad_phi_lv, debug_dir, "/grad_phi_lv.txt");
        debug_print_to_file(grad_phi_rv, debug_dir, "/grad_phi_rv.txt");
        debug_print_to_file(grad_psi_ab, debug_dir, "/grad_psi_ab.txt");
        debug_print_to_file(F,           debug_dir, "/F.txt");
        debug_print_to_file(S,           debug_dir, "/S.txt");
        debug_print_to_file(T,           debug_dir, "/T.txt");

        GridFunction grad_phi_epi_gf, grad_phi_lv_gf, grad_phi_rv_gf, grad_psi_ab_gf;
        vertex_vector_to_grid_function(&mesh, grad_phi_epi, &grad_phi_epi_gf);
        vertex_vector_to_grid_function(&mesh, grad_phi_lv,  &grad_phi_lv_gf);
        vertex_vector_to_grid_function(&mesh, grad_phi_rv,  &grad_phi_rv_gf);
        vertex_vector_to_grid_function(&mesh, grad_psi_ab,  &grad_psi_ab_gf);

        GridFunction F_gf, S_gf, T_gf;
        vertex_vector_to_grid_function(&mesh, F, &F_gf);
        vertex_vector_to_grid_function(&mesh, S, &S_gf);
        vertex_vector_to_grid_function(&mesh, T, &T_gf);

#endif
        // Save the MFEM mesh
        std::string mesh_out(mfem_output_dir);
        mesh_out += "/" + opts.mesh_basename + ".mesh";
        std::ofstream mesh_ofs(mesh_out.c_str());
        mesh_ofs.precision(8);
        mesh.Print(mesh_ofs);

        // Save the solutions
        save_solution(&x_phi_epi, mfem_output_dir, opts.mesh_basename, "_phi_epi.gf");
        save_solution(&x_phi_lv,  mfem_output_dir, opts.mesh_basename, "_phi_lv.gf");
        save_solution(&x_phi_rv,  mfem_output_dir, opts.mesh_basename, "_phi_rv.gf");
        save_solution(&x_psi_ab,  mfem_output_dir, opts.mesh_basename, "_psi_ab.gf");

#ifdef DEBUG
        save_solution(&grad_phi_epi_gf, mfem_output_dir, opts.mesh_basename, "_grad_phi_epi.gf");
        save_solution(&grad_phi_lv_gf,  mfem_output_dir, opts.mesh_basename, "_grad_phi_lv.gf");
        save_solution(&grad_phi_rv_gf,  mfem_output_dir, opts.mesh_basename, "_grad_phi_rv.gf");
        save_solution(&grad_psi_ab_gf,  mfem_output_dir, opts.mesh_basename, "_grad_psi_ab.gf");
#endif

        // Save in paraview as well
        ParaViewDataCollection *pd = NULL;
        if (opts.paraview) {
            std::string paraview_path(opts.output_dir);
            paraview_path += "/paraview";
            mksubdir(paraview_path);

            pd = new ParaViewDataCollection(opts.mesh_basename, &mesh);
            pd->SetPrefixPath(paraview_path);
            pd->RegisterField("phi epi", &x_phi_epi);
#ifdef DEBUG
            pd->RegisterField("grad phi epi", &grad_phi_epi_gf);
            pd->RegisterField("grad phi lv",  &grad_phi_lv_gf);
            pd->RegisterField("grad phi rv",  &grad_phi_rv_gf);
            pd->RegisterField("grad psi ab",  &grad_psi_ab_gf);
            pd->RegisterField("F", &F_gf);
            pd->RegisterField("S", &S_gf);
            pd->RegisterField("T", &T_gf);
#endif
            pd->RegisterField("phi lv",  &x_phi_lv);
            pd->RegisterField("phi rv",  &x_phi_rv);
            pd->RegisterField("psi ab",  &x_psi_ab);
            pd->SetLevelsOfDetail(1);
            pd->SetDataFormat(VTKFormat::BINARY);
            pd->SetHighOrderOutput(false);
            pd->Save();
        }

    }
}

