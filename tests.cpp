// Copyright (C) 2022 Iver Håkonsen
//
// cardiac-fibers is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// cardiac-fibers is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License along with
// cardiac-fibers.  If not, see <https://www.gnu.org/licenses/>.
//
// Authors: Iver Håkonsen <hakonseniver@yahoo.no

#include "mfem.hpp"

#include "fem.hpp"
#include "calculus.hpp"

using namespace mfem;

static bool vec3_equal(vec3& u, vec3& v) {
    return u[0] == v[0]
        && u[1] == v[1]
        && u[2] == v[2];
}

static void test_vec3_equal(
        vec3& u,
        std::string const& u_name,
        vec3& v,
        std::string const& v_name)
{
    if (!vec3_equal(u, v)) {
        std::cout << "[FAILED] ("
                  << u_name << " = {" << u[0] << ", " << u[1] << ", " << u[2] << "}"
                  << ", "
                  << v_name << " = {" << v[0] << ", " << v[1] << ", " << v[2] << "}"
                  << ")" << std::endl;
    } else {
        std::cout << "[PASSED]" << std::endl;
    }
}

bool test_axis(void)
{

    bool passed = true;
    std::cout << "test_axis:" << std::endl;

    vec3 u;
    u[0] = 1.0; u[1] = 0.0; u[2] = 0.0;

    vec3 v;
    v[0] = 0.0; v[1] = 0.5; v[2] = 0.0;

    mat3x3 Q;
    axis(Q, u, v);

    vec3 e0, e1, e2;

    e0[0] = Q[0][0];
    e0[1] = Q[1][0];
    e0[2] = Q[2][0];

    e1[0] = Q[0][1];
    e1[1] = Q[1][1];
    e1[2] = Q[2][1];

    e2[0] = Q[0][2];
    e2[1] = Q[1][2];
    e2[2] = Q[2][2];

    // Test 1: e1 == u
    std::cout << " Test 1: e1 == u: ";
    test_vec3_equal(e1, "e1", u, "u");


    // Test 2: v == ||v|| * e2
    vec3 v_norm_e2;
    {
        const double v_norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        v_norm_e2 = e2;
        v_norm_e2 *= v_norm;
    }
    std::cout << " Test 2: v == ||v|| * e2: ";
    test_vec3_equal(v,"v", v_norm_e2, "||v||*e2");

    // Test 3: e0 dot e1 == 0
    std::cout << " Test 3: e0 dot e1 == 0: ";
    const double e0_dot_e1 = vec3_dot(e0, e1);
    if (e0_dot_e1 != 0.0) {
        std::cout << "[FAILED] (e0 dot e1 = " << e0_dot_e1 << ")" << std::endl;
    } else {
        std::cout << "[PASSED]" << std::endl;
    }

    // Test 4: e0 dot e2 == 0
    std::cout << " Test 4: e0 dot e2 == 0: ";
    const double e0_dot_e2 = vec3_dot(e0, e2);
    if (e0_dot_e2 != 0.0) {
        std::cout << "[FAILED] (e0 dot e2 = " << e0_dot_e2 << ")" << std::endl;
    } else {
        std::cout << "[PASSED]" << std::endl;
    }

    return passed;
}

static bool mat3x3equal(mat3x3& A, mat3x3& B, double tol)
{
    return abs(A[0][0] - B[0][0]) < tol
        && abs(A[0][1] - B[0][1]) < tol
        && abs(A[0][2] - B[0][2]) < tol
        && abs(A[1][0] - B[1][0]) < tol
        && abs(A[1][1] - B[1][1]) < tol
        && abs(A[1][2] - B[1][2]) < tol
        && abs(A[2][0] - B[2][0]) < tol
        && abs(A[2][1] - B[2][1]) < tol
        && abs(A[2][2] - B[2][2]) < tol;
}

static void mat3x3print(mat3x3& A)
{
    std::cout << std::setw(6) << std::setprecision(3)
              << "\t" << A[0][0] << " " << A[0][1] << " " << A[0][2] << std::endl
              << "\t" << A[1][0] << " " << A[1][1] << " " << A[1][2] << std::endl
              << "\t" << A[2][0] << " " << A[2][1] << " " << A[2][2] << std::endl;
}

bool test_quaternions(void)
{

    std::cout << "test_quaternions" << std::endl;
    double a = 30.0;

    const double sina = sin(a*PI/180);
    const double cosa = cos(a*PI/180);

    mat3x3 A;
    A[0][0] = cosa; A[0][1] = -sina; A[0][2] = 0;
    A[1][0] = sina; A[1][1] =  cosa; A[1][2] = 0;
    A[2][0] = 0;    A[2][1] =  0;    A[2][2] = 1;

    mat3x3 B;
    quat q;

    rot2quat(q, A);
    quat2rot(B, q);

    std::cout << " Test 1: A == quat2rot(rot2quat(A)): ";
    if (!mat3x3equal(A, B, 1e-12)) {
        std::cout << "[FAILED]" << std::endl;
        std::cout << " A = " << std::endl;
        mat3x3print(A);
        std::cout << " q = {" << q[0] << ", " << q[1] << ", " << q[2] << ", " << q[3] << "}" << std::endl;
        std::cout << " B = " << std::endl;
        mat3x3print(B);
    } else {
        std::cout << "[PASSED]" << std::endl;
    }

    return true;
}

void solve_laplace_epi(GridFunction& x, Mesh& mesh)
{
    int nattr = mesh.bdr_attributes.Max();
    Array<int> ess_bdr(nattr);          // Essential boundaries
    Array<int> nonzero_ess_bdr(nattr);  // Essential boundaries with value 1.0
    Array<int> zero_ess_bdr(nattr);     // Essential boundaries with value 0.0

    const int epi  = 2;
    const int lv   = 3;
    const int rv   = 4;

    ess_bdr = 0;
    ess_bdr[epi-1] = 1;
    ess_bdr[lv -1] = 1;
    ess_bdr[rv -1] = 1;

    nonzero_ess_bdr = 0;
    nonzero_ess_bdr[epi-1] = 1;

    zero_ess_bdr = 0;
    zero_ess_bdr[lv-1] = 1;
    zero_ess_bdr[rv-1] = 1;

    ConstantCoefficient zero(0.0);
    ConstantCoefficient one(1.0);

    x = 0.0;
    x.ProjectBdrCoefficient(zero, zero_ess_bdr);
    x.ProjectBdrCoefficient(one, nonzero_ess_bdr);

    Array<int> ess_tdof_list;
    x.FESpace()->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    LinearForm b(x.FESpace());
    b.AddDomainIntegrator(new DomainLFIntegrator(zero));
    b.Assemble();

    BilinearForm a(x.FESpace());
    a.AddDomainIntegrator(new DiffusionIntegrator(one));
    a.Assemble();

    SparseMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

    GSSmoother M(A);
    PCG(A, M, B, X, 0, 200, 1e-12, 0.0);

    a.RecoverFEMSolution(X, b, x);
}

void test_gradients(void)
{
    std::cout << "test_gradients: ";

    Mesh mesh("./mesh/gmsh/heart01.msh", 1, 1);

    H1_FECollection h1_fec(1, mesh.Dimension());
    L2_FECollection l2_fec(0, mesh.Dimension());

    FiniteElementSpace fespace_scalar_h1(&mesh, &h1_fec);
    FiniteElementSpace fespace_vector_l2(&mesh, &l2_fec, 3, Ordering::byVDIM);

    // Solve the Laplacian from epi to endo on heart01
    GridFunction x(&fespace_scalar_h1);
    solve_laplace_epi(x, mesh);

    GridFunction grad_manual(&fespace_vector_l2);
    {
        Vector vertices;
        mesh.GetVertices(vertices);
        const double *vert = vertices.Read();

        const Table& h1_element_to_dof = fespace_scalar_h1.GetElementToDofTable();
        const Table& l2_element_to_dof = fespace_vector_l2.GetElementToDofTable();

        const int *h1_I = h1_element_to_dof.ReadI();
        const int *h1_J = h1_element_to_dof.ReadJ();

        const int *l2_I = l2_element_to_dof.ReadI();
        const int *l2_J = l2_element_to_dof.ReadJ();

        const int num_elements = mesh.GetNE();
        const int num_vertices = mesh.GetNV();

        const double *x_vals = x.Read();
        double *grad_vals = grad_manual.Write();

        compute_gradient(grad_vals, x_vals, vert,
                         num_elements, num_vertices,
                         h1_I, h1_J, l2_I, l2_J);

    }

    GridFunction grad_mfem(&fespace_vector_l2);
    {
        GradientGridFunctionCoefficient ggfc(&x);
        grad_mfem.ProjectCoefficient(ggfc);
    }

    const double *grad_manual_vals = grad_manual.Read();
    const double *grad_mfem_vals = grad_mfem.Read();

    bool passed = true;
    double max_diff = 0.0;
    for (int i = 0; i < grad_manual.Size(); i++) {
        double diff = abs(grad_manual_vals[i] - grad_mfem_vals[i]);
        if (diff > 1e-12) {
            passed = false;
            break;
        }
        if (max_diff < diff) {
            max_diff = diff;
        }
    }
    if (passed) {
        std::cout << "[PASSED] (Maximal elementwise diff: " << max_diff << ")"<< std::endl;

    } else {
        std::cout << "[FAILED]" << std::endl;
    }
}

void test_projection(void)
{
    std::cout << "test_projection: ";

    Mesh mesh("./mesh/gmsh/heart01.msh", 1, 1);

    H1_FECollection h1_fec(1, mesh.Dimension());
    L2_FECollection l2_fec(0, mesh.Dimension());

    FiniteElementSpace fespace_scalar_h1(&mesh, &h1_fec);
    FiniteElementSpace fespace_scalar_l2(&mesh, &l2_fec);

    // Solve the Laplacian from epi to endo on heart01
    GridFunction x(&fespace_scalar_h1);
    solve_laplace_epi(x, mesh);

    GridFunction x_l2_manual(&fespace_scalar_l2);
    {
        const Table& h1_element_to_dof = fespace_scalar_h1.GetElementToDofTable();
        const Table& l2_element_to_dof = fespace_scalar_l2.GetElementToDofTable();

        const int *h1_I = h1_element_to_dof.ReadI();
        const int *h1_J = h1_element_to_dof.ReadJ();

        const int *l2_I = l2_element_to_dof.ReadI();
        const int *l2_J = l2_element_to_dof.ReadJ();

        const int num_elements = mesh.GetNE();

        const double *x_h1 = x.Read();
        double *x_l2 = x_l2_manual.Write();

        project_h1_to_l2(x_l2, x_h1, num_elements,
                         h1_I, h1_J, l2_I, l2_J);
    }

    GridFunction x_l2_mfem(&fespace_scalar_l2);
    {
        GridFunctionCoefficient gfc(&x);
        x_l2_mfem.ProjectCoefficient(gfc);
    }

    const double *x_l2_manual_vals = x_l2_manual.Read();
    const double *x_l2_mfem_vals = x_l2_mfem.Read();

    bool passed = true;
    double max_diff = 0.0;
    for (int i = 0; i < x_l2_manual.Size(); i++) {
        double diff = abs(x_l2_manual_vals[i] - x_l2_mfem_vals[i]);
        if (diff > 1e-12) {
            passed = false;
            break;
        }
        if (max_diff < diff) {
            max_diff = diff;
        }
    }
    if (passed) {
        std::cout << "[PASSED] (Maximal elementwise diff: " << max_diff << ")"<< std::endl;

    } else {
        std::cout << "[FAILED]" << std::endl;
    }
}


int main(void)
{
    test_axis();
    test_quaternions();
    test_gradients();
    test_projection();
    return 0;
}
