#include "mfem.hpp"
#include "calculus.hpp"

static bool vector3D_equal(mfem::Vector u, mfem::Vector v) {
    MFEM_ASSERT(u.Size() == v.Size() && u.Size() == 3, "u and v are not the same size");
    for (int i = 0; i < u.Size(); i++) {
        if (u(i) != v(i)) {
            return false;
        }
    }
    return true;
}

static bool test_vector3D_equal(
        mfem::Vector u,
        std::string const& u_name,
        mfem::Vector v,
        std::string const& v_name)
{
    if (!vector3D_equal(u, v)) {
        std::cout << "[FAILED] ("
                  << u_name << " = {" << u(0) << ", " << u(1) << ", " << u(2) << "}"
                  << ", "
                  << v_name << " = {" << v(0) << ", " << v(1) << ", " << v(2) << "}"
                  << ")" << std::endl;
        return false;
    } else {
        std::cout << "[PASSED]" << std::endl;
        return true;
    }
}

bool test_axis(void)
{

    bool passed = true;
    std::cout << "test_axis:" << std::endl;

    mfem::Vector u(3);
    u(0) = 1.0; u(1) = 0.0; u(2) = 0.0;

    mfem::Vector v(3);
    v(0) = 0.0; v(1) = 0.5; v(2) = 0.0;

    mfem::DenseMatrix Q(3,3);
    axis(Q, u, v);

    mfem::Vector e0(3), e1(3), e2(3);

    Q.GetColumn(0, e0);
    Q.GetColumn(1, e1);
    Q.GetColumn(2, e2);

    // Test 1: e1 == u
    std::cout << " Test 1: e1 == u: ";
    passed = passed && test_vector3D_equal(e1, "e1", u, "u");


    // Test 2: v == ||v|| * e2
    mfem::Vector v_norm_e2(3);
    {
        const double v_norm = v.Norml2();
        v_norm_e2 = e2;
        v_norm_e2 *= v_norm;
    }
    std::cout << " Test 2: v == ||v|| * e2: ";
    passed = passed && test_vector3D_equal(v,"v", v_norm_e2, "||v||*e2");

    // Test 3: e0 dot e1 == 0
    std::cout << " Test 3: e0 dot e1 == 0: ";
    const double e0_dot_e1 = e0*e1;
    if (e0_dot_e1 != 0.0) {
        std::cout << "[FAILED] (e0 dot e1 = " << e0_dot_e1 << ")" << std::endl;
    } else {
        std::cout << "[PASSED]" << std::endl;
    }

    // Test 4: e0 dot e2 == 0
    std::cout << " Test 4: e0 dot e2 == 0: ";
    const double e0_dot_e2 = e0*e2;
    if (e0_dot_e2 != 0.0) {
        std::cout << "[FAILED] (e0 dot e2 = " << e0_dot_e2 << ")" << std::endl;
    } else {
        std::cout << "[PASSED]" << std::endl;
    }

    return passed;
}

static bool Matrix3D_equal(mfem::DenseMatrix& A, mfem::DenseMatrix& B, double tol)
{
    MFEM_ASSERT(A.Height() == A.Width() && A.Height() == 3, "A is wrong size.")
    MFEM_ASSERT(B.Height() == B.Width() && B.Height() == 3, "A is wrong size.")
    return abs(A(0,0) - B(0,0)) < tol
        && abs(A(0,1) - B(0,1)) < tol
        && abs(A(0,2) - B(0,2)) < tol
        && abs(A(1,0) - B(1,0)) < tol
        && abs(A(1,1) - B(1,1)) < tol
        && abs(A(1,2) - B(1,2)) < tol
        && abs(A(2,0) - B(2,0)) < tol
        && abs(A(2,1) - B(2,1)) < tol
        && abs(A(2,2) - B(2,2)) < tol;
}

static void Matrix3D_print(mfem::DenseMatrix& A)
{
    std::cout << std::setw(6) << std::setprecision(3)
              << "\t" << A(0,0) << " " << A(0,1) << " " << A(0,2) << std::endl
              << "\t" << A(1,0) << " " << A(1,1) << " " << A(1,2) << std::endl
              << "\t" << A(2,0) << " " << A(2,1) << " " << A(2,2) << std::endl;
}

bool test_quaternions(void)
{

    std::cout << "test_quaternions" << std::endl;
    double a = 30.0;

    const double sina = sin(a*PI/180);
    const double cosa = cos(a*PI/180);

    const double A_vals[3][3] = {
        {cosa, -sina, 0},
        {sina,  cosa, 0},
        {0,     0,    1}
    };

    mfem::DenseMatrix A(A_vals);
    mfem::DenseMatrix B(3, 3);
    mfem::Vector q(4);

    rot2quat(q, A);
    quat2rot(B, q);

    std::cout << " Test 1: A == quat2rot(rot2quat(A)): ";
    if (!Matrix3D_equal(A, B, 1e-12)) {
        std::cout << "[FAILED]" << std::endl;
        std::cout << " A = " << std::endl;
        Matrix3D_print(A);
        std::cout << " q = {" << q(0) << ", " << q(1) << ", " << q(2) << ", " << q(3) << "}" << std::endl;
        std::cout << " B = " << std::endl;
        Matrix3D_print(B);
    } else {
        std::cout << "[PASSED]" << std::endl;
    }

    return true;
}

int main(void)
{
    test_axis();
    test_quaternions();
    return 0;
}
