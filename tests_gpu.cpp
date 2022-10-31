#include "mfem.hpp"
#include "calculus_gpu.hpp"

static bool vec3equal(vec3& u, vec3& v) {
    return u[0] == v[0]
        && u[1] == v[1]
        && u[2] == v[2];
}

static void test_vec3equal(
        vec3& u,
        std::string const& u_name,
        vec3& v,
        std::string const& v_name)
{
    if (!vec3equal(u, v)) {
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
    test_vec3equal(e1, "e1", u, "u");


    // Test 2: v == ||v|| * e2
    vec3 v_norm_e2;
    {
        const double v_norm = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        veccopy(v_norm_e2, e2);
        vecmul(v_norm_e2, v_norm);
    }
    std::cout << " Test 2: v == ||v|| * e2: ";
    test_vec3equal(v,"v", v_norm_e2, "||v||*e2");

    // Test 3: e0 dot e1 == 0
    std::cout << " Test 3: e0 dot e1 == 0: ";
    const double e0_dot_e1 = vecdot(e0, e1);
    if (e0_dot_e1 != 0.0) {
        std::cout << "[FAILED] (e0 dot e1 = " << e0_dot_e1 << ")" << std::endl;
    } else {
        std::cout << "[PASSED]" << std::endl;
    }

    // Test 4: e0 dot e2 == 0
    std::cout << " Test 4: e0 dot e2 == 0: ";
    const double e0_dot_e2 = vecdot(e0, e2);
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

int main(void)
{
    test_axis();
    test_quaternions();
    return 0;
}
