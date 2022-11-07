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
#include "mfem/general/forall.hpp"

#include "calculus.hpp"
#include "util.hpp"

#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

using namespace mfem;

// Cross product of two 3D vectors a and b, store in c.
MFEM_HOST_DEVICE
static void cross(vec3& c, vec3& a, vec3& b)
{
    // c = a x b
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

// Dot product of two quaternions
MFEM_HOST_DEVICE
static double quatdot(quat& q1, quat& q2)
{
    return q1[0] * q2[0]
         + q1[1] * q2[1]
         + q1[2] * q2[2]
         + q1[3] * q2[3];
}

// Set q1 = q2
MFEM_HOST_DEVICE
static void quatcopy(quat& q1, quat& q2)
{
    q1[0] = q2[0];
    q1[1] = q2[1];
    q1[2] = q2[2];
    q1[3] = q2[3];
}

// Negate the values of quaternion q
MFEM_HOST_DEVICE
static void quatnegate(quat& q)
{
    q[0] = -(q[0]);
    q[1] = -(q[1]);
    q[2] = -(q[2]);
    q[3] = -(q[3]);
}

// Multiple quaternion q with scalar a
MFEM_HOST_DEVICE
static void quatmul(quat& q, double a)
{
    q[0] *= a;
    q[1] *= a;
    q[2] *= a;
    q[3] *= a;
}

// Elementwise multiplication q2 into q1
MFEM_HOST_DEVICE
static void quatelemmul(quat& q1, quat& q2)
{
    q1[0] *= q2[0];
    q1[1] *= q2[1];
    q1[2] *= q2[2];
    q1[3] *= q2[3];
}

// Add quaternion q2 to q1
MFEM_HOST_DEVICE
static void quatadd(quat& q1, quat& q2)
{
    q1[0] += q2[0];
    q1[1] += q2[1];
    q1[2] += q2[2];
    q1[3] += q2[3];
}

// Normalize quaternion
MFEM_HOST_DEVICE
static void quatnormalize(quat& q)
{
    double sum = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    if (sum == 0.0) return;
    double m = 1.0 / sqrt(sum);
    q[0] = q[0] * m;
    q[1] = q[1] * m;
    q[2] = q[2] * m;
    q[3] = q[3] * m;
}

// Set vector a = b
MFEM_HOST_DEVICE
void veccopy(vec3& a, vec3& b)
{
    a[0] = b[0];
    a[1] = b[1];
    a[2] = b[2];
}


MFEM_HOST_DEVICE
void vecset(vec3& a, const double *b)
{
    a[0] = b[0];
    a[1] = b[1];
    a[2] = b[2];
}


// Normalize vector a, a = a / ||a||
MFEM_HOST_DEVICE
static void vecnormalize(vec3& a)
{
    const double sum = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (sum == 0.0) return;
    const double m = 1.0 / sqrt(sum);
    a[0] = a[0] * m;
    a[1] = a[1] * m;
    a[2] = a[2] * m;
}

// Dot product of vectors a and b
MFEM_HOST_DEVICE
double vecdot(vec3& a, vec3& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Subtract vector b from a; a = a - b
MFEM_HOST_DEVICE
static void vecsub(vec3& a, vec3& b)
{
    a[0] = a[0] - b[0];
    a[1] = a[1] - b[1];
    a[2] = a[2] - b[2];
}

// Scale vector a; a_i = a_i * b
MFEM_HOST_DEVICE
void vecmul(vec3& a, double b)
{
    a[0] = a[0] * b;
    a[1] = a[1] * b;
    a[2] = a[2] * b;
}

// Negate vector a; a = -a
MFEM_HOST_DEVICE
static void vecnegate(vec3& a)
{
    a[0] = -a[0];
    a[1] = -a[1];
    a[2] = -a[2];
}

static bool vecisnonzero(vec3& a)
{
    const double sum = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (sum > 0.0) {
        return true;
    }
    return false;
}

// Initialize matrix A values to zero
MFEM_HOST_DEVICE
static void matzero(mat3x3& A)
{
    A[0][0] = 0.0; A[0][1] = 0.0; A[0][2] = 0.0;
    A[1][0] = 0.0; A[1][1] = 0.0; A[1][2] = 0.0;
    A[2][0] = 0.0; A[2][1] = 0.0; A[2][2] = 0.0;
}

// Copy matrix B into A
MFEM_HOST_DEVICE
static void matcopy(mat3x3& A, mat3x3& B)
{
    A[0][0] = B[0][0]; A[0][1] = B[0][1]; A[0][2] = B[0][2];
    A[1][0] = B[1][0]; A[1][1] = B[1][1]; A[1][2] = B[1][2];
    A[2][0] = B[2][0]; A[2][1] = B[2][1]; A[2][2] = B[2][2];
}

// Set matrix C to the matrix multiplication of A and B; C = A x B
MFEM_HOST_DEVICE
static void matmul(mat3x3& C, mat3x3& A, mat3x3& B )
{
    C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
    C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
    C[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2];

    C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
    C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
    C[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2];

    C[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0];
    C[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1];
    C[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2];
}

// Convert the 3x3 rotation matrix Q to a quaternion q, using the algorithm
// described in Appendix 1 of
//
//   Shoemake, K. (1985, July). Animating rotation with quaternion curves. In
//   Proceedings of the 12th annual conference on Computer graphics and
//   interactive techniques (pp. 245-254).
MFEM_HOST_DEVICE
void rot2quat(quat& q, mat3x3& Q)
{
    const double M11=Q[0][0], M12=Q[1][0], M13=Q[2][0];
    const double M21=Q[0][1], M22=Q[1][1], M23=Q[2][1];
    const double M31=Q[0][2], M32=Q[1][2], M33=Q[2][2];

    const double w2 = 0.25 * (1 + M11 + M22 + M33);
    const double err = 1e-15;

    double w, x, y, z;

    if (w2 > err) {
        w = sqrt(w2);
        x = (M23 - M32) / (4.0*w);
        y = (M31 - M13) / (4.0*w);
        z = (M12 - M21) / (4.0*w);
    } else {
        w = 0.0;
        const double x2 = -0.5*(M22 + M33);
        if (x2 > err) {
            x = sqrt(x2);
            y = M12 / (2.0*x);
            z = M13 / (2.0*x);
        } else {
            x = 0.0;
            const double y2 = 0.5*(1-M33);
            if (y2 > err) {
                y = sqrt(y2);
                z = M23 / (2.0*y);
            } else {
                y = 0.0;
                z = 1.0;
            }
        }
    }

    q[0] = w;
    q[1] = x;
    q[2] = y;
    q[3] = z;

    // Make sure the quaternion is normalized
    quatnormalize(q);
}

// Convert the quaternion q to a 3x3 rotation matrix Q, using the algorithm
// described in Appendix 1 of
//
//   Shoemake, K. (1985, July). Animating rotation with quaternion curves. In
//   Proceedings of the 12th annual conference on Computer graphics and
//   interactive techniques (pp. 245-254).
MFEM_HOST_DEVICE
void quat2rot(mat3x3& Q, quat& q)
{
    const double w = q[0], x = q[1], y = q[2], z = q[3];

    {
        const double dot = w*w + x*x + y*y + z*z;
        MFEM_ASSERT_KERNEL(dot - 1.0 <= 1e-12, "Quaternion not normalized");
        MFEM_CONTRACT_VAR(dot);
    }

    const double x2 = x*x;
    const double y2 = y*y;
    const double z2 = z*z;

    const double wx = w*x;
    const double wy = w*y;
    const double wz = w*z;

    const double xy = x*y;
    const double xz = x*z;

    const double yz = y*z;

    Q[0][0] = 1.0 - 2.0*y2 - 2.0*z2;
    Q[1][0] =       2.0*xy + 2.0*wz;
    Q[2][0] =       2.0*xz - 2.0*wy;

    Q[0][1] =       2.0*xy - 2.0*wz;
    Q[1][1] = 1.0 - 2.0*x2 - 2.0*z2;
    Q[2][1] =       2.0*yz + 2.0*wx;

    Q[0][2] =       2.0*xz + 2.0*wy;
    Q[1][2] =       2.0*yz - 2.0*wx;
    Q[2][2] = 1.0 - 2.0*x2 - 2.0*y2;

}

// Spherical Linear intERPolation
MFEM_HOST_DEVICE
static void slerp(quat& q, quat& q1, quat& q2, double t)
{
    double dot = quatdot(q1, q2);
    quatcopy(q, q2);

    if (dot < 0) {
        dot = -dot;
        quatnegate(q);
    }

    // Slerp(q1, q2, t) = ((sin(1-t)*theta)/sin(theta))q1 + ((sin(t)*theta)/sin(theta))q2
    // where theta = acos(q1 dot q2)
    if (dot < 1.0 - 1e-4) {
        const double angle = acos(dot);
        const double a = sin(angle * (1-t))/sin(angle);
        const double b = sin(angle * t)/sin(angle);

        quat q1a = {0};
        quatcopy(q1a, q1);
        quatmul(q1a, a);
        quatmul(q, b);
        quatadd(q, q1a);
    } else {
        // Linear interpolation: q = q1(1-t) + q2*t

        quatcopy(q, q1);
        quatmul(q, (1-t));

        quat qt;
        quatcopy(qt, q2);
        quatmul(qt, t);

        quatadd(q, qt);
    }
}

// Given two vectors u and v, u being a vector pointing in
// the apicobasal direction and v a vector pointing in the transmural
// direction, return a 3x3 orthogonal matrix representing the coordinate system
// for assigning fiber orientation within the myocardium.
//
// As defined in Function 2 in the supplementary material of Bayer2012.
MFEM_HOST_DEVICE
void axis(mat3x3& Q, vec3& u, vec3& v)
{

    vec3 e0, e1, e2;

    // e1 = u / ||u||
    veccopy(e1, u);
    vecnormalize(e1);

    // e2 = u / || u ||
    // where u = v - (e0*v)e0
    //
    // Normalize v as an initial guess for e0
    veccopy(e2, v);
    vecnormalize(e2);


    vec3 e1_dot_e2_e1;
    {
        const double e1_dot_e2 = vecdot(e1, e2);
        veccopy(e1_dot_e2_e1, e1);
        vecmul(e1_dot_e2_e1, e1_dot_e2);
    }
    vecsub(e2, e1_dot_e2_e1);
    vecnormalize(e2);

    // e0 = e1 x e2
    cross(e0, e1, e2);

    vecnormalize(e0);

    Q[0][0] = e0[0];
    Q[1][0] = e0[1];
    Q[2][0] = e0[2];

    Q[0][1] = e1[0];
    Q[1][1] = e1[1];
    Q[2][1] = e1[2];

    Q[0][2] = e2[0];
    Q[1][2] = e2[1];
    Q[2][2] = e2[2];
}

// Take the coordinate system Q, in the form of a 3x3 matrix, and the fiber
// orientation angles a(lpha) and b(eta) at a given point in the mesh, and
// return an orthonormal coordinate system (F S T), in the form of a 3x3
// matrix, where F is the longitudinal direction, S is the sheet normal, and T
// is the transverse direction.
//
// As defined in Function 3 in the supplementary material of Bayer2012.
MFEM_HOST_DEVICE
void orient(mat3x3& Q_out, mat3x3& Q, double a, double b)
{
    const double sina = sin(a*PI/180.0);
    const double sinb = sin(b*PI/180.0);
    const double cosa = cos(a*PI/180.0);
    const double cosb = cos(b*PI/180.0);

    // The algorithm outlines the product of three matrices:
    //
    //         | cosa -sina  0 |   | 1    0    0   |
    //     Q x | sina  cosa  0 | x | 0   cosb sinb |
    //         |   0    0    1 |   | 0  -sinb cosb |
    //
    // Since the two rightmost matrices are constant (depending only on the
    // angles a and b), we instead precompute it, and end up with a single
    // matmul:
    //
    //         | cosa  -sina*cosb  -sina*sinb |
    //     Q x | sina   cosa*cosb   cosa*sinb |
    //         | 0        -sinb       cosb   |


    mat3x3 A;

    A[0][0] = cosa; A[0][1] = -sina*cosb; A[0][2] = -sina*sinb;
    A[1][0] = sina; A[1][1] =  cosa*cosb; A[1][2] = cosa*sinb;
    A[2][0] = 0;    A[2][1] = -sinb;      A[2][2] = cosb;

    matmul(Q_out, Q, A);
}

// BIdirectional SLERP
// Linearly interpolate two orhtogonal matrices Qa and Qb to produce a new
// orthogonal matrix Qab, which is determined by the interpolation factor t.
// When t = 0, Qab = Qa, and when t = 1, Qab = Qb.
//
// As defined in Function 4 in the supplementary material of Bayer2012.
MFEM_HOST_DEVICE
void bislerp(mat3x3& Qab, mat3x3& Qa, mat3x3& Qb, double t, double tol)
{
    if (t <= tol) {
        matcopy(Qab, Qa);
        return;
    }

    if (t >= 1.0 - tol) {
        matcopy(Qab, Qb);
        return;
    }

    // Translate the rotation matrices Qa and Qb into quaternions
    quat qa = {0}, qb = {0};
    rot2quat(qa, Qa);
    rot2quat(qb, Qb);

    // Find qm in { ±qa, ±i*qa, ±j*qa, ±k*qa} that maximizes ||qm*qb||

    quat i = {0}; i[1] = 1.0;
    quat j = {0}; j[2] = 1.0;
    quat k = {0}; k[3] = 1.0;

    quat i_qa = {0}; quatcopy(i_qa, i); quatelemmul(i_qa, qa);
    quat j_qa = {0}; quatcopy(j_qa, j); quatelemmul(j_qa, qa);
    quat k_qa = {0}; quatcopy(k_qa, k); quatelemmul(k_qa, qa);

    quat qa_minus;   quatcopy(qa_minus,   qa);   quatmul(qa_minus, -1.0);
    quat i_qa_minus; quatcopy(i_qa_minus, i_qa); quatmul(i_qa_minus, -1.0);
    quat j_qa_minus; quatcopy(j_qa_minus, j_qa); quatmul(j_qa_minus, -1.0);
    quat k_qa_minus; quatcopy(k_qa_minus, k_qa); quatmul(k_qa_minus, -1.0);

    quat *quat_array[8];
    quat_array[0] = &qa;
    quat_array[1] = &qa_minus;
    quat_array[2] = &i_qa;
    quat_array[3] = &i_qa_minus;
    quat_array[4] = &j_qa;
    quat_array[5] = &j_qa_minus;
    quat_array[6] = &k_qa;
    quat_array[7] = &k_qa_minus;

    double max_abs_dot  = 0.0;
    quat qm = {0};
    for (int i = 0; i < 8; i++) {
        quat *v = quat_array[i];
        const double abs_dot = abs(quatdot(qb, (*v)));
        if (abs_dot > max_abs_dot) {
            max_abs_dot = abs_dot;
            quatcopy(qm, (*v)) ;
        }
    }

    // If the angle is very small, i.e. max_dot is very close to one, return Qb.
    if (max_abs_dot > 1-1e-12) {
        matcopy(Qab, Qb);
        return;
    }

    quat q = {0};
    slerp(q, qm, qb, t);
    quatnormalize(q);
    quat2rot(Qab, q);
}

// Calculate the fiber orientation at vertex i
MFEM_HOST_DEVICE
static void calculate_fiber(
    int i,
    const double phi_epi,
    const double phi_lv,
    const double phi_rv,
    vec3& grad_phi_epi,
    vec3& grad_phi_lv,
    vec3& grad_phi_rv,
    vec3& grad_psi_ab,
    double alpha_endo,
    double alpha_epi,
    double beta_endo,
    double beta_epi,
    double *F,
    double *S,
    double *T,
    double tol)
{

    double depth;
    if (phi_lv + phi_rv < 1e-15) {
        depth = 0.5;
    } else {
        depth = phi_rv / (phi_lv + phi_rv);
    }

    const double alpha_s_d = alpha_endo*(1.0-depth) - alpha_endo*depth;
    const double beta_s_d  = beta_endo *(1.0-depth) - beta_endo *depth;

    const double alpha_w_epi = alpha_endo*(1.0-phi_epi) + alpha_epi * phi_epi;
    const double beta_w_epi  = beta_endo *(1.0-phi_epi) + beta_epi  * phi_epi;

    const bool grad_phi_epi_nonzero = vecisnonzero(grad_phi_epi);
    const bool grad_phi_lv_nonzero  = vecisnonzero(grad_phi_lv);
    const bool grad_phi_rv_nonzero  = vecisnonzero(grad_phi_rv);


    mat3x3 Q_lv = {0};
    if (grad_phi_lv_nonzero) {
        vec3 grad_phi_lv_neg;
        {
            veccopy(grad_phi_lv_neg, grad_phi_lv);
            vecnegate(grad_phi_lv_neg);
        }

        mat3x3 T;
        axis(T, grad_psi_ab, grad_phi_lv_neg);
        orient(Q_lv, T, alpha_s_d, beta_s_d);
    }

    mat3x3 Q_rv = {0};
    if (grad_phi_rv_nonzero) {
        mat3x3 T;
        axis(T, grad_psi_ab, grad_phi_rv);
        orient(Q_rv, T, alpha_s_d, beta_s_d);
    }

    mat3x3 Q_epi = {0};
    if (grad_phi_epi_nonzero) {
        mat3x3 T;
        axis(T, grad_psi_ab, grad_phi_epi);
        orient(Q_epi, T, alpha_w_epi, beta_w_epi);
    }

    mat3x3 Q_endo = {0};
    mat3x3 FST = {0};

    const bool in_rv_outer_wall = !grad_phi_lv_nonzero
                               &&  grad_phi_rv_nonzero
                               &&  grad_phi_epi_nonzero;

    const bool in_lv_outer_wall =  grad_phi_lv_nonzero
                               && !grad_phi_rv_nonzero
                               &&  grad_phi_epi_nonzero;

    const bool in_septum        =  grad_phi_lv_nonzero
                               &&  grad_phi_rv_nonzero
                               && !grad_phi_epi_nonzero;

    if (in_rv_outer_wall) {
        // Since the gradient of phi_lv is zero, we know that we are somewhere
        // in the outer wall of the right ventricle. phi_lv should also be zero here.
        MFEM_ASSERT_KERNEL(phi_lv < 1e-15, "phi_lv is not zero.");

        // As phi_lv is zero and we can skip the first bislerp because we
        // already know that it will result in Q_rv.

        bislerp(FST, Q_rv, Q_epi, phi_epi, tol);

    } else if (in_lv_outer_wall) {
        // Since the gradient of phi_rv is zero, we know that we are somewhere
        // in the outer wall of the left ventricle. phi_rv should also be zero here.
        MFEM_ASSERT_KERNEL(phi_rv < 1e-15, "phi_rv is not zero.");

        bislerp(FST, Q_lv, Q_epi, phi_epi, tol);
    } else if (in_septum) {
        // phi_epi should be approximately zero in the septum
        MFEM_ASSERT_KERNEL(phi_epi < 1e-15, "phi_epi is not zero.");

        // As phi_epi is zero, we can skip the second bislerp, because we know
        // it will results in Q_endo.
        bislerp(FST, Q_lv, Q_rv, depth, tol);

    } else if (grad_phi_lv_nonzero && grad_phi_rv_nonzero && grad_phi_epi_nonzero) {

        // All the gradients are nonzero; use the algorithm in the paper
        bislerp(Q_endo, Q_lv, Q_rv, depth, tol);
        bislerp(FST, Q_endo, Q_epi, phi_epi, tol);
    } else {
        // TODO: Eigenvector
#ifdef DEBUG
        std::cout << "\ti=" << i << "\tEIGEN" << std::endl;
#endif
        F[3*i+0] = 1.0;
        F[3*i+1] = 0.0;
        F[3*i+2] = 0.0;

        S[3*i+0] = 0.0;
        S[3*i+1] = 1.0;
        S[3*i+2] = 0.0;

        T[3*i+0] = 0.0;
        T[3*i+1] = 0.0;
        T[3*i+2] = 1.0;
    }


    F[3*i+0] = FST[0][0];
    F[3*i+1] = FST[1][0];
    F[3*i+2] = FST[2][0];

    S[3*i+0] = FST[0][1];
    S[3*i+1] = FST[1][1];
    S[3*i+2] = FST[2][1];

    T[3*i+0] = FST[0][2];
    T[3*i+1] = FST[1][2];
    T[3*i+2] = FST[2][2];

}

// Calculate the fiber orientations in each vertex.
//
// As defined in Algorithm DefineFibers in the supplementary material of Bayer2012.
void define_fibers(
    int n,
    const double *phi_epi,
    const double *phi_lv,
    const double *phi_rv,
    const double *grad_phi_epi,
    const double *grad_phi_lv,
    const double *grad_phi_rv,
    const double *grad_psi_ab,
    double alpha_endo,
    double alpha_epi,
    double beta_endo,
    double beta_epi,
    double *F,
    double *S,
    double *T,
    double tol)
{
    MFEM_FORALL(i, n,
    {
        const double phi_epi_i = CLAMP(phi_epi[i], 0.0, 1.0);
        const double phi_lv_i  = CLAMP(phi_lv[i],  0.0, 1.0);
        const double phi_rv_i  = CLAMP(phi_rv[i],  0.0, 1.0);

        MFEM_ASSERT_KERNEL(abs(phi_epi_i + phi_lv_i + phi_rv_i - 1.0) < 1e-3,
                    "The laplacians epi, lv and rv do not add up to 1.");

        vec3 grad_phi_epi_i;
        vecset(grad_phi_epi_i, &grad_phi_epi[3*i]);

        vec3 grad_phi_lv_i;
        vecset(grad_phi_lv_i, &grad_phi_lv[3*i]);

        vec3 grad_phi_rv_i;
        vecset(grad_phi_rv_i, &grad_phi_rv[3*i]);

        vec3 grad_psi_ab_i;
        vecset(grad_psi_ab_i, &grad_psi_ab[3*i]);

        MFEM_ASSERT_KERNEL(abs(grad_phi_epi_i[0] + grad_phi_lv_i[0] + grad_phi_rv_i[0]) < 1e-3
                        && abs(grad_phi_epi_i[1] + grad_phi_lv_i[1] + grad_phi_rv_i[1]) < 1e-3
                        && abs(grad_phi_epi_i[2] + grad_phi_lv_i[2] + grad_phi_rv_i[2]) < 1e-3,
                        "The gradients do not add up to zero");

        calculate_fiber(i, phi_epi_i, phi_lv_i, phi_rv_i,
                        grad_phi_epi_i, grad_phi_lv_i, grad_phi_rv_i, grad_psi_ab_i,
                        alpha_endo, alpha_epi, beta_endo, beta_epi,
                        F, S, T, tol);
    });
}

