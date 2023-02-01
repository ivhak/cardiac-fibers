// Copyright (C) 2022 Iver Håkonsen
//
// cardiac-fibers is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// cardiac-fibers is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// cardiac-fibers.  If not, see <https://www.gnu.org/licenses/>.
//
// Authors: Iver Håkonsen <hakonseniver@yahoo.no

#include "mfem.hpp"
#include "mfem/general/forall.hpp"

#include "calculus.hpp"
#include "util.hpp"

#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

using namespace mfem;

// Dot product of two quaternions
MFEM_HOST_DEVICE
static double quat_dot(quat& q1, quat& q2)
{
    return q1[0] * q2[0]
         + q1[1] * q2[1]
         + q1[2] * q2[2]
         + q1[3] * q2[3];
}


// Normalize quaternion
MFEM_HOST_DEVICE
static void quat_normalize(quat& q)
{
    double sum = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    if (sum == 0.0) return;
    double m = 1.0 / sqrt(sum);
    q[0] = q[0] * m;
    q[1] = q[1] * m;
    q[2] = q[2] * m;
    q[3] = q[3] * m;
}

// vec3_cross product of two 3D vectors a and b, store in c.
MFEM_HOST_DEVICE
void vec3_cross(vec3& c, vec3& a, vec3& b)
{
    // c = a x b
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

// Dot product of vectors a and b
MFEM_HOST_DEVICE
double vec3_dot(vec3& a, vec3& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

MFEM_HOST_DEVICE
static void vec3_set_from_ptr(vec3& a, const double *b)
{
    a[0] = b[0];
    a[1] = b[1];
    a[2] = b[2];
}


// Normalize vector a, a = a / ||a||
MFEM_HOST_DEVICE
static void vec3_normalize(vec3& a)
{
    const double sum = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    if (sum == 0.0) return;
    const double m = 1.0 / sqrt(sum);
    a[0] = a[0] * m;
    a[1] = a[1] * m;
    a[2] = a[2] * m;
}


// Negate vector a; a = -a
MFEM_HOST_DEVICE
static void vec3_negate(vec3& a)
{
    a[0] = -a[0];
    a[1] = -a[1];
    a[2] = -a[2];
}

// Calculate the magnitude of a vec3
MFEM_HOST_DEVICE
static double vec3_magnitude(vec3& v){
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

// Set matrix C to the matrix multiplication of A and B; C = A x B
MFEM_HOST_DEVICE
static void mat3x3_mul(mat3x3& C, mat3x3& A, mat3x3& B )
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
    quat_normalize(q);
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
    double dot = quat_dot(q1, q2);
    // If the dot product is close to one, the angle approaches zero, and slerp
    // reduces to a regular linear interpolation.
    if (dot > 1-1e-12) {
        // Slerp(q1, q2, t) = ((sin(1-t)*theta)/sin(theta))q1 + ((sin(t)*theta)/sin(theta))q2
        // where theta = acos(q1 dot q2)
        const double angle = acos(dot);
        const double a = sin(angle * (1-t))/sin(angle);
        const double b = sin(angle * t)/sin(angle);

        q = q2;
        q *= b;
        quat q1a = q1;
        q1a *= a;
        q += q1a;
    } else {
        q = q1;
        q *= (1-t);
        quat q2t = q2;
        q2t *= t;
        q += q2t;
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
    e1 = u;
    vec3_normalize(e1);

    // e2 = v - (e1*v)e1 / || v - (e1*v)e1 ||
    e2 = v;
    vec3 e1_dot_v_e1;
    {
        const double e1_dot_v = vec3_dot(e1, v);
        e1_dot_v_e1 = e1;
        e1_dot_v_e1 *= e1_dot_v;
    }
    e2 -= e1_dot_v_e1;
    vec3_normalize(e2);

    // e0 = e1 x e2
    vec3_cross(e0, e1, e2);

    vec3_normalize(e0);

    Q[0][0] = e0[0]; Q[0][1] = e1[0]; Q[0][2] = e2[0];
    Q[1][0] = e0[1]; Q[1][1] = e1[1]; Q[1][2] = e2[1];
    Q[2][0] = e0[2]; Q[2][1] = e1[2]; Q[2][2] = e2[2];
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
    // mat3x3_mul:
    //
    //         | cosa  -sina*cosb  -sina*sinb |
    //     Q x | sina   cosa*cosb   cosa*sinb |
    //         | 0        -sinb       cosb   |


    mat3x3 A;

    A[0][0] = cosa; A[0][1] = -sina*cosb; A[0][2] = -sina*sinb;
    A[1][0] = sina; A[1][1] =  cosa*cosb; A[1][2] = cosa*sinb;
    A[2][0] = 0;    A[2][1] = -sinb;      A[2][2] = cosb;

    mat3x3_mul(Q_out, Q, A);
}

// BIdirectional SLERP
// Linearly interpolate two orhtogonal matrices Qa and Qb to produce a new
// orthogonal matrix Qab, which is determined by the interpolation factor t.
// When t = 0, Qab = Qa, and when t = 1, Qab = Qb.
//
// As defined in Function 4 in the supplementary material of Bayer2012.
MFEM_HOST_DEVICE
void bislerp(mat3x3& Qab, mat3x3& Qa, mat3x3& Qb, double t)
{
    // Translate the rotation matrices Qa and Qb into quaternions
    quat qa = {0}, qb = {0};
    rot2quat(qa, Qa);
    rot2quat(qb, Qb);
    quat_normalize(qa);
    quat_normalize(qb);

    // Find qm in { ±qa, ±i*qa, ±j*qa, ±k*qa} that maximizes ||qm*qb||
    const double a = qa[0], b = qa[1], c = qa[2], d = qa[3];

    quat i_qa = {-b,  a, -d,  c};
    quat j_qa = {-c,  d,  a, -b};
    quat k_qa = {-d, -c,  b,  a};

    quat minus_qa   = {-a, -b, -c, -d};
    quat minus_i_qa = {-i_qa[0], -i_qa[1], -i_qa[2], -i_qa[3]};
    quat minus_j_qa = {-j_qa[0], -j_qa[1], -j_qa[2], -j_qa[3]};
    quat minus_k_qa = {-k_qa[0], -k_qa[1], -k_qa[2], -k_qa[3]};

    quat *quat_array[8] = {
        &qa,   &minus_qa,
        &i_qa, &minus_i_qa,
        &j_qa, &minus_j_qa,
        &k_qa, &minus_k_qa,
    };

    quat qm = {0};
    double max_abs_dot  = -1.0;
    for (int i = 0; i < 8; i++) {
        quat *v = quat_array[i];
        const double abs_dot = abs(quat_dot(*v, qb));
        if (abs_dot > max_abs_dot) {
            max_abs_dot = abs_dot;
            qm = *v;
        }
    }

#if 1
    if (max_abs_dot > 1-1e-12) {
        quat2rot(Qab, qb);
        return;
    }
#endif

    // We have found the candiate qm that requires the smallest rotation angle.
    // Interpolate with slerp.
    quat q = {0};
    slerp(q, qm, qb, t);
    quat_normalize(q);
    quat2rot(Qab, q);
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
    const double tol_lv,
    const double tol_rv,
    const double tol_epi)
{
    MFEM_FORALL(i, n,
    {

        const double epi = CLAMP(phi_epi[i], 0.0, 1.0);
        const double lv  = CLAMP(phi_lv[i],  0.0, 1.0);
        const double rv  = CLAMP(phi_rv[i],  0.0, 1.0);

        vec3 grad_epi, grad_lv, grad_rv, grad_ab;
        vec3_set_from_ptr(grad_epi, &grad_phi_epi[3*i]);
        vec3_set_from_ptr(grad_lv,  &grad_phi_lv[3*i]);
        vec3_set_from_ptr(grad_rv,  &grad_phi_rv[3*i]);
        vec3_set_from_ptr(grad_ab,  &grad_psi_ab[3*i]);

        double depth = (lv > 0 || rv > 0) ? rv / (lv + rv) : 0.0;
        const double alpha_s = alpha_endo*(1.0-depth) - alpha_endo*depth;
        const double alpha_w = alpha_endo*(1.0-epi)   + alpha_epi*epi;

        const double beta_s  = beta_endo*(1.0-depth) - beta_endo*depth;
        const double beta_w  = beta_endo*(1.0-epi)   + beta_epi*epi;

        const double grad_epi_mag = vec3_magnitude(grad_epi);
        const double grad_lv_mag  = vec3_magnitude(grad_lv);
        const double grad_rv_mag  = vec3_magnitude(grad_rv);

        mat3x3 Q_lv = {{{0}}};
        {
            vec3 grad_lv_neg = grad_lv;
            vec3_negate(grad_lv_neg);

            mat3x3 T = {{{0}}};
            axis(T, grad_ab, grad_lv_neg);
            orient(Q_lv, T, alpha_s, beta_s);
        }


        mat3x3 Q_epi = {{{0}}};
        {
            mat3x3 T = {{{0}}};
            axis(T, grad_ab, grad_epi);
            orient(Q_epi, T, alpha_w, beta_w);
        }


        mat3x3 Q_fiber = {{{0}}};
        if (grad_epi_mag < tol_epi) {
            // We are in the septum
            Q_fiber = Q_lv;
        } else if (grad_epi_mag >= tol_epi && grad_lv_mag >= tol_lv && grad_rv_mag < tol_rv) {
            // We are in the LV free wall
            Q_fiber = Q_epi;
        } else if (grad_epi_mag >= tol_epi && grad_lv_mag < tol_lv  && grad_rv_mag >= tol_rv) {
            // We are in the RV free wall
            Q_fiber = Q_epi;
        } else {
            // We are in the junction between the septum and the LV and RV free
            // walls. Here we have to apply the full algorithm.
            mat3x3 Q_rv = {{{0}}};
            {
                mat3x3 T = {{{0}}};
                axis(T, grad_ab, grad_rv);
                orient(Q_rv, T, alpha_s, beta_s);
            }
            mat3x3 Q_endo = {{{0}}};
            bislerp(Q_endo, Q_lv, Q_rv, depth);
            bislerp(Q_fiber, Q_endo, Q_epi, epi);
        }

        F[3*i+0] = Q_fiber[0][0]; S[3*i+0] = Q_fiber[0][1]; T[3*i+0] = Q_fiber[0][2];
        F[3*i+1] = Q_fiber[1][0]; S[3*i+1] = Q_fiber[1][1]; T[3*i+1] = Q_fiber[1][2];
        F[3*i+2] = Q_fiber[2][0]; S[3*i+2] = Q_fiber[2][1]; T[3*i+2] = Q_fiber[2][2];

    });
    cf_device_barrier();
}

void define_fibers_single_ventricle(
        int n,
        const double *phi_epi,
        const double *grad_phi_epi,
        const double *grad_psi_ab,
        double alpha_endo,
        double alpha_epi,
        double beta_endo,
        double beta_epi,
        double *F,
        double *S,
        double *T)
{
    MFEM_FORALL(i, n, {
        // In the case where there is a single ventricle, we can simply skip
        // all the considerations we have to make for biventricular meshes. Use
        // the gradient `grad_psi_ab` as the apicobasal, as usual, and set the
        // transmural vector to be `grad_phi_epi`. Orient using alpha_w and
        // beta_w.

        double epi = phi_epi[i];
        vec3 grad_epi, grad_ab;
        vec3_set_from_ptr(grad_epi, &grad_phi_epi[3*i]);
        vec3_set_from_ptr(grad_ab,  &grad_psi_ab[3*i]);
        const double alpha_w = alpha_endo*(1.0-epi) + alpha_epi*epi;
        const double beta_w  = beta_endo*(1.0-epi)  + beta_epi*epi;

        mat3x3 Q_epi = {{{0}}};
        {
            mat3x3 T = {{{0}}};
            axis(T, grad_ab, grad_epi);
            orient(Q_epi, T, alpha_w, beta_w);
        }

        F[3*i+0] = Q_epi[0][0]; S[3*i+0] = Q_epi[0][1]; T[3*i+0] = Q_epi[0][2];
        F[3*i+1] = Q_epi[1][0]; S[3*i+1] = Q_epi[1][1]; T[3*i+1] = Q_epi[1][2];
        F[3*i+2] = Q_epi[2][0]; S[3*i+2] = Q_epi[2][1]; T[3*i+2] = Q_epi[2][2];

    });
    cf_device_barrier();
}


