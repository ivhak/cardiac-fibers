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
#include "calculus_gpu.hpp"


#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

using namespace mfem;

// Cross product of two 3D vectors a and b, store in c.
MFEM_HOST_DEVICE static void cross(Vector3D& c, Vector3D& a, Vector3D& b)
{
    // c = a x b
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

// Dot product of two quaternions
MFEM_HOST_DEVICE static double quatdot(Quaternion& q1, Quaternion& q2)
{
    return q1[0] * q2[0]
         + q1[1] * q2[1]
         + q1[2] * q2[2]
         + q1[3] * q2[3];
}

// Set q1 = q2
MFEM_HOST_DEVICE static void quatcopy(Quaternion& q1, Quaternion& q2)
{
    q1[0] = q2[0];
    q1[1] = q2[1];
    q1[2] = q2[2];
    q1[3] = q2[3];
}

// Negate the values of quaternion q
MFEM_HOST_DEVICE static void quatnegate(Quaternion& q)
{
    q[0] = -(q[0]);
    q[1] = -(q[1]);
    q[2] = -(q[2]);
    q[3] = -(q[3]);
}

// Multiple quaternion q with scalar a
MFEM_HOST_DEVICE static void quatmul(Quaternion& q, double a)
{
    q[0] *= a;
    q[1] *= a;
    q[2] *= a;
    q[3] *= a;
}

// Elementwise multiplication q2 into q1
MFEM_HOST_DEVICE static void quatelemmul(Quaternion& q1, Quaternion& q2)
{
    q1[0] *= q2[0];
    q1[1] *= q2[1];
    q1[2] *= q2[2];
    q1[3] *= q2[3];
}

// Add quaternion q2 to q1
MFEM_HOST_DEVICE static void quatadd(Quaternion& q1, Quaternion& q2)
{
    q1[0] += q2[0];
    q1[1] += q2[1];
    q1[2] += q2[2];
    q1[3] += q2[3];
}

// Normalize quaternion
MFEM_HOST_DEVICE static void quatnormalize(Quaternion& q)
{
    double sum = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    double sum_sqr = sqrt(sum);
    if (sum_sqr == 0.0) return;
    q[0] = q[0] / sum_sqr;
    q[1] = q[1] / sum_sqr;
    q[2] = q[2] / sum_sqr;
    q[3] = q[3] / sum_sqr;
}

// Set vector a = b
MFEM_HOST_DEVICE void veccopy(Vector3D& a, Vector3D& b)
{
    a[0] = b[0];
    a[1] = b[1];
    a[2] = b[2];
}

// Normalize vector a, a = a / ||a||
MFEM_HOST_DEVICE static void vecnormalize(Vector3D& a)
{
    const double sum = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    const double sum_sqr = sqrt(sum);
    if (sum == 0.0) return;
    a[0] = a[0] / sum_sqr;
    a[1] = a[1] / sum_sqr;
    a[2] = a[2] / sum_sqr;
}

// Dot product of vectors a and b
MFEM_HOST_DEVICE double vecdot(Vector3D& a, Vector3D& b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Subtract vector b from a; a = a - b
MFEM_HOST_DEVICE static void vecsub(Vector3D& a, Vector3D& b)
{
    a[0] = a[0] - b[0];
    a[1] = a[1] - b[1];
    a[2] = a[2] - b[2];
}

// Scale vector a; a_i = a_i * b
MFEM_HOST_DEVICE void vecmul(Vector3D& a, double b)
{
    a[0] = a[0] * b;
    a[1] = a[1] * b;
    a[2] = a[2] * b;
}

// Negate vector a; a = -a
MFEM_HOST_DEVICE static void vecnegate(Vector3D& a)
{
    a[0] = -a[0];
    a[1] = -a[1];
    a[2] = -a[2];
}

// Initialize matrix A values to zero
MFEM_HOST_DEVICE static void matzero(Matrix3x3& A)
{
    A[0][0] = 0.0; A[0][1] = 0.0; A[0][2] = 0.0;
    A[1][0] = 0.0; A[1][1] = 0.0; A[1][2] = 0.0;
    A[2][0] = 0.0; A[2][1] = 0.0; A[2][2] = 0.0;
}

// Copy matrix B into A
MFEM_HOST_DEVICE static void matcopy(Matrix3x3& A, Matrix3x3& B)
{
    A[0][0] = B[0][0]; A[0][1] = B[0][1]; A[0][2] = B[0][2];
    A[1][0] = B[1][0]; A[1][1] = B[1][1]; A[1][2] = B[1][2];
    A[2][0] = B[2][0]; A[2][1] = B[2][1]; A[2][2] = B[2][2];
}

// Set matrix C to the matrix multiplication of A and B; C = A x B
MFEM_HOST_DEVICE static void matmul(Matrix3x3& C, Matrix3x3& A, Matrix3x3& B )
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
MFEM_HOST_DEVICE void rot2quat(Quaternion& q, Matrix3x3& Q)
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
MFEM_HOST_DEVICE void quat2rot(Matrix3x3& Q, Quaternion& q)
{
    const double w = q[0], x = q[1], y = q[2], z = q[3];

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
MFEM_HOST_DEVICE static void slerp(Quaternion q, Quaternion q1, Quaternion q2, double t)
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

        Quaternion q1a = {0};
        quatcopy(q1a, q1);
        quatmul(q1a, a);
        quatmul(q, b);
        quatadd(q, q1a);
    } else {
        // Linear interpolation: q = q1(1-t) + q2*t

        quatcopy(q, q1);
        quatmul(q, (1-t));

        Quaternion qt;
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
MFEM_HOST_DEVICE void axis(Matrix3x3& Q, Vector3D& u, Vector3D& v)
{

    Vector3D e0, e1, e2;

    // e1 = u / ||u||
    veccopy(e1, u);
    vecnormalize(e1);

    // e2 = u / || u ||
    // where u = v - (e0*v)e0
    //
    // Normalize v as an initial guess for e0
    veccopy(e2, v);
    vecnormalize(e2);


    Vector3D e1_dot_e2_e1;
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
MFEM_HOST_DEVICE void orient(Matrix3x3& Q_out, Matrix3x3& Q, double a, double b)
{


    const double sina = sin(a*PI/180);
    const double sinb = sin(b*PI/180);
    const double cosa = cos(a*PI/180);
    const double cosb = cos(b*PI/180);

    Matrix3x3 A;

    A[0][0] = cosa;
    A[0][1] = -sina*cosb;
    A[0][2] = -sina*sinb;

    A[1][0] = sina;
    A[1][1] = cosa*cosb;
    A[1][2] = cosa*sinb;

    A[2][0] = 0;
    A[2][1] = -sinb;
    A[2][2] = cosb;

    matmul(Q_out, Q, A);
}
//
// BIdirectional SLERP
// Linearly interpolate two orhtogonal matrices Qa and Qb to produce a new
// orthogonal matrix Qab, which is determined by the interpolation factor t.
// When t = 0, Qab = Qa, and when t = 1, Qab = Qb.
//
// As defined in Function 4 in the supplementary material of Bayer2012.
MFEM_HOST_DEVICE void bislerp(Matrix3x3& Qab, Matrix3x3& Qa, Matrix3x3& Qb, double t)
{
    const double tol = 1e-12;

    if (t <= tol) {
        matcopy(Qab, Qa);
        return;
    }

    if (t >= 1.0 - tol) {
        matcopy(Qab, Qb);
        return;
    }

    // Translate the rotation matrices Qa and Qb into quaternions
    Quaternion qa, qb;
    rot2quat(qa, Qa);
    rot2quat(qb, Qb);

    // Find qm in { ±qa, ±i*qa, ±j*qa, ±k*qa} that maximizes ||qm*qb||

    Quaternion i; i[0] = 0.0; i[1] = 1.0; i[2] = 0.0; i[3] = 0.0;
    Quaternion j; j[0] = 0.0; j[1] = 0.0; j[2] = 1.0; j[3] = 0.0;
    Quaternion k; k[0] = 0.0; k[1] = 0.0; k[2] = 0.0; k[3] = 1.0;

    Quaternion i_qa; quatcopy(i_qa, i); quatelemmul(i_qa, qa);
    Quaternion j_qa; quatcopy(j_qa, j); quatelemmul(j_qa, qa);
    Quaternion k_qa; quatcopy(k_qa, k); quatelemmul(k_qa, qa);

    Quaternion qa_minus;   quatcopy(qa_minus,   qa);   quatmul(qa_minus, -1.0);
    Quaternion i_qa_minus; quatcopy(i_qa_minus, i_qa); quatmul(i_qa_minus, -1.0);
    Quaternion j_qa_minus; quatcopy(j_qa_minus, j_qa); quatmul(j_qa_minus, -1.0);
    Quaternion k_qa_minus; quatcopy(k_qa_minus, k_qa); quatmul(k_qa_minus, -1.0);

    Quaternion* quat_array[8];
    quat_array[0] = &qa;
    quat_array[1] = &qa_minus;
    quat_array[2] = &i_qa;
    quat_array[3] = &i_qa_minus;
    quat_array[4] = &j_qa;
    quat_array[5] = &j_qa_minus;
    quat_array[6] = &k_qa;
    quat_array[7] = &k_qa_minus;

    double max_abs_dot  = 0.0;
    Quaternion qm;
    for (int i = 0; i < 8; i++) {
        Quaternion *v = quat_array[i];
        const double abs_dot = abs(quatdot(qb, (*v)));
        if (abs_dot > max_abs_dot) {
            max_abs_dot = abs_dot;
            quatcopy(qm, (*v)) ;
        }
    }

#if 1
    // If the angle is very small, i.e. max_dot is very close to one, return Qb.
    if (max_abs_dot > 1-tol) {
        matcopy(Qab, Qb);
        return;
    }
#endif

    Quaternion q;
    slerp(q, qm, qb, t);
    quatnormalize(q);
    quat2rot(Qab, q);
}

void define_fibers(
    int n,
    const double *phi_epi,
    const double *phi_lv,
    const double *phi_rv,
    const double *psi_ab,
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
    double *T)
{

    const double tol = 1e-12;
    MFEM_FORALL(i, n,
    {

        const double phi_epi_i = CLAMP(phi_epi[i], 0.0, 1.0);
        const double phi_lv_i  = CLAMP(phi_lv[i],  0.0, 1.0);
        const double phi_rv_i  = CLAMP(phi_rv[i],  0.0, 1.0);

        Vector3D grad_phi_epi_i;
        grad_phi_epi_i[0] = grad_phi_epi[3*i+0];
        grad_phi_epi_i[1] = grad_phi_epi[3*i+1];
        grad_phi_epi_i[2] = grad_phi_epi[3*i+2];

        Vector3D grad_phi_lv_i;
        grad_phi_lv_i[0] = grad_phi_lv[3*i+0];
        grad_phi_lv_i[1] = grad_phi_lv[3*i+1];
        grad_phi_lv_i[2] = grad_phi_lv[3*i+2];

        Vector3D grad_phi_rv_i;
        grad_phi_rv_i[0] = grad_phi_rv[3*i+0];
        grad_phi_rv_i[1] = grad_phi_rv[3*i+1];
        grad_phi_rv_i[2] = grad_phi_rv[3*i+2];

        Vector3D grad_psi_ab_i;
        grad_psi_ab_i[0] = grad_psi_ab[3*i+0];
        grad_psi_ab_i[1] = grad_psi_ab[3*i+1];
        grad_psi_ab_i[2] = grad_psi_ab[3*i+2];

        // TODO: What to do here? TOLERANCE
        double depth;
        if (phi_lv_i + phi_rv_i < tol) {
            depth = 0.5;
        } else {
            depth = phi_rv_i / (phi_lv_i + phi_rv_i);
        }

        const double alpha_s_d = alpha_endo*(1.0-depth) - alpha_endo*depth;
        const double beta_s_d  = beta_endo *(1.0-depth) - beta_endo *depth;

        const double alpha_w_epi = alpha_endo*(1.0-phi_epi_i) + alpha_epi * phi_epi_i;
        const double beta_w_epi  = beta_endo *(1.0-phi_epi_i) + beta_epi  * phi_epi_i;


        Matrix3x3 Q_lv;
        matzero(Q_lv);
        if (phi_lv_i > tol) {
            Matrix3x3 T;
            Vector3D grad_phi_lv_i_neg;
            {
                veccopy(grad_phi_lv_i_neg, grad_phi_lv_i);
                vecnegate(grad_phi_lv_i_neg);
            }
            axis(T, grad_psi_ab_i, grad_phi_lv_i_neg);
            orient(Q_lv, T, alpha_s_d, beta_s_d);
        }

        Matrix3x3 Q_rv;
        matzero(Q_rv);
        if (phi_rv_i > tol) {
            Matrix3x3 T;
            axis(T, grad_psi_ab_i, grad_phi_rv_i);
            orient(Q_rv, T, alpha_s_d, beta_s_d);
        }

        Matrix3x3 Q_endo;
        bislerp(Q_endo, Q_lv, Q_rv, depth);

        Matrix3x3 Q_epi;
        matzero(Q_epi);
        if (phi_epi_i > tol) {
            Matrix3x3 T;
            axis(T, grad_psi_ab_i, grad_phi_epi_i);
            orient(Q_epi, T, alpha_w_epi, beta_w_epi);
        }

        Matrix3x3 FST;
        bislerp(FST, Q_endo, Q_epi, phi_epi_i);

        F[3*i+0] = FST[0][0];
        F[3*i+1] = FST[1][0];
        F[3*i+2] = FST[2][0];

        S[3*i+0] = FST[0][1];
        S[3*i+1] = FST[1][1];
        S[3*i+2] = FST[2][1];

        T[3*i+0] = FST[0][2];
        T[3*i+1] = FST[1][2];
        T[3*i+2] = FST[2][2];
    });
}
// Set the gradient in each vertex to be the average of the gradient in the
// centers of the surrounding elements
void par_calculate_gradients(double* grads, ParGridFunction& x, ParMesh& mesh, Table* v2e)
{

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

            ElementTransformation *tr = x.ParFESpace()->GetElementTransformation(el_id);
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
        grads[3*i + 0] = vertex_gradient(0);
        grads[3*i + 1] = vertex_gradient(1);
        grads[3*i + 2] = vertex_gradient(2);
    }
}
