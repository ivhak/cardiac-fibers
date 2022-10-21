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
#include "calculus.hpp"

#define CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

using namespace mfem;

// Cross product of two 3D vectors a and b, store in c.
static void cross(Vector& c, Vector &a, Vector &b) {
    MFEM_ASSERT(a.Size() == 3, "a is of the wrong size, should be 3.");
    MFEM_ASSERT(b.Size() == 3, "b is of the wrong size, should be 3.");
    MFEM_ASSERT(c.Size() == 3, "c is of the wrong size, should be 3.");

    // c = a x b
    c(0) = a(1)*b(2) - a(2)*b(1);
    c(1) = a(2)*b(0) - a(0)*b(2);
    c(2) = a(0)*b(1) - a(1)*b(0);
}

// Convert the 3x3 rotation matrix Q to a quaternion q, using the algorithm
// described in Appendix 1 of
//
//   Shoemake, K. (1985, July). Animating rotation with quaternion curves. In
//   Proceedings of the 12th annual conference on Computer graphics and
//   interactive techniques (pp. 245-254).
void rot2quat(Vector& q, DenseMatrix& Q)
{
    MFEM_ASSERT(Q.Width() == 3 && Q.Height() == 3, "Q is of the wrong size, should be 3x3.");
    MFEM_ASSERT(q.Size() == 4, "q is the wrong size, should be 4.");

    const double M11=Q(0,0), M12=Q(1,0), M13=Q(2,0);
    const double M21=Q(0,1), M22=Q(1,1), M23=Q(2,1);
    const double M31=Q(0,2), M32=Q(1,2), M33=Q(2,2);

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

    q(0) = w;
    q(1) = x;
    q(2) = y;
    q(3) = z;

    // Make sure the quaternion is normalized
    q /= q.Norml2();
}

// Convert the quaternion q to a 3x3 rotation matrix Q, using the algorithm
// described in Appendix 1 of
//
//   Shoemake, K. (1985, July). Animating rotation with quaternion curves. In
//   Proceedings of the 12th annual conference on Computer graphics and
//   interactive techniques (pp. 245-254).
void quat2rot(DenseMatrix& Q, Vector& q)
{
    MFEM_ASSERT(q.Size() == 4, "q is of the wrong size, should be 4.");
    MFEM_ASSERT(Q.Width() == 3 && Q.Height() == 3, "Q is of the wrong size, should be 3x3.");

    const double w = q(0), x = q(1), y = q(2), z = q(3);

    {
        const double dot = w*w + x*x + y*y + z*z;
        MFEM_ASSERT(dot - 1.0 <= 1e-12, "Quaternion not normalized");
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

    Q(0,0) = 1.0 - 2.0*y2 - 2.0*z2;
    Q(1,0) =       2.0*xy + 2.0*wz;
    Q(2,0) =       2.0*xz - 2.0*wy;

    Q(0,1) =       2.0*xy - 2.0*wz;
    Q(1,1) = 1.0 - 2.0*x2 - 2.0*z2;
    Q(2,1) =       2.0*yz + 2.0*wx;

    Q(0,2) =       2.0*xz + 2.0*wy;
    Q(1,2) =       2.0*yz - 2.0*wx;
    Q(2,2) = 1.0 - 2.0*x2 - 2.0*y2;

}

// Spherical Linear intERPolation
static void slerp(Vector& q, Vector& q1, Vector &q2, double t)
{
    double dot = q1 * q2;
    q = q2;

    if (dot < 0) {
        dot = -dot;
        q.Neg();
    }

    // XXX: The LLNL/cardiod implementation does this. Verify why!
    if (dot < 1.0 - 1e-4) {
        // Slerp(q1, q2, t) = ((sin(1-t)*theta)/sin(theta))q1 + ((sin(t)*theta)/sin(theta))q2
        // where theta = acos(q1 dot q2)
        const double angle = acos(dot);
        const double a = sin(angle * (1-t))/sin(angle);
        const double b = sin(angle * t)/sin(angle);

        Vector q1a = q1;
        q1a *= a;
        q *= b;
        q += q1a;
    } else {
        // Use linear interpolation: q = q1(1-t) + q2*t
        Vector a = q1;
        q1 *= (1-t);

        Vector b = q2;
        b *= t;

        q1 += b;
        q = q1;
    }

}

// Given two vectors u and v, u being a vector pointing in
// the apicobasal direction and v a vector pointing in the transmural
// direction, return a 3x3 orthogonal matrix representing the coordinate system
// for assigning fiber orientation within the myocardium.
//
// As defined in Function 2 in the supplementary material of Bayer2012.
void axis(DenseMatrix& Q, Vector& u, Vector &v)
{

    MFEM_ASSERT(u.Size() == 3, "u is the wrong size");
    MFEM_ASSERT(v.Size() == 3, "v is the wrong size");
    MFEM_ASSERT(Q.Width() == 3 && Q.Height() == 3, "Q is of the wrong size, should be 3x3.");

    Vector e_0(3), e_1(3), e_2(3);

    // e_1 = u / ||u||
    e_1 = u;
    e_1 /= u.Norml2();

    // e_2 = u / || u ||
    // where u = v - (e_0*v)e_0
    //
    // Normalize v as an initial guess for e_0
    e_2 = v;
    e_2 /= v.Norml2();

    Vector e1_dot_e2_e1(3);
    {
        const double e1_dot_e2 = e_1 * e_2;
        e1_dot_e2_e1 = e_1;
        e1_dot_e2_e1 *= e1_dot_e2;
    }
    e_2 -= e1_dot_e2_e1;
    e_2 /= e_2.Norml2();


    // e_0 = e_1 x e_2
    cross(e_0, e_1, e_2);

    e_0 /= e_0.Norml2();

    Q.SetCol(0, e_0); // Circumferential direction
    Q.SetCol(1, e_1); // Apicobasal direction
    Q.SetCol(2, e_2); // Transmural direction
}

// Take the coordinate system Q, in the form of a 3x3 matrix, and the fiber
// orientation angles a(lpha) and b(eta) at a given point in the mesh, and
// return an orthonormal coordinate system (F S T), in the form of a 3x3
// matrix, where F is the longitudinal direction, S is the sheet normal, and T
// is the transverse direction.
//
// As defined in Function 3 in the supplementary material of Bayer2012.
void orient(DenseMatrix& Q_out, DenseMatrix& Q, double a, double b)
{
    MFEM_ASSERT(Q_out.Width() == 3 && Q_out.Height() == 3, "Q_out is the wrong size, should be 3x3.");
    MFEM_ASSERT(Q.Width() == 3 && Q.Width() == 3,          "Q is the wrong size, should be 3x3.");


    const double sina = sin(a*PI/180);
    const double sinb = sin(b*PI/180);
    const double cosa = cos(a*PI/180);
    const double cosb = cos(b*PI/180);

    const double A_vals[3][3] = {
        {cosa, -sina, 0},
        {sina,  cosa, 0},
        {0,     0,    1}
    };

    const double B_vals[3][3] = {
        {1,  0,    0   },
        {0,  cosb, sinb},
        {0, -sinb, cosb}
    };

    DenseMatrix A(A_vals);
    DenseMatrix B(B_vals);

    DenseMatrix Temp(3, 3);

    Mult(Q, A, Temp);
    Mult(Temp, B, Q_out);
}

// BIdirectional SLERP
// Linearly interpolate two orhtogonal matrices Qa and Qb to produce a new
// orthogonal matrix Qab, which is determined by the interpolation factor t.
// When t = 0, Qab = Qa, and when t = 1, Qab = Qb.
//
// As defined in Function 4 in the supplementary material of Bayer2012.
void bislerp(DenseMatrix& Qab, DenseMatrix& Qa, DenseMatrix& Qb, double t)
{
    const double tol = 1e-12;

    MFEM_ASSERT(t >= 0.0 && t <= 1.0, "t is not in [0, 1].");

    if (t <= tol) {
        Qab = Qa;
        return;
    }

    if (t >= 1.0 - tol) {
        Qab = Qb;
        return;
    }

    // Translate the rotation matrices Qa and Qb into quaternions
    Vector qa(4), qb(4);
    rot2quat(qa, Qa);
    rot2quat(qb, Qb);

    // Find qm in { ±qa, ±i*qa, ±j*qa, ±k*qa} that maximizes ||qm*qb||

    Vector i(4); i = 0.0; i(1) = 1.0;
    Vector j(4); j = 0.0; j(2) = 1.0;
    Vector k(4); k = 0.0; k(3) = 1.0;

    Vector i_qa = i; i_qa *= qa;
    Vector j_qa = j; j_qa *= qa;
    Vector k_qa = k; k_qa *= qa;

    Vector qa_minus   =   qa;   qa_minus *= -1.0;
    Vector i_qa_minus = i_qa; i_qa_minus *= -1.0;
    Vector j_qa_minus = j_qa; j_qa_minus *= -1.0;
    Vector k_qa_minus = k_qa; k_qa_minus *= -1.0;

    Array<Vector*> quat_array(8);
    quat_array[0] = &qa;
    quat_array[1] = &qa_minus;
    quat_array[2] = &i_qa;
    quat_array[3] = &i_qa_minus;
    quat_array[4] = &j_qa;
    quat_array[5] = &j_qa_minus;
    quat_array[6] = &k_qa;
    quat_array[7] = &k_qa_minus;

    double max_abs_dot  = 0.0;
    Vector qm(4);
    for (int i = 0; i < 8; i++) {
        Vector v = *quat_array[i];
        const double abs_dot = abs(qb * v);
        if (abs_dot > max_abs_dot) {
            max_abs_dot = abs_dot;
            qm = v;
        }
    }

#if 1
    // If the angle is very small, i.e. max_dot is very close to one, return Qb.
    if (max_abs_dot > 1-tol) {
        Qab = Qb;
        return;
    }
#endif

    Vector q(4);
    slerp(q, qm, qb, t);
    q /= q.Norml2();
    quat2rot(Qab, q);
}


// Set the gradient in each vertex to be the average of the gradient in the
// centers of the surrounding elements
void calculate_gradients(double* grads, GridFunction& x, Mesh& mesh, Table* v2e)
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
        grads[3*i + 0] = vertex_gradient(0);
        grads[3*i + 1] = vertex_gradient(1);
        grads[3*i + 2] = vertex_gradient(2);
    }
}

#if 0
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
#endif

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
    for (int i = 0; i < n; i++) {

        const double phi_epi_i = CLAMP(phi_epi[i], 0.0, 1.0);
        const double phi_lv_i  = CLAMP(phi_lv[i],  0.0, 1.0);
        const double phi_rv_i  = CLAMP(phi_rv[i],  0.0, 1.0);

        MFEM_ASSERT(abs(phi_epi_i + phi_lv_i + phi_rv_i - 1.0) < 1e-3,
                    "The laplacians epi, lv and rv do not add up to 1.");

        Vector grad_phi_epi_i((double *)&grad_phi_epi[3*i], 3);
        Vector grad_phi_lv_i((double *)&grad_phi_lv[3*i], 3);
        Vector grad_phi_rv_i((double *)&grad_phi_rv[3*i], 3);
        Vector grad_psi_ab_i((double *)&grad_psi_ab[3*i], 3);

        MFEM_ASSERT(abs(grad_phi_epi_i(0) + grad_phi_lv_i(0) + grad_phi_rv_i(0)) < 1e-3
                 && abs(grad_phi_epi_i(1) + grad_phi_lv_i(1) + grad_phi_rv_i(1)) < 1e-3
                 && abs(grad_phi_epi_i(2) + grad_phi_lv_i(2) + grad_phi_rv_i(2)) < 1e-3,
                    "The gradients do not add up to zero");

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


        DenseMatrix Q_lv(3,3);
        Q_lv = 0.0;
        if (phi_lv_i > tol) {
            DenseMatrix T(3,3);
            Vector grad_phi_lv_i_neg(3);
            {
                grad_phi_lv_i_neg = grad_phi_lv_i;
                grad_phi_lv_i_neg.Neg();
            }
            axis(T, grad_psi_ab_i, grad_phi_lv_i_neg);
            orient(Q_lv, T, alpha_s_d, beta_s_d);
        }

        DenseMatrix Q_rv(3,3);
        Q_rv = 0.0;
        if (phi_rv_i > tol) {
            DenseMatrix T(3,3);
            axis(T, grad_psi_ab_i, grad_phi_rv_i);
            orient(Q_rv, T, alpha_s_d, beta_s_d);
        }

        DenseMatrix Q_endo(3,3);
        bislerp(Q_endo, Q_lv, Q_rv, depth);

        DenseMatrix Q_epi(3,3);
        Q_epi = 0.0;
        if (phi_epi_i > tol) {
            DenseMatrix T(3,3);
            axis(T, grad_psi_ab_i, grad_phi_epi_i);
            orient(Q_epi, T, alpha_w_epi, beta_w_epi);
        }

        DenseMatrix FST(3,3);
        bislerp(FST, Q_endo, Q_epi, phi_epi_i);

        F[3*i+0] = FST(0,0);
        F[3*i+1] = FST(1,0);
        F[3*i+2] = FST(2,0);

        S[3*i+0] = FST(0,1);
        S[3*i+1] = FST(1,1);
        S[3*i+2] = FST(2,1);

        T[3*i+0] = FST(0,2);
        T[3*i+1] = FST(1,2);
        T[3*i+2] = FST(2,2);
    }
}
