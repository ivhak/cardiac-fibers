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

using namespace mfem;

// Cross product of two 3D vectors b and c, store in a.
static void cross(Vector& a, Vector &b, Vector &c) {
    MFEM_ASSERT(a.Size() == 3, "a is of the wrong size, should be 3.");
    MFEM_ASSERT(b.Size() == 3, "b is of the wrong size, should be 3.");
    MFEM_ASSERT(c.Size() == 3, "c is of the wrong size, should be 3.");

    // a = b x c
    a(0) = b(1)*c(2) - a(2)*b(1);
    a(1) = b(2)*c(0) - a(0)*b(2);
    a(2) = b(0)*c(1) - a(1)*b(0);
}

// Convert the 3x3 rotation matrix Q to a quaternion q, using the algorithm
// described in Appendix 1 of
//
//   Shoemake, K. (1985, July). Animating rotation with quaternion curves. In
//   Proceedings of the 12th annual conference on Computer graphics and
//   interactive techniques (pp. 245-254).
static void rot2quat(Vector& q, DenseMatrix& Q)
{
    // MFEM_ASSERT(Q.Size() == 3, "Q is of the wrong size, should be 3x3.");

    q.SetSize(4); // q = w + x*i + y*j + z*k

    double M11, M21, M31;
    double M12, M22, M32;
    double M13, M23, M33;

    M11=Q(0,0); M21=Q(1,0); M31=Q(2,0);
    M12=Q(0,1); M22=Q(1,1); M32=Q(2,1);
    M13=Q(0,2); M23=Q(1,2); M33=Q(2,2);

    double w2 = 0.25 * (1 + M11 + M22 + M33);
    double err = 1e-6;

    double w, x, y, z;

    if (w2 > err) {
        w = sqrt(w2);
        x = (M23 - M32) / (4.0*w);
        y = (M31 - M13) / (4.0*w);
        z = (M12 - M21) / (4.0*w);
    } else {
        w = 0.0;
        double x2 = -0.5*(M22 + M33);
        if (x2 > err) {
            x = sqrt(x2);
            y = M12 / (2.0*x);
            z = M13 / (2.0*x);
        } else {
            x = 0.0;
            double y2 = 0.5*(1-M33);
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
}

// Convert the quaternion q to a 3x3 rotation matrix Q, using the algorithm
// described in Appendix 1 of
//
//   Shoemake, K. (1985, July). Animating rotation with quaternion curves. In
//   Proceedings of the 12th annual conference on Computer graphics and
//   interactive techniques (pp. 245-254).
static void quat2rot(DenseMatrix& Q, Vector& q)
{
    MFEM_ASSERT(q.Size() == 4, "q is of the wrong size, should be 4.");

    Q.SetSize(3, 3);

    double w, x, y, z;
    w = q(0); x = q(1); y = q(2); z = q(3);

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
    Q(1,1) = 1.0 - 2.0*x2 - 2.0*y2;
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

    // Slerp(q1, q2, t) = ((sin(1-t)*theta)/sin(theta))q1 + ((sin(t)*theta)/sin(theta))q2
    // where theta = acos(q1 dot q2)
    const double angle = acos(dot);
    const double a = sin(angle * (1-t))/sin(angle);
    const double b = sin(angle * t)/sin(angle);

    Vector q1a = q1;
    q1a *= a;
    q *= b;
    q += q1a;

}

// Given two vectors grad_psi and grad_phi, grad_psi being a vector pointing in
// the apicobasal direction and grad_phi a vector pointing in the transmural
// direction, return a 3x3 orthogonal matrix representing the coordinate system
// for assigning fiber orientation within the myocardium.
//
// As defined in Function 2 in the supplementary material of Bayer2012.
void axis(DenseMatrix& Q, Vector& grad_psi, Vector &grad_phi)
{

    MFEM_ASSERT(grad_psi.Size() == 3, "grad_psi is the wrong size");
    MFEM_ASSERT(grad_phi.Size() == 3, "grad_phi is the wrong size");

    // MFEM_ASSERT(Q.Size() == 3, "Q is of the wrong size, should be 3x3.");
    Vector e_0(3), e_1(3), e_2(3);

    // e_1 = grad_psi / ||grad_psi||
    e_1 = grad_psi;
    e_1 /= grad_psi.Norml2();

    // e_2 = u / || u ||
    //
    // where u = grad_phi - (e_0*grad_phi)e_0
    //
    // Normalize grad_phi as an initial guess for e_0
    e_0 = grad_phi;
    e_0 /= grad_phi.Norml2();

    // Calculate u = grad_phi - t*e_0
    // where t = (e_0 * grad_phi) is a scalar
    Vector u(3);;
    {
        Vector e_0_t;
        const double t = e_1 * e_2;
        e_0_t = e_0;
        e_0_t *= t;
        u = grad_phi - e_0_t;
    }

    e_2 = u;
    e_2 /= u.Norml2();

    // e_0 = e_1 x e_2
    cross(e_0, e_1, e_2);

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
    // MFEM_ASSERT(Q_out.Size() == 3, "Q_out is the wrong size, should be 3x3.")
    // MFEM_ASSERT(Q.Size() == 3,     "Q_out is the wrong size, should be 3x3.")


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

    // Translate the rotation matrices Qa and Qb into quaternions
    Vector qa, qb;
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

    double max_dot = 0.0;
    Vector qm(4);
    for (int i = 0; i < 8; i++) {
        Vector v = *quat_array[i];
        const double dot = abs(qb * v);
        if (dot < max_dot) {
            max_dot = dot;
            qm = v;
        }
    }

    // If the angle is very small, i.e. max_dot is very close to one, return Qb.
    if (max_dot > 1-tol) {
        Qab = Qb;
        return;
    }

    Vector q;
    slerp(q, qm, qb, t);
    quat2rot(Qab, q);
}


