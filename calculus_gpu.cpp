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

typedef struct {
    double data[2];
    double& operator[](int i) {
        return data[i];
    }
} Vector3D;

typedef struct {
    double data[4]
    double& operator[](int i) {
        return data[i];
    }
} Quaternion;

typedef struct {
    double data[3][3]
    double* operator[](int i) {
        return (double *) &data[i];
    }
} Matrix3x3;



// Cross product of two 3D vectors a and b, store in c.
__device__ static void cross(Vector3D& c, Vector3D& a, Vector3D& b)
{
    // c = a x b
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

// Dot product of two quaternions
__device__ static double quatdot(Quaternion& q1, Quaternion& q2)
{
    return q1[0] * q2[0]
         + q1[1] * q2[1]
         + q1[2] * q2[2]
         + q1[3] * q2[3];
}

// Set q1 = q2
__device__ static double quatcopy(Quaternion& q1, Quaternion& q2)
{
    q1[0] = q2[0];
    q1[1] = q2[1];
    q1[2] = q2[2];
    q1[3] = q2[3];
}

// Negate the values of quaternion q
__device__ static double quatnegate(Quaternion& q)
{
    q[0] = -(q[0]);
    q[1] = -(q[1]);
    q[2] = -(q[2]);
    q[3] = -(q[3]);
}

// Multiple quaternion q with scalar a
__device__ static double quatmult(Quaternion& q, double a)
{
    q[0] *= a;
    q[1] *= a;
    q[2] *= a;
    q[3] *= a;
}

// Add quaternion q2 to q1
__device__ static double quatadd(Quaternion& q1, Quaternion& q2)
{
    q1[0] += q2[0];
    q1[1] += q2[1];
    q1[2] += q2[2];
    q1[3] += q2[3];
}

__device__ static double quatnormalize(Quaternion& q)
{
    double sum = q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    double sum_sqr = sqrt(sum);
    q[0] = q[0] / sum_sqr;
    q[1] = q[1] / sum_sqr;
    q[2] = q[2] / sum_sqr;
    q[3] = q[3] / sum_sqr;
}

// Convert the 3x3 rotation matrix Q to a quaternion q, using the algorithm
// described in Appendix 1 of
//
//   Shoemake, K. (1985, July). Animating rotation with quaternion curves. In
//   Proceedings of the 12th annual conference on Computer graphics and
//   interactive techniques (pp. 245-254).
__device__ void rot2quat(Vector3D& q, Matrix3x3& Q)
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
__device__ static void quat2rot(Matrix3x3 Q, Quaternion q)
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
__device__ static void slerp(Quaternion q, Quaternion q1, Quaternion q2, double t)
{
    double dot = quatdot(q1, q2);
    quatcopy(q, q2);

    if (dot < 0) {
        dot = -dot;
        quatnegate(&q);
    }

    // Slerp(q1, q2, t) = ((sin(1-t)*theta)/sin(theta))q1 + ((sin(t)*theta)/sin(theta))q2
    // where theta = acos(q1 dot q2)
    if (dot < 1.0 - 1e-4) {
        const double angle = acos(dot);
        const double a = sin(angle * (1-t))/sin(angle);
        const double b = sin(angle * t)/sin(angle);

        Quaternion q1a = {0};
        quatcopy(q1a, q1);
        quatmult(q1a, a);
        quatmult(q, b);
        quatadd(q, q1a);
    } else {
        Quaternion tmp;

        quatcopy(q, q1);
        quatmult(q, (1-t));

        quatcopy(tmp, q2);
        quatmult(tmp, t);

        quatadd(q, b);
    }
}
