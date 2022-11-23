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

#ifndef CALCULUS_HPP
#define CALCULUS_HPP
#include "mfem.hpp"

const double PI=3.14159265358979323846;

typedef struct vec3 {
    double data[3];
    MFEM_HOST_DEVICE
    double& operator[](int i) {
        return data[i];
    }

    MFEM_HOST_DEVICE
    vec3 operator-(vec3& b) {
        vec3 out = {data[0] -b.data[0], data[1] - b.data[1], data[2] - b.data[2] };
        return out;
    }

    MFEM_HOST_DEVICE
    void operator+=(vec3& b) {
        data[0] += b.data[0];
        data[1] += b.data[1];
        data[2] += b.data[2];
    }

    MFEM_HOST_DEVICE
    void operator-=(vec3& b) {
        data[0] -= b.data[0];
        data[1] -= b.data[1];
        data[2] -= b.data[2];
    }

    MFEM_HOST_DEVICE
    void operator=(vec3&b) {
        data[0] = b.data[0];
        data[1] = b.data[1];
        data[2] = b.data[2];
    }

    MFEM_HOST_DEVICE
    void operator*=(double b) {
        data[0] *= b;
        data[1] *= b;
        data[2] *= b;
    }
} vec3;

typedef struct quat {
    double data[4];

    MFEM_HOST_DEVICE double& operator[](int i) {
        return data[i];
    }

    MFEM_HOST_DEVICE
    void operator=(quat& b) {
        data[0] = b.data[0];
        data[1] = b.data[1];
        data[2] = b.data[2];
        data[3] = b.data[3];
    }

    MFEM_HOST_DEVICE
    void operator*=(double b) {
        data[0] *= b;
        data[1] *= b;
        data[2] *= b;
        data[3] *= b;
    }

    MFEM_HOST_DEVICE
    void operator+=(quat& b) {
        data[0] += b.data[0];
        data[1] += b.data[1];
        data[2] += b.data[2];
        data[3] += b.data[3];
    }
} quat;

typedef struct mat3x3 {
    double data[3][3];
    MFEM_HOST_DEVICE double* operator[](int i) {
        return (double *) &data[i];
    }
} mat3x3;

MFEM_HOST_DEVICE double vecdot(vec3& a, vec3& b);
MFEM_HOST_DEVICE void cross(vec3& c, vec3& a, vec3& b);

MFEM_HOST_DEVICE void quat2rot(mat3x3& Q, quat& q);
MFEM_HOST_DEVICE void rot2quat(quat& q, mat3x3& Q);
MFEM_HOST_DEVICE void orient(mat3x3& Q_out, mat3x3& Q, double a, double b);
MFEM_HOST_DEVICE void axis(mat3x3& Q, vec3& psi, vec3& phi);
MFEM_HOST_DEVICE void bislerp(mat3x3& Qab, mat3x3& Qa, mat3x3& Qb, double t);

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
    double tol);
#endif
