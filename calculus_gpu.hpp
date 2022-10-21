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

#ifndef CACLULUS_HPP
#define CACLULUS_HPP
#include "mfem.hpp"

const double PI=3.14159265;

typedef struct Vector3D {
    double data[2];
    MFEM_HOST_DEVICE double& operator[](int i) {
        return data[i];
    }
} Vector3D;

typedef struct Quaternion {
    double data[4];
    MFEM_HOST_DEVICE double& operator[](int i) {
        return data[i];
    }
} Quaternion;

typedef struct Matrix3x3 {
    double data[3][3];
    MFEM_HOST_DEVICE double* operator[](int i) {
        return (double *) &data[i];
    }
} Matrix3x3;

MFEM_HOST_DEVICE void veccopy(Vector3D& a, Vector3D& b);
MFEM_HOST_DEVICE void vecmul(Vector3D& a, double b);
MFEM_HOST_DEVICE double vecdot(Vector3D& a, Vector3D& b);

MFEM_HOST_DEVICE void quat2rot(Matrix3x3& Q, Quaternion& q);
MFEM_HOST_DEVICE void rot2quat(Quaternion& q, Matrix3x3& Q);
MFEM_HOST_DEVICE void orient(Matrix3x3& Q_out, Matrix3x3& Q, double a, double b);
MFEM_HOST_DEVICE void axis(Matrix3x3& Q, Vector3D& psi, Vector3D& phi);
MFEM_HOST_DEVICE void bislerp(Matrix3x3& Qab, Matrix3x3& Qa, Matrix3x3& Qb, double t);

void par_calculate_gradients(double* grads, mfem::ParGridFunction& x, mfem::ParMesh& mesh, mfem::Table* v2e);

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
    double *T);
#endif
