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

void quat2rot(mfem::DenseMatrix& Q, mfem::Vector& q);
void rot2quat(mfem::Vector& q, mfem::DenseMatrix& Q);
void orient(mfem::DenseMatrix& Q_out, mfem::DenseMatrix& Q, double a, double b);
void axis(mfem::DenseMatrix& Q, mfem::Vector& psi, mfem::Vector &phi);
void bislerp(mfem::DenseMatrix& Qab, mfem::DenseMatrix& Qa, mfem::DenseMatrix& Qb, double t);


void laplace(
    mfem::GridFunction *x,
    mfem::Mesh& mesh,
    mfem::Array<int> &ess_bdr,
    mfem::Array<int> &nonzero_ess_bdr,
    mfem::Array<int> &zero_ess_bdr,
    int apex,
    int verbose);

void calculate_gradients(double* grads, mfem::GridFunction& x, mfem::Mesh& mesh, mfem::Table* v2e);

void define_fibers(
    mfem::Mesh& mesh,
    const double *phi_epi,
    const double *phi_lv,
    const double *phi_rv,
    const double *psi_ab,
    double *grad_phi_epi,
    double *grad_phi_lv,
    double *grad_phi_rv,
    double *grad_psi_ab,
    double alpha_endo,
    double alpha_epi,
    double beta_endo,
    double beta_epi,
    double *F,
    double *S,
    double *T);
#endif
