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

void orient(mfem::DenseMatrix& Q_out, mfem::DenseMatrix& Q, double a, double b);
void axis(mfem::DenseMatrix& Q, mfem::Vector& psi, mfem::Vector &phi);
void bislerp(mfem::DenseMatrix& Qab, mfem::DenseMatrix& Qa, mfem::DenseMatrix& Qb, double t);
#endif
