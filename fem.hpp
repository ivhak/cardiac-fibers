// Copyright (C) 2022 Iver Håkonsen
//
// cardiac-fibers is free software: you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option) any
// later version.
//
// cardiac-fibers is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
// details.
//
// You should have received a copy of the GNU Lesser General Public License along with
// cardiac-fibers.  If not, see <https://www.gnu.org/licenses/>.
//
// Authors: Iver Håkonsen <hakonseniver@yahoo.no

#ifndef FEM_HPP
#define FEM_HPP

#include "mfem.hpp"

void project_h1_to_l2(
        double *l2_vals,
        const double *h1_vals,
        const int ne,
        const int *h1_table_col,
        const int *h1_table_row,
        const int *l2_table_col,
        const int *l2_table_row);

void compute_gradient(
        double *gradient,
        const double *laplace,
        const double *vert,
        const int ne,
        const int nv,
        const int *h1_table_col,
        const int *h1_table_row,
        const int *l2_table_col,
        const int *l2_table_row);

void interpolate_gradient_to_h1(
        double *h1_vals,
        const double *l2_vals,
        const int nv,
        const int *v2e_table_col,
        const int *v2e_table_row,
        const int *l2_table_col,
        const int *l2_table_row);
#endif
