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
