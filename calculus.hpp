#ifndef CACLULUS_HPP
#define CACLULUS_HPP
#include "mfem.hpp"

const double PI=3.14159265;

void orient(mfem::DenseMatrix& Q_out, mfem::DenseMatrix& Q, double a, double b);
void axis(mfem::DenseMatrix& Q, mfem::Vector& psi, mfem::Vector &phi);

#endif
