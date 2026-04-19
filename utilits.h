#ifndef LINAL1_UTILITS_H
#define LINAL1_UTILITS_H


#include "matrix.h"

double norm(const Vector& v);
double discrepancy(const Matrix& A,const Vector& x, const Vector& b);
double relative_error(const Vector& x, const Vector& x_true);


#endif //LINAL1_UTILITS_H
