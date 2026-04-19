#ifndef LINAL1_LU_H
#define LINAL1_LU_H

#include "matrix.h"

void lu_decomposition(const Matrix& A, Matrix& L, Matrix& U);

Vector forward_substitution(const Matrix& L, const Vector& b);
Vector backward_substitution(const Matrix& U, const Vector& y);

Vector solve_lu(const Matrix& L,const Matrix& U, const Vector& b);

#endif //LINAL1_LU_H
