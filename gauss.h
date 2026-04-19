#ifndef LINAL1_GAUSS_H
#define LINAL1_GAUSS_H

#include "matrix.h"

Vector gauss_without_pivot(Matrix A, Vector b);
Vector gauss_with_pivot(Matrix A, Vector b);


#endif //LINAL1_GAUSS_H
