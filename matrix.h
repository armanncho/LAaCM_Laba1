#ifndef LINAL1_MATRIX_H
#define LINAL1_MATRIX_H

#include <vector>

using Matrix = std::vector<std::vector<double>>;
using Vector = std::vector<double>;

Matrix generate_matrix(int n);
Vector generate_vector(int n);
Matrix generate_hilbert_matrix(int n);

Vector multiply(const Matrix& A, const Vector& x);


#endif //LINAL1_MATRIX_H
