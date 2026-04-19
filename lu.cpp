#include "lu.h"
#include <cmath>
#include <stdexcept>

// LU РАЗЛОЖЕНИЕ

void lu_decomposition(const Matrix &A, Matrix &L, Matrix &U) {
    int n = A.size();

    L = Matrix(n, Vector(n, 0.0));
    U = Matrix(n, Vector(n, 0.0));

    for (int i = 0; i < n; i++) {
        L[i][i] = 1.0;
    }

    for (int k = 0; k < n; k++) {

        // --- считаем U ---
        for (int j = k; j < n; j++) {
            double sum = 0.0;

            for (int p = 0; p < k; p++) {
                sum += L[k][p] * U[p][j];
            }

            U[k][j] = A[k][j] - sum;
        }
        if (std::abs(U[k][k]) < 1e-18)
            throw std::runtime_error("Zero on diagonal");

        // --- считаем L ---
        for (int i = k + 1; i < n; i++) {
            double sum = 0.0;

            for (int p = 0; p < k; p++) {
                sum += L[i][p] * U[p][k];
            }

            L[i][k] = (A[i][k] - sum) / U[k][k];
        }
    }
}


// ПРЯМОЙ ХОД (Ly = b)

Vector forward_substitution(const Matrix& L, const Vector& b) {
    int n = L.size();
    Vector y(n);

    for (int i = 0; i < n; i++) {
        double sum = 0.0;

        for (int j = 0; j < i; j++)
            sum += L[i][j] * y[j];

        if (std::abs(L[i][i]) < 1e-18)
            throw std::runtime_error("Zero on diagonal");

        y[i] = (b[i] - sum) / L[i][i];
    }

    return y;
}


// ОБРАТНЫЙ ХОД (Ux = y)

Vector backward_substitution(const Matrix& U, const Vector& y) {
    int n = U.size();
    Vector x(n);

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;

        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }

        if (std::abs(U[i][i]) < 1e-18)
            throw::std::runtime_error("Zero on diagonal");

        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}

// РЕШЕНИЕ Ax = b через LU

Vector solve_lu(const Matrix& L, const Matrix& U, const Vector& b) {
    Vector y = forward_substitution(L, b);
    return backward_substitution(U, y);
}