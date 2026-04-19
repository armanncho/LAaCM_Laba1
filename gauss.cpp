#include "gauss.h"
#include <cmath>
#include <stdexcept>

Vector gauss_without_pivot(Matrix A, Vector b) {
    int n = A.size();

    // прямой ход
    for (int k = 0; k < n; k++) {
        if (std::abs(A[k][k]) < 1e-18)
            throw std::runtime_error("Zero pivot encountered");

        for (int i = k + 1; i < n; i++) {
            double ratio =  A[i][k] / A[k][k];


            for (int j = k; j < n; j++) {
                A[i][j] -= ratio * A[k][j];
            }

            b[i] -= ratio * b[k];
        }
    }

    // обратный ход - находим корни
    Vector x(n);

    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;

        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }

        if (std::abs(A[i][i]) < 1e-18)
            throw std::runtime_error("Zero on diagonal during back substitution");

        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;

}

Vector gauss_with_pivot(Matrix A, Vector b) {
    int n = A.size();

    // прямой ход
    for (int k = 0; k < n; k++) {
        int index_max = k;

        for (int i = k + 1; i < n; i++) {
            if (std::abs(A[i][k]) > std::abs(A[index_max][k])) { // сравнение по модулю
                index_max = i;
            }
        }

        std::swap(A[k], A[index_max]);
        std::swap(b[k], b[index_max]);

        if (std::abs(A[k][k]) < 1e-18)
            throw std::runtime_error("Zero pivot encountered");

        for (int i = k + 1; i < n; i++) {
            double ratio = A[i][k] / A[k][k];

            for (int j = k; j < n; j++) {
                A[i][j] -= ratio * A[k][j];
            }

            b[i] -= ratio * b[k];
        }
    }

    Vector x(n);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;

        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }

        if (std::abs(A[i][i]) < 1e-18)
            throw std::runtime_error("Zero on diagonal");

        x[i] = (b[i] - sum)/ A[i][i];
    }

    return x;

}