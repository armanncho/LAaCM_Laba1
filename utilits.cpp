#include "utilits.h"
#include <chrono> // для времени
#include <cmath>
#include <functional>

double norm(const Vector &v) {
    int n = v.size();
    double sum = 0.0;

    for (int i = 0; i < n; i++) {
        sum += v[i] * v[i];
    }

    return std::sqrt(sum);
}

double discrepancy(const Matrix& A,const Vector& x, const Vector& b) { // невязка = несоотвествие
    int n = A.size();

    Vector Ax(n, 0.0);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Ax[i] += A[i][j] * x[j];
        }
    }

    Vector r(n);
    for (int i = 0; i < n; i++) {
        r[i] = Ax[i] - b[i];
    }

    return norm(r);
}

double relative_error(const Vector& x, const Vector& x_true) {
    int n = x.size();

    Vector diff(n); // x - x_true

    for (int i = 0; i < n; i++) {
        diff[i] = x[i] - x_true[i];
    }

    return norm(diff) / norm(x_true);
}