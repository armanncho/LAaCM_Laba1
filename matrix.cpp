#include <vector>
#include <random>

// Псевдонимы типов для удобства
using Vector = std::vector<double>;
using Matrix = std::vector<std::vector<double>>;

static std::mt19937 generator(52); // фиксированный seed для воспроизводимости

// Генерирует случайную матрицу n x n
Matrix generate_matrix(int n) {
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

    Matrix A(n, std::vector<double>(n));

    for(int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            A[i][j] = distribution(generator);

    return A;
}

// Генерирует случайную правую часть (вектор)
Vector generate_vector(int n) {
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

    Vector v(n);

    for (int i = 0; i < n; i++)
        v[i] = distribution(generator);

    return v;
}

// Генерирует матрицу Гильберта для проверки точности
Matrix generate_hilbert_matrix(int n) {
    Matrix H(n, std::vector<double>(n));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Формула H_ij = 1 / (i + j + 1) для индексов с 0
            H[i][j] = 1.0 / (i + j + 1);
        }
    }

    return H;
}

// Умножение матрицы на вектор (нужно для вычисления невязки и правой части b=Hx)
Vector multiply(const Matrix& A, const Vector& x) {
    int n = A.size();
    Vector b(n, 0.0);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            b[i] += A[i][j] * x[j];

    return b;
}