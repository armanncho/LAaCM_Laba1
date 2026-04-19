#include <iostream>
#include <vector>
#include <chrono>

#include "matrix.h"
#include "gauss.h"
#include "lu.h"
#include "utilits.h"

using clock_type = std::chrono::high_resolution_clock;

double time_now() {
    return std::chrono::duration<double>(clock_type::now().time_since_epoch()).count();
}

int main() {

    // Эксперимент — СРАВНЕНИЕ ОДНОЙ СИСТЕМЫ

    std::cout << "===== EXPERIMENT 4.1 =====\n";


    for (int n : sizes) {
        Matrix A = generate_matrix(n);
        Vector b = generate_vector(n);

        // Gauss без pivot
        auto start = clock_type::now();
        gauss_without_pivot(A, b);
        auto end = clock_type::now();
        double t_gauss_no = std::chrono::duration<double>(end - start).count();

        //  Gauss с pivot
        start = clock_type::now();
        gauss_with_pivot(A, b);
        end = clock_type::now();
        double t_gauss_pivot = std::chrono::duration<double>(end - start).count();

        //  LU decomposition
        Matrix L, U;

        start = clock_type::now();
        lu_decomposition(A, L, U);
        end = clock_type::now();
        double t_lu_decomp = std::chrono::duration<double>(end - start).count();

        //  LU solve
        start = clock_type::now();
        solve_lu(L, U, b);
        end = clock_type::now();
        double t_lu_solve = std::chrono::duration<double>(end - start).count();

        double t_lu_total = t_lu_decomp + t_lu_solve;

        std::cout << "n = " << n << "\n";
        std::cout << "Gauss (no pivot): " << t_gauss_no << "\n";
        std::cout << "Gauss (pivot): " << t_gauss_pivot << "\n";
        std::cout << "LU decomp: " << t_lu_decomp << "\n";
        std::cout << "LU solve: " << t_lu_solve << "\n";
        std::cout << "LU total: " << t_lu_total << "\n";
        std::cout << "--------------------------\n";
    }


    // =========================================
    // Эксперимент — МНОЖЕСТВЕННЫЕ ПРАВЫЕ ЧАСТИ
    // =========================================
    std::cout << "\n===== EXPERIMENT 4.2 =====\n";

    int n = 500;
    Matrix A = generate_matrix(n);


    for (int k : ks) {
        std::vector<Vector> bs;

        for (int i = 0; i < k; i++) {
            bs.push_back(generate_vector(n));
        }

        //  Gauss (pivot)
        auto start = clock_type::now();

        for (auto& b : bs) {
            gauss_with_pivot(A, b);
        }

        auto end = clock_type::now();
        double t_gauss = std::chrono::duration<double>(end - start).count();

        // LU
        Matrix L, U;

        auto start_lu = clock_type::now();

        lu_decomposition(A, L, U);

        for (auto& b : bs) {
            solve_lu(L, U, b);
        }

        auto end_lu = clock_type::now();
        double t_lu = std::chrono::duration<double>(end_lu - start_lu).count();

        std::cout << "k = " << k << "\n";
        std::cout << "Gauss (pivot): " << t_gauss << "\n";
        std::cout << "LU: " << t_lu << "\n";
        std::cout << "--------------------------\n";
    }


    // Эксперимент  — МАТРИЦА ГИЛЬБЕРТА
    std::cout << "\n===== EXPERIMENT 4.3 =====\n";

    std::vector<int> sizes_h = {5, 10, 15};

    for (int n : sizes_h) {
        Matrix H = generate_hilbert_matrix(n);

        Vector x_true(n, 1.0);
        Vector b = multiply(H, x_true);

        //  решения
        Vector x_gauss, x_pivot, x_lu;

        // Gauss без pivot
        bool gauss_ok = true;
        try {
            x_gauss = gauss_without_pivot(H, b);
        } catch (...) {
            gauss_ok = false;
        }

        // Pivot
        bool pivot_ok = true;
        try {
            x_pivot = gauss_with_pivot(H, b);
        } catch (...) {
            pivot_ok = false;
        }

        // LU
        bool lu_ok = true;
        Matrix L, U;
        try {
            lu_decomposition(H, L, U);
            x_lu = solve_lu(L, U, b);
        } catch (...) {
            lu_ok = false;
        }


        //  вывод
        std::cout << "n = " << n << "\n";

        if (gauss_ok) {
            std::cout << "Gauss error: "
                      << relative_error(x_gauss, x_true) << "\n";
            std::cout << "Gauss residual: "
                      << discrepancy(H, x_gauss, b) << "\n";
        } else {
            std::cout << "Gauss failed\n";
        }

        if (pivot_ok) {
            std::cout << "Pivot error: "
                      << relative_error(x_pivot, x_true) << "\n";
            std::cout << "Pivot residual: "
                      << discrepancy(H, x_pivot, b) << "\n";
        } else {
            std::cout << "Pivot failed\n";
        }

        if (lu_ok) {
            std::cout << "LU error: "
                      << relative_error(x_lu, x_true) << "\n";
            std::cout << "LU residual: "
                      << discrepancy(H, x_lu, b) << "\n";
        } else {
            std::cout << "LU failed\n";
        }

        std::cout << "--------------------------\n";
    }

    return 0;
}