// brute_force.cpp
#include "brute_force.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <iostream>

const double EPS = 1e-9;

// Решение СЛАУ Ax = b методом Гаусса с выбором главного элемента по столбцу.
// Возвращает true, если решение единственное, иначе false.
bool solveLinearSystem(std::vector<std::vector<double>> A,
                       std::vector<double> b,
                       std::vector<double>& x) {
    int n = A.size();
    // Расширенная матрица
    for (int i = 0; i < n; ++i) {
        A[i].push_back(b[i]);
    }

    // Прямой ход
    for (int col = 0; col < n; ++col) {
        // Поиск главного элемента в столбце col
        int maxRow = col;
        for (int row = col + 1; row < n; ++row) {
            if (std::fabs(A[row][col]) > std::fabs(A[maxRow][col])) {
                maxRow = row;
            }
        }
        if (std::fabs(A[maxRow][col]) < EPS) {
            return false; // Матрица вырождена
        }
        if (maxRow != col) {
            std::swap(A[col], A[maxRow]);
        }

        // Нормализация строки col
        double pivot = A[col][col];
        for (int j = col; j <= n; ++j) {
            A[col][j] /= pivot;
        }

        // Исключение в остальных строках
        for (int row = 0; row < n; ++row) {
            if (row != col && std::fabs(A[row][col]) > EPS) {
                double factor = A[row][col];
                for (int j = col; j <= n; ++j) {
                    A[row][j] -= factor * A[col][j];
                }
            }
        }
    }

    x.resize(n);
    for (int i = 0; i < n; ++i) {
        x[i] = A[i][n];
    }
    return true;
}

// Рекурсивная генерация сочетаний
void generateCombinations(int n, int m, int start, std::vector<int>& combination,
                          const std::vector<std::vector<double>>& A,
                          const std::vector<double>& b,
                          const std::vector<double>& c,
                          bool minimize,
                          BruteForceResult& best) {
    if (combination.size() == m) {
        // Собираем базисную матрицу B
        std::vector<std::vector<double>> B(m, std::vector<double>(m));
        for (int i = 0; i < m; ++i) {
            int colIdx = combination[i];
            for (int row = 0; row < m; ++row) {
                B[row][i] = A[row][colIdx];
            }
        }

        std::vector<double> xB;
        if (!solveLinearSystem(B, b, xB)) {
            return; // вырожденная матрица
        }

        // Проверка неотрицательности
        bool feasible = true;
        for (double val : xB) {
            if (val < -EPS) {
                feasible = false;
                break;
            }
        }
        if (!feasible) return;

        // Формируем полное решение
        int nVars = A[0].size();
        std::vector<double> x(nVars, 0.0);
        for (int i = 0; i < m; ++i) {
            x[combination[i]] = xB[i];
        }

        // Вычисляем значение целевой функции
        double value = 0.0;
        for (int j = 0; j < nVars; ++j) {
            value += c[j] * x[j];
        }

        // Обновляем лучшее решение
        if (!best.feasible) {
            best.feasible = true;
            best.bounded = true;
            best.optimalValue = value;
            best.solution = x;
        } else {
            if ((minimize && value < best.optimalValue - EPS) ||
                (!minimize && value > best.optimalValue + EPS)) {
                best.optimalValue = value;
                best.solution = x;
            }
        }
        return;
    }

    for (int i = start; i < n; ++i) {
        combination.push_back(i);
        generateCombinations(n, m, i + 1, combination, A, b, c, minimize, best);
        combination.pop_back();
    }
}

BruteForceResult bruteForceSolve(const CanonicalProblem& cp) {
    BruteForceResult result;
    result.feasible = false;
    result.bounded = true; // по умолчанию считаем, что область ограничена (если нет признаков неограниченности)

    int m = cp.m;   // число уравнений
    int n = cp.n;   // число переменных

    if (n < m) {
        // Не может быть допустимых решений (система переопределена)
        return result;
    }

    // Минимизация или максимизация? (в канонической форме всегда MINIMIZE)
    bool minimize = (cp.sense == MINIMIZE);

    std::vector<int> combination;
    generateCombinations(n, m, 0, combination, cp.A, cp.b, cp.c, minimize, result);

    // Если не нашли ни одного БДР – задача недопустима
    if (!result.feasible) {
        result.bounded = false;
    }

    return result;
}
