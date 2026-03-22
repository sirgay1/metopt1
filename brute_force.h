// brute_force.h
#ifndef BRUTE_FORCE_H
#define BRUTE_FORCE_H

#include "lp_types.h"
#include <vector>

struct BruteForceResult {
    bool feasible;          // допустима ли задача
    bool bounded;           // ограничена ли целевая функция
    double optimalValue;    // оптимальное значение
    std::vector<double> solution; // вектор решения (размер n)
};

BruteForceResult bruteForceSolve(const CanonicalProblem& cp);

#endif
