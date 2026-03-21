#ifndef LP_TYPES_H
#define LP_TYPES_H

#include <vector>

enum ConstraintSign { LE, GE, EQ };
enum VariableSign { NONNEGATIVE, NONPOSITIVE, FREE };
enum OptimizationSense { MINIMIZE, MAXIMIZE };

struct LinearProgram {
    int n;
    int m;
    std::vector<double> c;
    std::vector<VariableSign> varSigns;
    std::vector<std::vector<double>> A;
    std::vector<ConstraintSign> constrSigns;
    std::vector<double> b;
    OptimizationSense sense;   // добавлено
};

struct CanonicalProblem {
    int n;
    int m;
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<std::vector<std::pair<int, double>>> mapping;
    OptimizationSense sense;   // добавлено
};

struct SymmetricProblem {
    int n;
    int m;
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> c;
    OptimizationSense sense;   // добавлено
};
#endif
