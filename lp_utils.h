#ifndef LP_UTILS_H
#define LP_UTILS_H

#include "lp_types.h"
#include <string>

LinearProgram read_problem();
CanonicalProblem toCanonicalForm(const LinearProgram& lp);
LinearProgram buildDual(const LinearProgram& primal);
SymmetricProblem toSymmetricForm(const LinearProgram& lp);
std::vector<double> recoverSolution(const CanonicalProblem& cp, const std::vector<double>& x_canon);
void print_linear_program(const LinearProgram& lp, const std::string& title = "");
void print_canonical_problem(const CanonicalProblem& cp, const std::string& title = "");
void print_symmetric_problem(const SymmetricProblem& sp, const std::string& title = "");

#endif
