#include "lp_utils.h"
#include "simplex.h"
#include "brute_force.h"
#include <iostream>
#include <vector>

using namespace std;

int main() {
    LinearProgram lp = read_problem();

    cout << "\n=== Исходная задача ===\n";
    print_linear_program(lp, "Прямая задача");

    LinearProgram dual_lp = buildDual(lp);
    cout << "\n=== Двойственная задача ===\n";
    print_linear_program(dual_lp, "Двойственная задача");

    CanonicalProblem cp = toCanonicalForm(lp);
    cout << "\n=== Каноническая форма ===\n";
    print_canonical_problem(cp, "Каноническая форма");

    // ----- Решение прямой задачи -----
    cout << "\n=== Решение прямой задачи (симплекс) ===\n";
    SimplexSolver primalSolver(cp.A, cp.b, cp.c, true, true);
    auto primalRes = primalSolver.solve();

    if (!primalRes.feasible) {
        cout << "\nПрямая задача не имеет допустимых решений.\n";
    } else if (!primalRes.bounded) {
        cout << "\nПрямая задача неограничена.\n";
    } else {
        cout << "\nОптимум прямой задачи: " << primalRes.optimalValue << "\n";
        vector<double> x_orig = recoverSolution(cp, primalRes.solution);
        cout << "Решение прямой задачи (в исходных переменных):\n";
        for (size_t i = 0; i < x_orig.size(); ++i)
            cout << "x" << i+1 << " = " << x_orig[i] << "\n";
    }

    // ---- Далее остаётся существующий код для двойственной задачи и прочих форм ----
    LinearProgram cp_as_lp;
    cp_as_lp.n = cp.n;
    cp_as_lp.m = cp.m;
    cp_as_lp.A = cp.A;
    cp_as_lp.b = cp.b;
    cp_as_lp.c = cp.c;
    cp_as_lp.varSigns.assign(cp.n, NONNEGATIVE);
    cp_as_lp.constrSigns.assign(cp.m, EQ);

    LinearProgram dual_cp = buildDual(cp_as_lp);

    cout << "\n=== Двойственная к канонической ===\n";
    print_linear_program(dual_cp, "Двойственная к канонической");

    SymmetricProblem sp = toSymmetricForm(lp);

    cout << "\n=== Симметричная форма ===\n";
    print_symmetric_problem(sp, "Симметричная форма");

    LinearProgram sp_as_lp;
    sp_as_lp.n = sp.n;
    sp_as_lp.m = sp.m;
    sp_as_lp.A = sp.A;
    sp_as_lp.b = sp.b;
    sp_as_lp.c = sp.c;
    sp_as_lp.varSigns.assign(sp.n, NONNEGATIVE);
    sp_as_lp.constrSigns.assign(sp.m, LE);

    LinearProgram dual_sp = buildDual(sp_as_lp);

    cout << "\n=== Двойственная к симметричной ===\n";
    print_linear_program(dual_sp, "Двойственная к симметричной");

    cout << "\n=== Решение двойственной задачи (симплекс) ===\n";

    CanonicalProblem dual_cp2 = toCanonicalForm(dual_lp);

    SimplexSolver solver(dual_cp2.A, dual_cp2.b, dual_cp2.c, true, true);
    auto res = solver.solve();

    if (!res.feasible) {
        cout << "\nДвойственная задача не имеет допустимых решений.\n";
    } else if (!res.bounded) {
        cout << "\nДвойственная задача неограничена.\n";
    } else {
        double dualOpt = res.optimalValue;
        if (dual_lp.sense == MAXIMIZE) dualOpt = -dualOpt;   // корректировка знака
        cout << "\nОптимум двойственной: " << dualOpt << "\n";

        vector<double> x_orig = recoverSolution(dual_cp2, res.solution);

        cout << "Решение двойственной:\n";
        for (size_t i = 0; i < x_orig.size(); ++i)
            cout << "y" << i+1 << " = " << x_orig[i] << "\n";
    }
    
    // Решение перебором вершин
    BruteForceResult bfResult = bruteForceSolve(cp);
    if (bfResult.feasible) {
        std::cout << "Оптимум (перебор): " << bfResult.optimalValue << "\n";
        std::vector<double> x_orig = recoverSolution(cp, bfResult.solution);
        std::cout << "Решение: ";
        for (double val : x_orig) std::cout << val << " ";
        std::cout << "\n";
    } else {
        std::cout << "Задача недопустима (по результатам перебора).\n";
    }
    
    return 0;
}
