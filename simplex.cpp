#include "simplex.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

SimplexSolver::SimplexSolver(const vector<vector<double>>& A,
                             const vector<double>& b,
                             const vector<double>& c,
                             bool min,
                             bool verb)
    : m(A.size()), numVars(A[0].size() + m), minimize(min), verbose(verb) {
    if (!minimize) {
        original_c.resize(A[0].size());
        for (size_t j = 0; j < A[0].size(); ++j) original_c[j] = -c[j];
    } else {
        original_c = c;
    }

    int n_orig = A[0].size();
    tableau.assign(m+1, vector<double>(numVars + 1, 0.0));
    basis.resize(m);

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n_orig; ++j) {
            tableau[i][j] = A[i][j];
        }
        tableau[i][n_orig + i] = 1.0;
        tableau[i][numVars] = b[i];
    }

    for (int j = 0; j < n_orig; ++j) {
        double sumA = 0.0;
        for (int i = 0; i < m; ++i) sumA += A[i][j];
        tableau[m][j] = sumA;
    }
    double sumB = 0.0;
    for (int i = 0; i < m; ++i) sumB += b[i];
    tableau[m][numVars] = sumB;

    for (int i = 0; i < m; ++i) {
        basis[i] = n_orig + i;
    }
}

void SimplexSolver::printTableau(const string& phase) {
    cout << "\n=== " << phase << " ===\n";
    cout << "Итерация " << iteration << "\n";
    cout << "Базис: ";
    for (int i = 0; i < m; ++i) {
        if (basis[i] < numVars)
            cout << "x" << basis[i] + 1 << " ";
        else
            cout << "y" << basis[i] - (numVars - m) + 1 << " ";
    }
    cout << "\nТаблица:\n";
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= numVars; ++j) {
            cout << setw(10) << tableau[i][j] << " ";
        }
        cout << "\n";
    }
    vector<double> x(numVars, 0.0);
    for (int i = 0; i < m; ++i) {
        if (basis[i] < numVars)
            x[basis[i]] = tableau[i][numVars];
    }
    cout << "Текущее решение: ";
    for (int j = 0; j < numVars; ++j) {
        cout << "x" << j+1 << "=" << x[j] << " ";
    }
    cout << "\nЗначение целевой функции: " << tableau[m][numVars] << "\n";
}

int SimplexSolver::chooseEnteringColumnPhaseI() {
    for (int j = 0; j < numVars; ++j) {
        if (tableau[m][j] > 1e-9) {
            bool hasPositive = false;
            for (int i = 0; i < m; ++i) {
                if (tableau[i][j] > 1e-9) {
                    hasPositive = true;
                    break;
                }
            }
            if (hasPositive) {
                return j;
            }
        }
    }
    return -1;
}

int SimplexSolver::chooseEnteringColumnPhaseII() {
    for (int j = 0; j < numVars; ++j) {
        if (tableau[m][j] > 1e-9) {
            return j;
        }
    }
    return -1;
}

int SimplexSolver::chooseLeavingRow(int col) {
    int leaving = -1;
    double minTheta = 1e100;
    for (int i = 0; i < m; ++i) {
        if (tableau[i][col] > 1e-9) {
            double theta = tableau[i][numVars] / tableau[i][col];
            if (theta < minTheta - 1e-9) {
                minTheta = theta;
                leaving = i;
            } else if (fabs(theta - minTheta) < 1e-9 && leaving != -1 && i < leaving) {
                leaving = i;
            }
        }
    }
    return leaving;
}

void SimplexSolver::pivot(int row, int col) {
    double pivotVal = tableau[row][col];
    for (int j = 0; j <= numVars; ++j) {
        tableau[row][j] /= pivotVal;
    }
    for (int i = 0; i <= m; ++i) {
        if (i != row && fabs(tableau[i][col]) > 1e-12) {
            double factor = tableau[i][col];
            for (int j = 0; j <= numVars; ++j) {
                tableau[i][j] -= factor * tableau[row][j];
            }
        }
    }
    basis[row] = col;
}

bool SimplexSolver::phaseI() {
    if (verbose) printTableau("ФАЗА I (начало)");
    iteration = 0;
    while (iteration < 1000) {
        int entering = chooseEnteringColumnPhaseI();
        if (entering == -1) break;

        int leaving = chooseLeavingRow(entering);
        if (leaving == -1) {
            if (verbose) cout << "В фазе I нет положительных элементов в разрешающем столбце. Задача недопустима.\n";
            return false;
        }

        if (verbose) {
            cout << "Итерация " << iteration+1 << ": входит x" << entering+1
                 << ", выходит ";
            if (basis[leaving] < numVars)
                cout << "x" << basis[leaving]+1;
            else
                cout << "y" << basis[leaving] - (numVars - m) + 1;
            cout << "\n";
        }
        pivot(leaving, entering);
        iteration++;
        if (verbose) printTableau("ФАЗА I");
    }
    double W = tableau[m][numVars];
    if (W > 1e-6) {
        if (verbose) cout << "Минимальное значение W = " << W << " > 0, задача недопустима.\n";
        return false;
    }
    return true;
}

void SimplexSolver::setupPhaseII(int orig_n) {
    int newNumVars = orig_n;
    vector<vector<double>> newTableau(m+1, vector<double>(newNumVars+1, 0.0));
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j < newNumVars; ++j) {
            newTableau[i][j] = tableau[i][j];
        }
        newTableau[i][newNumVars] = tableau[i][numVars];
    }
    tableau = move(newTableau);
    numVars = newNumVars;

    for (int i = 0; i < m; ++i) {
        if (basis[i] >= numVars) {
            for (int j = 0; j < numVars; ++j) {
                if (fabs(tableau[i][j]) > 1e-9) {
                    pivot(i, j);
                    break;
                }
            }
        }
    }

    vector<double> cB(m, 0.0);
    for (int i = 0; i < m; ++i) {
        if (basis[i] < numVars)
            cB[i] = original_c[basis[i]];
    }
    double Z0 = 0.0;
    for (int i = 0; i < m; ++i) {
        Z0 += cB[i] * tableau[i][numVars];
    }
    vector<double> reducedCost(numVars, 0.0);
    for (int j = 0; j < numVars; ++j) {
        reducedCost[j] = original_c[j];
        for (int i = 0; i < m; ++i) {
            reducedCost[j] -= cB[i] * tableau[i][j];
        }
    }
    for (int j = 0; j < numVars; ++j) {
        tableau[m][j] = -reducedCost[j];
    }
    tableau[m][numVars] = Z0;
}

bool SimplexSolver::phaseII() {
    if (verbose) printTableau("ФАЗА II (начало)");
    iteration = 0;
    while (iteration < 1000) {
        int entering = chooseEnteringColumnPhaseII();
        if (entering == -1) break;

        int leaving = chooseLeavingRow(entering);
        if (leaving == -1) {
            if (verbose) cout << "В фазе II нет положительных элементов в разрешающем столбце. Задача неограничена.\n";
            return false;
        }

        if (verbose) {
            cout << "Итерация " << iteration+1 << ": входит x" << entering+1
                 << ", выходит x" << basis[leaving]+1 << "\n";
        }
        pivot(leaving, entering);
        iteration++;
        if (verbose) printTableau("ФАЗА II");
    }
    return true;
}

SimplexSolver::Result SimplexSolver::solve() {
    Result res;
    if (!phaseI()) {
        res.feasible = false;
        res.bounded = false;
        return res;
    }

    setupPhaseII(original_c.size());

    if (!phaseII()) {
        res.feasible = true;
        res.bounded = false;
        return res;
    }

    res.feasible = true;
    res.bounded = true;
    vector<double> x(numVars, 0.0);
    for (int i = 0; i < m; ++i) {
        if (basis[i] < numVars)
            x[basis[i]] = tableau[i][numVars];
    }
    res.solution = x;
    res.optimalValue = tableau[m][numVars];
    if (!minimize) res.optimalValue = -res.optimalValue;
    return res;
}
