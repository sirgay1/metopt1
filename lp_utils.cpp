#include "lp_utils.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>

using namespace std;

LinearProgram read_problem() {
    LinearProgram lp;
    cout << "Введите количество переменных (n): "; cin >> lp.n;
    cout << "Введите количество ограничений (m): "; cin >> lp.m;
    int sense;
    cout << "Направление оптимизации (0 - минимизация, 1 - максимизация): "; cin >> sense;
    lp.sense = static_cast<OptimizationSense>(sense);

    lp.varSigns.resize(lp.n);
    cout << "Введите знаки для каждой переменной (0 - >=0, 1 - <=0, 2 - свободная):\n";
    for (int i = 0; i < lp.n; ++i) {
        int s; cin >> s;
        lp.varSigns[i] = static_cast<VariableSign>(s);
    }

    lp.c.resize(lp.n);
    cout << "Введите коэффициенты целевой функции (c1 c2 ... cn):\n";
    for (int i = 0; i < lp.n; ++i) cin >> lp.c[i];

    lp.A.resize(lp.m, vector<double>(lp.n));
    lp.constrSigns.resize(lp.m);
    lp.b.resize(lp.m);
    cout << "Введите ограничения построчно.\n";
    for (int i = 0; i < lp.m; ++i) {
        cout << "Ограничение " << i + 1 << ":\n";
        cout << "  Введите " << lp.n << " коэффициентов: ";
        for (int j = 0; j < lp.n; ++j) cin >> lp.A[i][j];
        int sign; cout << "  Знак (0 - <=, 1 - >=, 2 - =): "; cin >> sign;
        lp.constrSigns[i] = static_cast<ConstraintSign>(sign);
        cout << "  Правая часть: "; cin >> lp.b[i];
    }
    return lp;
}

CanonicalProblem toCanonicalForm(const LinearProgram& lp) {
    CanonicalProblem cp;
    cp.m = lp.m;
    cp.mapping.resize(lp.n);
    int newVarCount = 0;
    vector<int> origToNewStart(lp.n);
    for (int j = 0; j < lp.n; ++j) {
        origToNewStart[j] = newVarCount;
        switch (lp.varSigns[j]) {
            case NONNEGATIVE:
                cp.mapping[j].push_back({newVarCount, 1.0}); newVarCount++; break;
            case NONPOSITIVE:
                cp.mapping[j].push_back({newVarCount, -1.0}); newVarCount++; break;
            case FREE:
                cp.mapping[j].push_back({newVarCount, 1.0});
                cp.mapping[j].push_back({newVarCount+1, -1.0});
                newVarCount += 2; break;
        }
    }
    cp.n = newVarCount;
    cp.A.assign(lp.m, vector<double>(cp.n, 0.0));
    cp.c.assign(cp.n, 0.0);
    cp.b = lp.b;

    for (int i = 0; i < lp.m; ++i)
        for (int j = 0; j < lp.n; ++j)
            if (lp.A[i][j] != 0)
                for (auto& p : cp.mapping[j])
                    cp.A[i][p.first] += lp.A[i][j] * p.second;

    // Преобразование целевой функции с учётом направления
    vector<double> coeff = lp.c;
    if (lp.sense == MAXIMIZE) {
        for (int j = 0; j < lp.n; ++j) coeff[j] = -coeff[j];
    }
    for (int j = 0; j < lp.n; ++j)
        if (coeff[j] != 0)
            for (auto& p : cp.mapping[j])
                cp.c[p.first] += coeff[j] * p.second;

    int slackCount = 0;
    for (int i = 0; i < lp.m; ++i)
        if (lp.constrSigns[i] == LE || lp.constrSigns[i] == GE) slackCount++;
    int old_n = cp.n;
    cp.n += slackCount;
    for (int i = 0; i < cp.m; ++i) cp.A[i].resize(cp.n, 0.0);
    cp.c.resize(cp.n, 0.0);
    int col = old_n;
    for (int i = 0; i < lp.m; ++i) {
        if (lp.constrSigns[i] == LE) cp.A[i][col++] = 1.0;
        else if (lp.constrSigns[i] == GE) cp.A[i][col++] = -1.0;
    }
    for (int i = 0; i < cp.m; ++i) {
        if (cp.b[i] < 0) {
            for (int j = 0; j < cp.n; ++j) cp.A[i][j] = -cp.A[i][j];
            cp.b[i] = -cp.b[i];
        }
    }
    cp.sense = MINIMIZE;
    return cp;
}

LinearProgram buildDual(const LinearProgram& primal) {
    LinearProgram dual;
    dual.n = primal.m;
    dual.m = primal.n;
    dual.c = primal.b;
    dual.varSigns.resize(dual.n);

    // Знаки переменных двойственной (определяются знаками ограничений прямой)
    if (primal.sense == MINIMIZE) {
        for (int i = 0; i < primal.m; ++i) {
            if (primal.constrSigns[i] == LE)
                dual.varSigns[i] = NONPOSITIVE;   // y_i <= 0
            else if (primal.constrSigns[i] == GE)
                dual.varSigns[i] = NONNEGATIVE;   // y_i >= 0
            else
                dual.varSigns[i] = FREE;          // y_i свободна
        }
    } else { // MAXIMIZE
        for (int i = 0; i < primal.m; ++i) {
            if (primal.constrSigns[i] == LE)
                dual.varSigns[i] = NONNEGATIVE;   // y_i >= 0
            else if (primal.constrSigns[i] == GE)
                dual.varSigns[i] = NONPOSITIVE;   // y_i <= 0
            else
                dual.varSigns[i] = FREE;
        }
    }

    // Транспонированная матрица
    dual.A.assign(dual.m, vector<double>(dual.n));
    for (int i = 0; i < primal.m; ++i)
        for (int j = 0; j < primal.n; ++j)
            dual.A[j][i] = primal.A[i][j];

    dual.b = primal.c;
    dual.constrSigns.resize(dual.m);

    // Знаки ограничений двойственной (определяются знаками переменных прямой)
    if (primal.sense == MINIMIZE) {
        for (int j = 0; j < primal.n; ++j) {
            if (primal.varSigns[j] == NONNEGATIVE)
                dual.constrSigns[j] = LE;      // (A^T y)_j <= c_j
            else if (primal.varSigns[j] == NONPOSITIVE)
                dual.constrSigns[j] = GE;      // (A^T y)_j >= c_j
            else
                dual.constrSigns[j] = EQ;      // (A^T y)_j = c_j
        }
    } else { // MAXIMIZE
        for (int j = 0; j < primal.n; ++j) {
            if (primal.varSigns[j] == NONNEGATIVE)
                dual.constrSigns[j] = GE;      // (A^T y)_j >= c_j
            else if (primal.varSigns[j] == NONPOSITIVE)
                dual.constrSigns[j] = LE;      // (A^T y)_j <= c_j
            else
                dual.constrSigns[j] = EQ;
        }
    }

    dual.sense = (primal.sense == MINIMIZE) ? MAXIMIZE : MINIMIZE;
    return dual;
}

SymmetricProblem toSymmetricForm(const LinearProgram& lp) {
    SymmetricProblem sp;

    // Замена переменных: свободные → разность двух неотрицательных
    vector<vector<int>> newIndices(lp.n);
    int new_n = 0;
    for (int j = 0; j < lp.n; ++j) {
        switch (lp.varSigns[j]) {
            case NONNEGATIVE:
                newIndices[j].push_back(new_n++);
                break;
            case NONPOSITIVE:
                newIndices[j].push_back(new_n++);
                break;
            case FREE:
                newIndices[j].push_back(new_n++);
                newIndices[j].push_back(new_n++);
                break;
        }
    }
    sp.n = new_n;
    sp.m = lp.m;
    sp.A.assign(sp.m, vector<double>(sp.n, 0.0));
    sp.b = lp.b;
    sp.c.assign(sp.n, 0.0);

    // Заполнение A и c с учётом замен
    for (int i = 0; i < lp.m; ++i) {
        for (int j = 0; j < lp.n; ++j) {
            double coeff = lp.A[i][j];
            if (coeff == 0) continue;
            if (lp.varSigns[j] == NONNEGATIVE) {
                sp.A[i][newIndices[j][0]] += coeff;
            } else if (lp.varSigns[j] == NONPOSITIVE) {
                sp.A[i][newIndices[j][0]] -= coeff;
            } else { // FREE
                sp.A[i][newIndices[j][0]] += coeff;
                sp.A[i][newIndices[j][1]] -= coeff;
            }
        }
    }
    for (int j = 0; j < lp.n; ++j) {
        double coeff = lp.c[j];
        if (coeff == 0) continue;
        if (lp.varSigns[j] == NONNEGATIVE) {
            sp.c[newIndices[j][0]] += coeff;
        } else if (lp.varSigns[j] == NONPOSITIVE) {
            sp.c[newIndices[j][0]] -= coeff;
        } else { // FREE
            sp.c[newIndices[j][0]] += coeff;
            sp.c[newIndices[j][1]] -= coeff;
        }
    }

    // Преобразование ограничений: все должны быть типа ≤
    vector<vector<double>> newA;
    vector<double> newB;
    for (int i = 0; i < lp.m; ++i) {
        if (lp.constrSigns[i] == LE) {
            newA.push_back(sp.A[i]);
            newB.push_back(sp.b[i]);
        } else if (lp.constrSigns[i] == GE) {
            vector<double> row(sp.n);
            for (int j = 0; j < sp.n; ++j) row[j] = -sp.A[i][j];
            newA.push_back(row);
            newB.push_back(-sp.b[i]);
        } else { // EQ
            vector<double> row1(sp.n), row2(sp.n);
            for (int j = 0; j < sp.n; ++j) {
                row1[j] = sp.A[i][j];
                row2[j] = -sp.A[i][j];
            }
            newA.push_back(row1);
            newB.push_back(sp.b[i]);
            newA.push_back(row2);
            newB.push_back(-sp.b[i]);
        }
    }
    sp.m = newA.size();
    sp.A = newA;
    sp.b = newB;

    // Направление: всегда минимизация (max → min через изменение знака c)
    if (lp.sense == MAXIMIZE) {
        for (int j = 0; j < sp.n; ++j) sp.c[j] = -sp.c[j];
    }
    sp.sense = MINIMIZE;
    return sp;
}

vector<double> recoverSolution(const CanonicalProblem& cp, const vector<double>& x_canon) {
    vector<double> x_orig(cp.mapping.size(), 0.0);
    for (size_t j = 0; j < cp.mapping.size(); ++j)
        for (auto& p : cp.mapping[j])
            x_orig[j] += p.second * x_canon[p.first];
    return x_orig;
}

void print_linear_program(const LinearProgram& lp, const string& title) {
    if (!title.empty()) cout << title << "\n";
    cout << "Количество переменных: " << lp.n << "\n";
    cout << "Количество ограничений: " << lp.m << "\n";
    cout << "Знаки переменных:\n";
    for (int j = 0; j < lp.n; ++j) {
        cout << "  x" << j+1 << " : ";
        switch (lp.varSigns[j]) {
            case NONNEGATIVE: cout << ">=0\n"; break;
            case NONPOSITIVE: cout << "<=0\n"; break;
            case FREE:        cout << "свободная\n"; break;
        }
    }
    cout << "Целевая функция: ";
    bool first = true;
    for (int j = 0; j < lp.n; ++j) {
        if (fabs(lp.c[j]) < 1e-12) continue;
        if (!first && lp.c[j] > 0) cout << "+ ";
        else if (!first && lp.c[j] < 0) cout << "- ";
        else if (first && lp.c[j] < 0) cout << "-";
        if (first || fabs(lp.c[j]) != 1.0)
            cout << fabs(lp.c[j]);
        cout << " x" << j+1 << " ";
        first = false;
    }
    if (lp.sense == MINIMIZE) cout << "-> min\n";
    else cout << "-> max\n";

    cout << "Ограничения:\n";
    for (int i = 0; i < lp.m; ++i) {
        first = true;
        for (int j = 0; j < lp.n; ++j) {
            if (fabs(lp.A[i][j]) < 1e-12) continue;
            if (!first && lp.A[i][j] > 0) cout << "+ ";
            else if (!first && lp.A[i][j] < 0) cout << "- ";
            else if (first && lp.A[i][j] < 0) cout << "-";
            if (first || fabs(lp.A[i][j]) != 1.0)
                cout << fabs(lp.A[i][j]);
            cout << " x" << j+1 << " ";
            first = false;
        }
        if (first) cout << "0 ";
        switch (lp.constrSigns[i]) {
            case LE: cout << "<= "; break;
            case GE: cout << ">= "; break;
            case EQ: cout << "= "; break;
        }
        cout << lp.b[i] << "\n";
    }
}

void print_canonical_problem(const CanonicalProblem& cp, const string& title) {
    if (!title.empty()) cout << title << "\n";
    cout << "Количество переменных: " << cp.n << "\n";
    cout << "Количество ограничений: " << cp.m << "\n";
    cout << "Все переменные неотрицательные.\n";
    cout << "Целевая функция: ";
    bool first = true;
    for (int j = 0; j < cp.n; ++j) {
        if (fabs(cp.c[j]) < 1e-12) continue;
        if (!first && cp.c[j] > 0) cout << "+ ";
        else if (!first && cp.c[j] < 0) cout << "- ";
        else if (first && cp.c[j] < 0) cout << "-";
        if (first || fabs(cp.c[j]) != 1.0)
            cout << fabs(cp.c[j]);
        cout << " x" << j+1 << " ";
        first = false;
    }
    if (cp.sense == MINIMIZE) cout << "-> min\n";
    else cout << "-> max\n";

    cout << "Ограничения (все – равенства):\n";
    for (int i = 0; i < cp.m; ++i) {
        first = true;
        for (int j = 0; j < cp.n; ++j) {
            if (fabs(cp.A[i][j]) < 1e-12) continue;
            if (!first && cp.A[i][j] > 0) cout << "+ ";
            else if (!first && cp.A[i][j] < 0) cout << "- ";
            else if (first && cp.A[i][j] < 0) cout << "-";
            if (first || fabs(cp.A[i][j]) != 1.0)
                cout << fabs(cp.A[i][j]);
            cout << " x" << j+1 << " ";
            first = false;
        }
        if (first) cout << "0 ";
        cout << "= " << cp.b[i] << "\n";
    }
}

void print_symmetric_problem(const SymmetricProblem& sp, const string& title) {
    if (!title.empty()) cout << title << "\n";
    cout << "Количество переменных: " << sp.n << "\n";
    cout << "Количество ограничений: " << sp.m << "\n";
    cout << "Все переменные неотрицательные.\n";
    cout << "Целевая функция: ";
    bool first = true;
    for (int j = 0; j < sp.n; ++j) {
        if (fabs(sp.c[j]) < 1e-12) continue;
        if (!first && sp.c[j] > 0) cout << "+ ";
        else if (!first && sp.c[j] < 0) cout << "- ";
        else if (first && sp.c[j] < 0) cout << "-";
        if (first || fabs(sp.c[j]) != 1.0)
            cout << fabs(sp.c[j]);
        cout << " x" << j+1 << " ";
        first = false;
    }
    if (sp.sense == MINIMIZE) cout << "-> min\n";
    else cout << "-> max\n";

    cout << "Ограничения (все – типа <=):\n";
    for (int i = 0; i < sp.m; ++i) {
        first = true;
        for (int j = 0; j < sp.n; ++j) {
            if (fabs(sp.A[i][j]) < 1e-12) continue;
            if (!first && sp.A[i][j] > 0) cout << "+ ";
            else if (!first && sp.A[i][j] < 0) cout << "- ";
            else if (first && sp.A[i][j] < 0) cout << "-";
            if (first || fabs(sp.A[i][j]) != 1.0)
                cout << fabs(sp.A[i][j]);
            cout << " x" << j+1 << " ";
            first = false;
        }
        if (first) cout << "0 ";
        cout << "<= " << sp.b[i] << "\n";
    }
}
