#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <vector>

class SimplexSolver {
public:
    SimplexSolver(const std::vector<std::vector<double>>& A,
                  const std::vector<double>& b,
                  const std::vector<double>& c,
                  bool min = true,
                  bool verb = true);

    struct Result {
        bool feasible;
        bool bounded;
        double optimalValue;
        std::vector<double> solution;
    };

    Result solve();

private:
    int m;
    int numVars;
    std::vector<std::vector<double>> tableau;
    std::vector<int> basis;
    std::vector<double> original_c;
    bool minimize;
    int iteration;
    bool verbose;

    void printTableau(const std::string& phase);
    int chooseEnteringColumnPhaseI();
    int chooseEnteringColumnPhaseII();
    int chooseLeavingRow(int col);
    void pivot(int row, int col);
    bool phaseI();
    void setupPhaseII(int orig_n);
    bool phaseII();
};

#endif
