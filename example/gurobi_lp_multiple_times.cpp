#include <lp_wrappers/problem.hpp>
#include <lp_wrappers/gurobi.hpp>
#include <iostream>

int main() {
    LPWrappers::Problem<double> problem(0);
    using Vector = LPWrappers::Problem<double>::Vector;
    std::cin >> problem;
//    std::cin >> problem;

    LPWrappers::GUROBI::Engine<double> solver;
    LPWrappers::GUROBI::Engine<double>::Vector soln;
    auto ret = solver.init(problem, soln);
    std::cout << ret << std::endl << soln << std::endl;
    for(unsigned int i = 0; i < 10000; i++) {
        auto ret = solver.init(problem, soln);
    }

    return 0;
}