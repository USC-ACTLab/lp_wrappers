#include <lp_wrappers/problem.hpp>
#include <lp_wrappers/cplex.hpp>
#include <iostream>

int main() {
    LPWrappers::Problem<double> problem(2);
    using Vector = LPWrappers::Problem<double>::Vector;
    problem.set_var_limits(0, 0, LPWrappers::Problem<double>::infinity);
    problem.set_var_limits(1, 0, LPWrappers::Problem<double>::infinity);
    Vector coeff(2);
    coeff(0) = 1; coeff(1) = 1;
    problem.add_constraint(coeff, LPWrappers::Problem<double>::minus_infinity, 6);
    coeff(0) = 1; coeff(1) = -1;
    problem.add_constraint(coeff, LPWrappers::Problem<double>::minus_infinity, 4);
    Vector c(2);
    c(0) = 2; c(1) = 1;
    problem.add_c(c);

//    std::cin >> problem;

    LPWrappers::CPLEX::Engine<double> solver;
    LPWrappers::CPLEX::Engine<double>::Vector soln;
    auto ret = solver.init(problem, soln);
    std::cout << ret << std::endl;
    if(ret == LPWrappers::OptReturnType::Optimal) {
        std::cout << soln << std::endl;
    }
    return 0;
}