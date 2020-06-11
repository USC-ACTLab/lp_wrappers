#ifndef LP_WRAPPERS_CPLEX_HPP
#define LP_WRAPPERS_CPLEX_HPP

#include <lp_wrappers/types.hpp>
#include <lp_wrappers/problem.hpp>
#include <ilcplex/ilocplex.h>

namespace LPWrappers {
namespace CPLEX {
    template<typename T>
    class Engine {
    public:
        using Problem_ = Problem<T>;
        using Vector = typename Problem_::Vector;
        using Index = typename Problem_::Index;

        OptReturnType init(const Problem_& problem, Vector& result) {
            IloEnv env;
            env.setOut(env.getNullStream());

            IloModel model(env);
            IloNumVarArray variables(env);

            for(Index i = 0; i < problem.num_vars(); i++) {
                IloNumVar var(env, problem.lbx()(i), problem.ubx()(i), ILOFLOAT);
                variables.add(var);
            }

            for(Index i = 0; i < problem.num_constraints(); i++) {
                IloExpr expr(env);
                for(Index j = 0; j < problem.num_vars(); j++) {
                    expr += problem.A()(i, j) * variables[j];
                }
                IloRange range(env, problem.lb()(i), expr, problem.ub()(i));
                model.add(range);
            }

            IloExpr linear_cost(env);
            for(Index i = 0; i < problem.num_vars(); i++) {
                linear_cost += variables[i] * problem.c()(i);
            }

            IloObjective obj(env, linear_cost, IloObjective::Minimize);
            model.add(obj);

            IloCplex cplex(model);
            cplex.setOut(env.getNullStream());
            cplex.setWarning(env.getNullStream());

            cplex.solve();

            auto status = cplex.getStatus();

            if(status == IloAlgorithm::Status::Unknown) {
                env.end();
                return OptReturnType::Unknown;
            } else if(status == IloAlgorithm::Status::Feasible) {
                loadResult(env, cplex, variables, result);
                env.end();
                return OptReturnType::Feasible;
            } else if(status == IloAlgorithm::Status::Optimal) {
                loadResult(env, cplex, variables, result);
                env.end();
                return OptReturnType::Optimal;
            } else if(status == IloAlgorithm::Status::Infeasible) {
                env.end();
                return OptReturnType::Infeasible;
            } else if(status == IloAlgorithm::Status::Unbounded) {
                loadResult(env, cplex, variables, result);
                env.end();
                return OptReturnType::Unbounded;
            } else if(status == IloAlgorithm::Status::InfeasibleOrUnbounded) {
                env.end();
                return OptReturnType::InfeasibleOrUnbounded;
            } else if(status == IloAlgorithm::Status::Error) {
                env.end();
                return OptReturnType::Error;
            }

            env.end();
            return OptReturnType::Unknown;

        }

        OptReturnType next(const Problem_& problem, Vector& result) {
            return init(problem, result);
        }

        OptReturnType next(const Problem_& problem, Vector& result,
                const Vector& initial_guess) {
            return init(problem, result);
        }
    private:
        void loadResult(const IloEnv& env, const IloCplex& cplex, const IloNumVarArray& variables, typename Problem<T>::Vector& result) {
            IloNumArray sol(env);
            cplex.getValues(sol, variables);
            result.resize(variables.getSize());
            for(typename Problem<T>::Index i = 0; i < result.rows(); i++) {
                result(i) = sol[i];
            }
        }

    }; // class CPLEX::Engine<T>

} // namespace CPLEX
} // namespace LPWrappers

#endif // LP_WRAPPERS_CPLEX_HPP