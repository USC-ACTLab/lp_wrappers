#ifndef LP_WRAPPERS_PROBLEM_HPP
#define LP_WRAPPERS_PROBLEM_HPP

#include <Eigen/Dense>
#include <limits>

namespace LPWrappers {

/*
 * Defines a linear program in the form
 *      minimize c^T x
 *      subject to
 *          lb <= Ax <= ub
 *          lbx <= x <= ubx
 */
template<typename T>
class Problem {
public:
    using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    using Row = Eigen::Matrix<T, 1, Eigen::Dynamic>;
    using Index = Eigen::Index;

    static constexpr T minus_infinity = std::numeric_limits<T>::lowest();
    static constexpr T infinity = std::numeric_limits<T>::max();


    Problem(Index N, Index M = 0) :
        A_mtr(M, N),
        c_mtr(N),
        lb_mtr(M),
        ub_mtr(M),
        lbx_mtr(N),
        ubx_mtr(N) {
        A_mtr.setZero();
        c_mtr.setZero();
        lbx_mtr.setConstant(N, minus_infinity);
        ubx_mtr.setConstant(N, infinity);
    }

    inline Index num_vars() const {
        return c_mtr.rows();
    }

    inline Index num_constraints() const {
        return A_mtr.rows();
    }

    inline bool is_lbx_unbounded(Index var_idx) const {
        return lbx_mtr(var_idx) == minus_infinity;
    }

    inline bool is_ubx_unbounded(Index var_idx) const {
        return ubx_mtr(var_idx) == infinity;
    }
private:
    Matrix A_mtr;
    Vector c_mtr, lb_mtr, ub_mtr, lbx_mtr, ubx_mtr;
}; // class Problem

} // namespace LPWrappers


#endif // LP_WRAPPERS_PROBLEM_HPP