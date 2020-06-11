#ifndef LP_WRAPPERS_PROBLEM_HPP
#define LP_WRAPPERS_PROBLEM_HPP

#include <Eigen/Dense>
#include <limits>
#include <absl/strings/str_cat.h>
#include <stdexcept>
#include <iostream>

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
        ubx_mtr(N),
        soft_weights(M),
        soft_convertible(M) {
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

    void reset() {
        A_mtr = Matrix(0, num_vars());
        lb_mtr = Vector(0);
        ub_mtr = Vector(0);
        c_mtr.setZero();
        lbx_mtr.setConstant(num_vars(), minus_infinity);
        ubx_mtr.setConstant(num_vars(), infinity);
    }

    void set_constraint(
            Index constraint_idx,
            const Row& coeff,
            T low,
            T up,
            bool is_soft_convertible = false,
            T soft_weight = T(1)
    ) {
        this->constraint_idx_check(constraint_idx);
        this->vector_size_check(coeff, num_vars());

        A_mtr.row(constraint_idx) = coeff;
        lb_mtr(constraint_idx) = low;
        ub_mtr(constraint_idx) = up;
        soft_convertible[constraint_idx] = is_soft_convertible;
        soft_weights(constraint_idx) = soft_weight;
    }

    void add_constraint(
        const Row& coeff,
        T low,
        T up,
        bool is_soft_convertible = false,
        T soft_weight = T(1)
    ) {
        this->vector_size_check(coeff, num_vars());
        A_mtr.conservativeResize(A_mtr.rows() + 1, Eigen::NoChange);
        lb_mtr.conservativeResize(lb_mtr.rows() + 1, Eigen::NoChange);
        ub_mtr.conservativeResize(ub_mtr.rows() + 1, Eigen::NoChange);
        soft_weights.conservativeResize(soft_weights.rows() + 1, Eigen::NoChange);

        A_mtr.row(A_mtr.rows() - 1) = coeff;
        lb_mtr(lb_mtr.rows() - 1) = low;
        ub_mtr(ub_mtr.rows() - 1) = up;
        soft_weights(soft_weights.rows() - 1) = soft_weight;
        soft_convertible.push_back(is_soft_convertible);
    }

    void set_var_limits(Index var_idx, T low, T up) {
        this->var_idx_check(var_idx);

        lbx_mtr(var_idx) = low;
        ubx_mtr(var_idx) = up;
    }

    void add_c(const Vector& c) {
        this->vector_size_check(c, num_vars());
        c_mtr += c;
    }

    void add_c_block(Index i, const Vector& c) {
        this->var_idx_check(c.rows() + i - 1);
        c_mtr.block(i, 0, c.rows(), 1) += c;
    }

    const Matrix& A() const {
        return A_mtr;
    }

    const Vector& lb() const {
        return lb_mtr;
    }

    const Vector& ub() const {
        return ub_mtr;
    }

    const Vector& lbx() const {
        return lbx_mtr;
    }

    const Vector& ubx() const {
        return ubx_mtr;
    }

    const Vector& c() const {
        return c_mtr;
    }

    template<typename S>
    Problem<S> cast() const {
        Problem<S> new_problem(0);
        new_problem.A_mtr = A_mtr.template cast<S>();
        new_problem.c_mtr = c_mtr.template cast<S>();
        new_problem.lb_mtr = lb_mtr.template cast<S>();
        new_problem.ub_mtr = ub_mtr.template cast<S>();
        new_problem.lbx_mtr = lbx_mtr.template cast<S>();
        new_problem.ubx_mtr = ubx_mtr.template cast<S>();
        new_problem.soft_convertible = soft_convertible;
        new_problem.soft_weights = soft_weights.template cast<S>();
        return new_problem;
    }

    bool is_consistent() const {
        for(Index i = 0; i < num_vars(); i++) {
            if(ubx_mtr(i) < lbx_mtr(i)) {
                return false;
            }
        }

        for(Index i = 0; i < num_constraints(); i++) {
            if(ub_mtr(i) < lb_mtr(i)) {
                return false;
            }
        }
    }

    bool verify(const Vector& solution, T tolerance = 0) const {
        this->vector_size_check(solution, num_vars());

        for(Index i = 0; i < num_vars(); i++) {
            if(solution(i) < lbx_mtr(i) - tolerance) {
                return false;
            }

            if(solution(i) > ubx_mtr(i) + tolerance) {
                return false;
            }
        }

        Vector constraint_evals = A_mtr * solution;
        for(Index i = 0; i < num_constraints(); i++) {
            if(constraint_evals(i) < lb_mtr(i) - tolerance) {
                return false;
            }

            if(constraint_evals(i) > ub_mtr(i) + tolerance) {
                return false;
            }
        }

        return true;
    }

    T objective(const Vector& solution) const {
        this->vector_size_check(solution, num_vars());
        return (c_mtr.transpose() * solution)(0, 0);
    }

    template<typename S>
    friend std::ostream& operator<<(std::ostream&, const Problem<S>&);

    template<typename S>
    friend std::istream& operator>>(std::istream&, Problem<S>&);
private:
    Matrix A_mtr;
    Vector c_mtr, lb_mtr, ub_mtr, lbx_mtr, ubx_mtr;

    std::vector<bool> soft_convertible;
    Vector soft_weights;

    void constraint_idx_check(Index cidx) const {
        if(cidx < 0 || cidx >= num_constraints()) {
            throw std::domain_error (
                absl::StrCat(
                    "constraint index out of range. given idx: ",
                    cidx,
                    " num_constraints: ",
                    num_constraints()
                )
            );
        }
    }

    void var_idx_check(Index vidx) const {
        if(vidx < 0 || vidx >= num_vars()) {
            throw std::domain_error (
                absl::StrCat(
                    "var idx out of range. given idx: ",
                    vidx,
                    " num_vars: ",
                    num_vars()
                )
            );
        }
    }

    void vector_size_check(const Row& row, Index sz) const {
        if(row.cols() != sz) {
            throw std::domain_error (
                absl::StrCat(
                    "vector size check failed. given vector size: ",
                    row.cols(),
                    ", required size: ",
                    sz
                )
            );
        }
    }


    void vector_size_check(const Vector& col, Index sz) const {
        if(col.rows() != sz) {
            throw std::domain_error (
                    absl::StrCat(
                        "vector size check failed."
                    )
            );
        }
    }
}; // class Problem

template<typename T>
std::ostream& operator<<(std::ostream& os, const Problem<T>& problem) {
    using Problem_ = Problem<T>;
    using Index = typename Problem_::Index;

    auto before_precision = os.precision();

    os.precision(std::numeric_limits<T>::max_digits10);

    Index n = problem.num_vars();
    Index m = problem.num_constraints();
    os << std::fixed << n << " " << m << std::endl;

    for(Index i = 0; i < n; i++) {
        os << problem.c_mtr(i) << " ";
    }
    os << std::endl;

    for(Index i = 0; i < n; i++) {
        os << problem.lbx_mtr(i) << " ";
    }
    os << std::endl;

    for(Index i = 0; i < n; i++) {
        os << problem.ubx_mtr(i) << " ";
    }

    for(Index i = 0; i < m; i++) {
        for(Index j = 0; j < n; j++) {
            os << problem.A_mtr(i, j) << " ";
        }
        os << std::endl;
    }

    for(Index i = 0; i < m; i++) {
        os << problem.lb_mtr(i) << " ";
    }
    os << std::endl;

    for(Index i = 0; i < m; i++) {
        os << problem.ub_mtr(i) << " ";
    }
    os << std::endl;

    for(Index i = 0; i < m; i++) {
        os << problem.soft_weights(i) << " ";
    }
    os << std::endl;

    for(Index i = 0; i < m; i++) {
        os << problem.soft_convertible[i] << " ";
    }
    os << std::endl;

    os.precision(before_precision);

    return os;
}

template<typename T>
std::istream& operator>>(std::istream& is, Problem<T>& problem) {
    using Problem_ = Problem<T>;
    using Index = typename Problem_::Index;

    Index n, m;
    is >> n >> m;
    problem.c_mtr.resize(n);
    problem.lbx_mtr.resize(n);
    problem.ubx_mtr.resize(n);
    problem.A_mtr.resize(m, n);
    problem.lb_mtr.resize(m);
    problem.ub_mtr.resize(m);
    problem.soft_weights.resize(m);
    problem.soft_convertible.resize(m);

    for(Index i = 0; i < n; i++) {
        is >> problem.c_mtr(i);
    }

    for(Index i = 0; i < n; i++) {
        is >> problem.lbx_mtr(i);
    }

    for(Index i = 0; i < n; i++) {
        is >> problem.ubx_mtr(i);
    }

    for(Index i = 0; i < m; i++) {
        for(Index j = 0; j < n; j++) {
            is >> problem.A_mtr(i, j);
        }
    }

    for(Index i = 0; i < m; i++) {
        is >> problem.lb_mtr(i);
    }

    for(Index i = 0; i < m; i++) {
        is >> problem.ub_mtr(i);
    }

    for(Index i = 0; i < m; i++) {
        is >> problem.soft_weights(i);
    }

    for(Index i = 0; i < m; i++) {
        bool a;
        is >> a;
        problem.soft_convertible[i] = a;
    }

    return is;
}

} // namespace LPWrappers


#endif // LP_WRAPPERS_PROBLEM_HPP