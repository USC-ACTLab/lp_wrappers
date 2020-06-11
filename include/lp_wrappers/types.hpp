#ifndef LP_WRAPPERS_TYPES_HPP
#define LP_WRAPPERS_TYPES_HPP

#include <Eigen/Dense>

namespace LPWrappers {
template<typename T>
using Row = Eigen::Matrix<T, 1, Eigen::Dynamic>;

template<typename T>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

enum class OptReturnType {
    Optimal,
    Feasible,
    Unbounded,
    Infeasible,
    Error,
    Unknown,
    InfeasibleOrUnbounded
};

std::ostream& operator<<(std::ostream& os, OptReturnType ret_type) {
    switch(ret_type) {
        case OptReturnType::Optimal:
            os << "Optimal";
            break;
        case OptReturnType::Feasible:
            os << "Feasible";
            break;
        case OptReturnType::Unbounded:
            os << "Unbounded";
            break;
        case OptReturnType::Infeasible:
            os << "Infeasible";
            break;
        case OptReturnType::Error:
            os << "Error";
            break;
        case OptReturnType::Unknown:
            os << "Unknown";
            break;
        case OptReturnType::InfeasibleOrUnbounded:
            os << "InfeasibleOrUnbounded";
            break;
    }

    return os;
}


}

#endif // LP_WRAPPERS_TYPES_HPP