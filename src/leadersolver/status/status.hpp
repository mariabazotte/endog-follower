#ifndef STATUS_HPP
#define STATUS_HPP

#include "gurobi_c++.h"

enum Status{
    Not_Solved = 0,
    Unknown = 1,
    Optimal = 2,
    Time_Limit = 3,
    Node_Limit = 4,
    Iteration_Limit = 5,
    Unbounded = 6,
    Infeasible = 7,
    Infeasible_or_Unbounded = 8
};

Status statusFromGurobi(int gurobi_status);
std::ostream& operator<<(std::ostream& lhs, const Status & status);

#endif