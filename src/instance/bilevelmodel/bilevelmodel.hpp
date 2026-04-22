#ifndef BILEVEL_MODEL_HPP
#define BILEVEL_MODEL_HPP

#include "../../input/input.hpp"

#include <map>
#include <cmath>
#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include <algorithm>

struct BilevelVariable {
    int id;
    std::string name;
    double lb, ub;
    double obj_leader;
    double obj_follower;
    char type;
    
    BilevelVariable(int id, std::string name, double lb, double ub, 
            double obj_leader, double obj_follower, char type): id(id),
                                                            name(name), 
                                                            lb(lb), ub(ub), 
                                                            obj_leader(obj_leader), 
                                                            obj_follower(obj_follower), 
                                                            type(type) {}
};

struct BilevelConstraint {
    std::string name;
    char sense;
    double rhs;
    std::map<int,double> coeffs;
    double bound;

    BilevelConstraint(std::string name, char sense, double rhs): name(name), 
                                                                 sense(sense), 
                                                                 rhs(rhs),bound(0.0) {}
};

struct BilevelModel {
    double inf = 1e20;

    // Leader
    bool leader_obj_min;
    int nb_leader_vars;
    int nb_leader_constrs;
    double leader_lb, leader_ub;
    std::vector<BilevelVariable> leader_vars;
    std::vector<BilevelConstraint> leader_constrs;

    // Follower
    bool follower_obj_min;
    int nb_follower_vars;
    int nb_follower_constrs;
    double follower_lb, follower_ub;
    std::vector<BilevelVariable> follower_vars;
    std::vector<BilevelConstraint> follower_constrs;

    int nb_follower_eq_constrs;

    // alignment
    double alignment;

    BilevelModel(std::string);

    private:
        void computeFollowerConstrsBounds();
        void computeAlignment();
};

#endif