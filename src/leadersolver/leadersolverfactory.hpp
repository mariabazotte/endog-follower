#ifndef LEADER_SOLVER_FACTORY_HPP
#define LEADER_SOLVER_FACTORY_HPP

#include "../instance/instance.hpp"
#include "../input/input.hpp"
#include "leadersolver.hpp"
#include "gurobi_c++.h"

class LeaderSolverFactory{
    public:
        inline AbstractLeaderSolver* createLeaderSolver(const Input & input, const Instance & instance, std::string name){
            
            Input::SolverApproach chosenSolverApproach = input.getSolverApproach();
            switch (chosenSolverApproach)
            {   
                case Input::SolverApproach::TR_SAA :{
                    return new SAALeaderSolver(input,instance,name);
                    break;
                }
                case Input::SolverApproach::DEP :{
                    return new LeaderSolver(input,instance,name);
                    break;
                }
                default:{
                    throw std::string("ERROR: Invalid Solver Approach option in command line.");
                    exit(0);
                    break;
                }
            }
            return NULL;
        }
};

#endif