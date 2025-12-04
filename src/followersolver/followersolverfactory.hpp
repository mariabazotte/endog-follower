#ifndef FOLLOWER_SOLVER_FACTORY_HPP
#define FOLLOWER_SOLVER_FACTORY_HPP

#include "gurobi_c++.h"
#include "../input/input.hpp"
#include "../instance/instance.hpp"

#include "followersolver.hpp"
#include "fixstrongweak/fixstrongweakfollower.hpp"
#include "depstrongweak/saadepstrongweakfollower.hpp"
#include "depstrongweak/dedepstrongweakfollower.hpp"
#include "depgeneral/depgeneralfollower.hpp"

#include "../leadersolver/leadersolver.hpp"

class FollowerSolverFactory{
    public:
        inline AbstractFollowerSolver* createFollowerSolver(const Input & input, Instance & instance, int pr, LeaderSolver *leader){
            
            Input::FollowerBehavior chosenBehavior = input.getFollowerBehavior();
            Input::SolverApproach chosenApproach = input.getSolverApproach();
            switch (chosenBehavior)
            {   
                case Input::FollowerBehavior::FixStrongWeak :{
                    if(chosenApproach != Input::SolverApproach::DEP) {
                        throw std::string("ERROR: Invalid Follower Behavior and Solver Approach options in command line.");
                        exit(0);
                        break;
                    }
                    return new FixStrongWeakFollowerSolver(input,instance,leader);
                    break;
                }
                case Input::FollowerBehavior::DepStrongWeak :{
                    switch (chosenApproach)
                    {
                        case Input::SolverApproach::TR_SAA :{
                            return new SAADepStrWkFollowerSolver(input,instance,leader,pr);
                            break;
                        }
                        case Input::SolverApproach::DEP :{
                            return new DEDepStrWkFollowerSolver(input,instance,leader);
                            break;
                        }
                        default:{
                            throw std::string("ERROR: Invalid Solver Approach option in command line.");
                            exit(0);
                            break;
                        }
                    }
                }
                case Input::FollowerBehavior::DepGeneral :{
                    if(chosenApproach != Input::SolverApproach::TR_SAA){
                        throw std::string("ERROR: Invalid Follower Behavior and Solver Approach options in command line.");
                        exit(0);
                        break;
                    }
                    return new DepGeneralFollowerSolver(input,instance,leader,pr);
                    break;
                }
                default:{
                    throw std::string("ERROR: Invalid Follower Behavior option in command line.");
                    exit(0);
                    break;
                }
            }
            return NULL;
        }
};

#endif