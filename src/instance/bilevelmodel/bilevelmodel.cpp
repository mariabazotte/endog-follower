#include "bilevelmodel.hpp"

BilevelModel::BilevelModel(std::string instance_file) {
    // Read bilevel model from file.
    std::ifstream file(instance_file.c_str());
    if(!file.fail()){ 
        nb_follower_eq_constrs = 0;

        file >> nb_leader_vars;
        file >> nb_leader_constrs;
        file >> leader_lb;
        file >> leader_ub;

        file >> nb_follower_vars;
        file >> nb_follower_constrs;
        file >> follower_lb;
        file >> follower_ub;

        // Reading Variables
        std::string name; char ty;
        double lb, ub, obju, objl;
        if(nb_leader_vars > 0){
            for(int id = 0; id < nb_leader_vars; ++id){
                file >> name;
                file >> lb;
                file >> ub;
                file >> obju;
                file >> ty;
                if(ty != 'C' && ty != 'B' && ty != 'I') throw std::runtime_error(std::string("Incorrect instance file. Leader variables (C,B,I)")); 
                leader_vars.push_back(BilevelVariable(id,name,lb,ub,obju,0.0,ty));
            }
        }
        if(nb_follower_vars > 0){
            for(int id = 0; id < nb_follower_vars; ++id){
                file >> name;
                file >> lb;
                file >> ub;
                file >> obju;
                file >> objl;
                file >> ty;
                if(ty != 'C') throw std::runtime_error(std::string("Incorrect instance file. Follower variables (C)")); 
                follower_vars.push_back(BilevelVariable(id+nb_leader_vars,name,lb,ub,obju,objl,ty));
            }
        }

        // Reading Constraints
        std::string name_constr;
        char sense;
        double rhs;
        double coeff;
        if(nb_leader_constrs > 0){
            for(int i = 0; i < nb_leader_constrs; ++i){
                file >> name_constr;
                file >> sense;
                file >> rhs;
                if(ty != '=' && ty != '<' && ty != '>') 
                    throw std::runtime_error(std::string("Incorrect instance file. Leader constraints (=,<,>)")); 
                leader_constrs.push_back(BilevelConstraint(name_constr,sense,rhs));

                for(int id = 0; id < nb_leader_vars; ++id) { // only leader variables in these constraints
                    file >> coeff;
                    leader_constrs[i].coeffs[id] = coeff;
                }
            }
        }

        for(int i = 0; i < nb_follower_constrs; ++i){
            file >> name_constr;
            file >> sense;
            file >> rhs;
            if(sense != '=' && sense != '<' && sense != '>')  
                throw std::runtime_error(std::string("Incorrect instance file. Follower constraints (=,<,>)")); 
            follower_constrs.push_back(BilevelConstraint(name_constr,sense,rhs));
            if(sense == '=') nb_follower_eq_constrs += 1;

            for(int id = 0; id < nb_leader_vars; ++id) {
                file >> coeff;
                follower_constrs[i].coeffs[id] = coeff;
            }
            for(int id = 0; id < nb_follower_vars; ++id) {
                file >> coeff;
                follower_constrs[i].coeffs[id+nb_leader_vars] = coeff;
            }
        }
        file.close();
    }else{
        std::string errorMessage = std::string("Error: Could not open file ") + instance_file + std::string("\n");
        throw std::runtime_error(errorMessage);
        exit(0);
    }

    // Verify bilevel model.
    // Leader problem is not bounded: problem defining general decision-dependent case.
    if(leader_lb <= -inf || leader_ub >= inf){
        throw std::runtime_error("Leader problem is not bounded.");
        exit(0);
    }
    // Follower problem is not bounded: problem defining strong-weak and general decision-dependent case.
    if(follower_lb <= -inf || follower_ub >= inf){
        throw std::runtime_error("Follower problem is not bounded.");
        exit(0);
    }
}