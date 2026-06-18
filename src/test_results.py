import math
import numpy as np
import pandas as pd
import openpyxl
import argparse
import sys
from pandas.errors import EmptyDataError, ParserError

############################################################################################################
#                                              Instances
############################################################################################################

instances = ["denegre-miblp_20_15_50_0110_10_",\
             "denegre-miblp_20_20_50_0110_5_",\
             "denegre-miblp_20_20_50_0110_5_e1_",\
             "denegre-miblp_20_20_50_0110_10_",\
             "denegre-miblp_20_20_50_0110_15_",\
             "denegre-milp_30_30_30_2310_5_",\
             "denegre-milp_30_30_30_2310_5_e1_",\
             "denegre-milp_30_30_30_2310_5_e2_",\
             "denegre-milp_30_30_30_2310_10_",\
             "denegre-milp_30_30_30_2310_10_e1_",\
             "denegre-milp_30_30_30_2310_10_e2_",\
             "xuwang-bmilplib_10_",\
             "xuwang-bmilplib_60_",\
             "xuwang-bmilplib_110_",\
             "xuwang-bmilplib_160_",\
             "xuwang-bmilplib_210_",\
             "xuwang-bmilplib_260_"] 
seeds     = [i for i in range(1,11)]

upd_instances = True

############################################################################################################
#                                              Behaviors
############################################################################################################

# Follower near/optimality
near_opt     = 1
eps_near_opt = 0.01
mul_near_opt = 0

# Follower behaviors
behaviors  = [0, 1, 2] 
spec_types = {0: [0], 1: [0,1,2,3,4,5], 2: [0,1,2,3,4]} 

name_behaviors  = ["fixstrongweak","depstrongweak","depgeneral"]
name_spec_types = { 0: [""], 
                    1: ["proportional","threshold","strong","fragile","strong_power","fragile_power"], 
                    2: ["neutral","proportional","fragile","fragile_power","strong_power"]}

# Parameters for each behavior type.
fix_strwk_coop  = [1.0, 0.7, 0.5, 0.3, 0.0]
dep_strwk_thr = [0.5, 2.0, 5.0, 10.0]
dep_strwk_str = [0.5, 2.0, 5.0, 10.0]
dep_strwk_frg = [0.5, 2.0, 5.0, 10.0]
dep_strwk_str_pow = [2.0, 5.0, 10.0]
dep_strwk_frg_pow = [2.0, 5.0, 10.0]
dep_gen_nbint  = [-1]
dep_gen_scal_frg = [0.5, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0]
dep_gen_scal_frg_pow = [1.0, 2.0, 5.0, 10.0, 20.0]
dep_gen_scal_str_pow = [1.0, 2.0, 5.0, 10.0, 20.0]

behavior_params = {
    0: {0: fix_strwk_coop},
    1: {0: [0.0], 1: dep_strwk_thr, 2: dep_strwk_str, 3: dep_strwk_frg, 4: dep_strwk_str_pow, 5: dep_strwk_frg_pow},
    2: {0: [[0.0, 0.0]], 1: [[i, 1.0] for i in dep_gen_nbint], 2: [[i, s] for i in dep_gen_nbint for s in dep_gen_scal_frg],3: [[i, s] for i in dep_gen_nbint for s in dep_gen_scal_frg_pow], 4: [[i, s] for i in dep_gen_nbint for s in dep_gen_scal_str_pow]}
}

# Complete configuration of behaviors.
behaviors_config = [(bv,tp,pr) for bv in behaviors for tp in spec_types[bv] for pr in behavior_params[bv][tp]]

def get_name_behavior_config(bv_config):
    bv = bv_config[0]
    tp = bv_config[1]
    pr = bv_config[2]
    return name_behaviors[bv] + "_" + name_spec_types[bv][tp] + "_" + str(pr)

############################################################################################################
#                                              Scenarios
############################################################################################################

# SAA problem defined for behaviors 1 (Decision-dependent strong-weak) and 2 (Decision-dependent general).
# We solve behavior 0 (Fixed strong-weak) only with deterministic equivalent formulation.
nb_saa_problems  = {1: 1, 2: 1} 
nb_saa_scenarios = {1: 100, 2: 1000} 
nb_saa_thinning  = {1: 1, 2: 5} 
if upd_instances:
    nb_saa_scenarios[2] = 1500
    nb_saa_thinning[2] = 10

# Evaluation with sampling only defined for behavior 2 (Decision-dependent general).
# For behaviors 1 and 2 we can enumerate the endogenous scenarios and compute the 
# endogenous probability, while for behavior 3 we need to use MCMC (sampling from convex set).
nb_eval_problems  = 4
nb_eval_scenarios = 10000
nb_eval_thinning  = 50
coord_har = 1
metrp_has = 0

# Method dep gen
method_dep_gen = "callback" # "callback", "ls", "sglevel_explstep1_grbminmax0", "penalty", "relax" 

time = 10800

############################################################################################################
#                                              Auxiliars
############################################################################################################

def gettstudent():
    if(len(seeds)==1):
        return 0
    if(len(seeds)==2):
        return 12.71
    if(len(seeds)==3):
        return 4.30
    if(len(seeds)==4):
        return 3.18
    if(len(seeds)==5):
        return 2.78
    if(len(seeds)==6):
        return 2.57
    if(len(seeds)==7):
        return 2.45
    if(len(seeds)==8):
        return 2.36
    if(len(seeds)==7):
        return 2.31
    if(len(seeds)==10):
        return 2.26
    return 0

name_folder = "../results/"
def get_file_name(instance,seed,bv_cfg,sol=""):
    bv = bv_cfg[0]
    tp = bv_cfg[1]
    pr = bv_cfg[2]
    
    # General folders.
    name = ''
    if bv == 0:
        name = name_folder + name_behaviors[bv] + "/"
    if bv == 1:
        name = name_folder + name_behaviors[bv] + "/" + sol + "/" + name_spec_types[bv][tp] + "/"
    if bv == 2:
        name = name_folder + name_behaviors[bv] + "/" + name_spec_types[bv][tp] + "/"

    # Parameters specific behavior.
    if bv == 0:
        name += "param" + f"{pr:.2f}" + "_"
    if bv == 1:
        # Proportional does not have extra parameters.
        if tp != 0:
            name += "param" + f"{pr:.2f}" + "_"
    if bv == 2:
        # Neutral does not have extra parameters.
        if tp != 0:
            name += "nbintv" + str(int(pr[0])) + "_" + "scal" f"{pr[1]:.2f}" + "_"
        name += method_dep_gen + "_"

    # Near optimiality parameters. 
    name += "nearopt" + str(near_opt) 
    if near_opt:
        name += "_epsnearopt" + f"{eps_near_opt:.2f}" 
        name += "_multnearopt" + str(mul_near_opt)

    # Instance parameters.
    name += "_" + instance + str(seed)
    if upd_instances:
        name += "_mod"

    # SAA parameters.
    if (bv == 1 and sol == 'saa'):
        name += "_nbsaapr" + str(nb_saa_problems[bv]) + "_nbsaasc" + str(nb_saa_scenarios[bv]) + "_nbsaath" + str(nb_saa_thinning[bv])
    if (bv == 2):
        name += "_nbsaapr" + str(nb_saa_problems[bv]) + "_nbsaasc" + str(nb_saa_scenarios[bv]) + "_nbsaath" + str(nb_saa_thinning[bv])
    
    if(bv == 1 and sol == 'dep'):
        name += "_errpwl0.000100"
        # name += "_intpwl100"

    # Evaluation parameters.
    if bv == 2:
        name += "_nbvalpr" + str(nb_eval_problems) + "_nbvalsc" + str(nb_eval_scenarios) + "_nbvalth" + str(nb_eval_thinning)
        name +=  "_coordhar" + str(coord_har) + "_metrophas" + str(metrp_has)

    return name

def get_sol_file_name(instance,seed,bv_cfg,sol=""):
    name = get_file_name(instance,seed,bv_cfg,sol)

    # Time parameter.
    name += "_time" + str(time) + "_solution.csv"

    return name

def get_cmp_file_name(instance,seed,bv_cfg,sol=""):
    name = get_file_name(instance,seed,bv_cfg,sol)

    if bv_cfg[0] != 2:
        name += "_nbvalpr" + str(nb_eval_problems) + "_nbvalsc" + str(nb_eval_scenarios) + "_nbvalth" + str(nb_eval_thinning)
        name +=  "_coordhar" + str(coord_har) + "_metrophas" + str(metrp_has)

    # Time parameter.
    name += "_time" + str(time) + "_comparison.csv"

    return name 

############################################################################################################
#                                                 Main
############################################################################################################

def main():
    # General instances.
    for instance in instances:
        # General definition of behaviors.
        for bv in behaviors:
            for tp in spec_types[bv]:
                for pr in behavior_params[bv][tp]:
                    # Seeds for instances.
                    for seed in seeds:
                        if bv == 0 or bv == 2:
                            name_sol_file = get_sol_file_name(instance,seed,(bv,tp,pr))
                            name_cmp_file = get_cmp_file_name(instance,seed,(bv,tp,pr))

                            try:
                                df_sol = pd.read_csv(name_sol_file,delimiter=';')
                                if df_sol.empty:
                                    raise ValueError("File contains only header")
                            except (EmptyDataError, FileNotFoundError, ParserError, ValueError):
                                with open("error" + str(eps_near_opt) + ".txt", "a") as f:
                                    f.write(f"{name_sol_file}\n")

                            try:
                                df_cmp = pd.read_csv(name_cmp_file,delimiter=';')
                                if df_cmp.empty:
                                    raise ValueError("File contains only header")
                            except (EmptyDataError, FileNotFoundError, ParserError, ValueError):
                                with open("error" + str(eps_near_opt) + ".txt", "a") as f:
                                    f.write(f"{name_cmp_file}\n")
                                
                        if bv == 1:
                            name_sol_file_dep = get_sol_file_name(instance,seed,(bv,tp,pr),"dep")
                            name_sol_file_saa = get_sol_file_name(instance,seed,(bv,tp,pr),"saa")
                            name_cmp_file = get_cmp_file_name(instance,seed,(bv,tp,pr),"saa")

                            try:
                                df_sol_dep = pd.read_csv(name_sol_file_dep,delimiter=';')
                                if df_sol_dep.empty:
                                    raise ValueError("File contains only header")
                            except (EmptyDataError, FileNotFoundError, ParserError, ValueError):
                                with open("error" + str(eps_near_opt) + ".txt", "a") as f:
                                    f.write(f"{name_sol_file_dep}\n")

                            try:
                                df_sol_saa = pd.read_csv(name_sol_file_saa,delimiter=';')
                                if df_sol_saa.empty:
                                    raise ValueError("File contains only header")
                            except (EmptyDataError, FileNotFoundError, ParserError, ValueError):
                                with open("error" + str(eps_near_opt) + ".txt", "a") as f:
                                    f.write(f"{name_sol_file_saa}\n")

                            try:
                                df_cmp = pd.read_csv(name_cmp_file,delimiter=';')
                                if df_cmp.empty:
                                    raise ValueError("File contains only header")
                            except (EmptyDataError, FileNotFoundError, ParserError, ValueError):
                                with open("error" + str(eps_near_opt) + ".txt", "a") as f:
                                    f.write(f"{name_cmp_file}\n")



 

if __name__ == "__main__":
    main()

