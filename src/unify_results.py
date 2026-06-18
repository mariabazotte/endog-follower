import math
import numpy as np
import pandas as pd
pd.set_option('future.no_silent_downcasting', True)
import openpyxl
import argparse
import sys
from scipy import stats
import warnings
warnings.filterwarnings("error")

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

instance_configs = [
    [[5, 55, 10], ["denegre-miblp_20_15_50_0110_10_"]],\
    [[5, 55, 15], ["denegre-miblp_20_20_50_0110_15_"]],\
    [[10, 110, 10], ["denegre-miblp_20_20_50_0110_10_"]],\
    [[15, 165, 5], ["denegre-miblp_20_20_50_0110_5_", "denegre-miblp_20_20_50_0110_5_e1_"]],\
    [[20, 180, 10], ["denegre-milp_30_30_30_2310_10_", "denegre-milp_30_30_30_2310_10_e1_", "denegre-milp_30_30_30_2310_10_e2_"]],\
    [[25, 225, 5], ["denegre-milp_30_30_30_2310_5_", "denegre-milp_30_30_30_2310_5_e1_", "denegre-milp_30_30_30_2310_5_e2_"]],\
    [[10, 40, 10], ["xuwang-bmilplib_10_"]],\
    [[60, 240, 60], ["xuwang-bmilplib_60_"]],\
    [[110, 440, 110], ["xuwang-bmilplib_110_"]],\
    [[160, 640, 160], ["xuwang-bmilplib_160_"]],\
    [[210, 840, 210], ["xuwang-bmilplib_210_"]],\
    [[260, 1040, 260], ["xuwang-bmilplib_260_"]],\
]

upd_instances = True

name_folder = "../results/"
# 
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
dep_gen_nbint  = [-1.0]
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
method_dep_gen = "callback"

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
    if upd_instances == True:
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
#                                              Dictionaries
############################################################################################################

dict_sol_all = dict()
for bv in behaviors:
    dict_sol_all[bv] = {"instance":[], "seed":[], "align":[], "type":[], "params":[]}
    if bv != 2:
        sol = "_dep"
        for opt in ["status","lb","ub","gap","eval","obj_follower","nb_bnb","time"]:
            dict_sol_all[bv][opt+sol] = []
        if near_opt:
            for opt in ["obj_near_opt_follower","obj_near_pes_follower","obj_near_mean_follower"]:
                dict_sol_all[bv][opt+sol] = []
    if bv != 0:
        sol = "_saa"
        for opt in ["lb_1","ub_1","gap_1","status","lb","std_lb","perc_std_lb","stat_lb","ub","std_ub","perc_std_ub","stat_ub","gap","abs_gap","std_gap","perc_std_gap","stat_gap","obj_follower","nb_bnb","var_nb_bnb","nb_pr_not_opt","time","time_lb","time_eval"]:
            dict_sol_all[bv][opt+sol] = []
        if near_opt:
            for opt in ["obj_near_opt_follower","obj_near_pes_follower","obj_near_mean_follower","obj_near_std_follower"]:
                dict_sol_all[bv][opt+sol] = []
    if bv == 1:
        dict_sol_all[bv]["gap_lb"] = []
        dict_sol_all[bv]["gap_ub"] = []

def init_dict_sol_all(bv,tp,params,instance,seed):
    dict_sol_all[bv]["instance"].append(instance) 
    dict_sol_all[bv]["seed"].append(seed)
    dict_sol_all[bv]["type"].append(name_spec_types[bv][tp])
    dict_sol_all[bv]["params"].append(params)

def add_dict_sol_all_align(bv,df):
    if df.shape[0] < 1:
        dict_sol_all[bv]["align"].append(np.nan)
        return np.nan

    dict_sol_all[bv]["align"].append(float(df["ALIGNMENT"][0]))
    return float(df["ALIGNMENT"][0])

def add_dict_sol_all_dep(bv,df):
    sol = "_dep"

    # Could not find a feasible solution within the time limit.
    if df.shape[0] < 1:
        dict_sol_all[bv]["status"+sol].append("-")
        dict_sol_all[bv]["lb"+sol].append("-")
        dict_sol_all[bv]["ub"+sol].append("-")
        dict_sol_all[bv]["gap"+sol].append(100.0)
        dict_sol_all[bv]["eval"+sol].append("-")
        dict_sol_all[bv]["obj_follower"+sol].append("-")
        if near_opt:
            dict_sol_all[bv]["obj_near_opt_follower"+sol].append("-")
            dict_sol_all[bv]["obj_near_pes_follower"+sol].append("-")
            dict_sol_all[bv]["obj_near_mean_follower"+sol].append("-")
        dict_sol_all[bv]["nb_bnb"+sol].append("-")
        dict_sol_all[bv]["time"+sol].append(float(time))
        return

    dict_sol_all[bv]["status"+sol].append(int(df["STATUS"][0]))
    dict_sol_all[bv]["lb"+sol].append(float(df["LB"][0]))
    dict_sol_all[bv]["ub"+sol].append(float(df["UB"][0]))
    dict_sol_all[bv]["gap"+sol].append(min(100.0,100*float(df["GAP"][0])))
    dict_sol_all[bv]["eval"+sol].append(float(df["EVAL"][0]))
    dict_sol_all[bv]["obj_follower"+sol].append(float(df["OBJ_FOLLOWER"][0]))
    
    if near_opt:
        dict_sol_all[bv]["obj_near_opt_follower"+sol].append(float(df["OBJ_OPT_FOLLOWER"][0]))
        dict_sol_all[bv]["obj_near_pes_follower"+sol].append(float(df["OBJ_PES_FOLLOWER"][0]))
        dict_sol_all[bv]["obj_near_mean_follower"+sol].append(float(df["OBJ_MEAN_FOLLOWER"][0]))
    
    dict_sol_all[bv]["nb_bnb"+sol].append(int(df["NB_BNB"][0]))
    dict_sol_all[bv]["time"+sol].append(float(df["TIME"][0]))

def add_dict_sol_all_saa(bv,df):
    sol = "_saa"

    # Could not find a feasible solution within the time limit.
    if df.shape[0] < 1:
        dict_sol_all[bv]["lb_1"+sol].append("-")
        dict_sol_all[bv]["ub_1"+sol].append("-")
        dict_sol_all[bv]["gap_1"+sol].append(100.0)
        dict_sol_all[bv]["status"+sol].append(3)
        dict_sol_all[bv]["lb"+sol].append("-")
        dict_sol_all[bv]["std_lb"+sol].append("-")
        dict_sol_all[bv]["perc_std_lb"+sol].append("-")
        dict_sol_all[bv]["stat_lb"+sol].append("-")
        dict_sol_all[bv]["ub"+sol].append("-")
        dict_sol_all[bv]["std_ub"+sol].append("-")
        dict_sol_all[bv]["perc_std_ub"+sol].append("-")
        dict_sol_all[bv]["stat_ub"+sol].append("-")
        dict_sol_all[bv]["gap"+sol].append(100.0)
        dict_sol_all[bv]["abs_gap"+sol].append("-")
        dict_sol_all[bv]["std_gap"+sol].append("-")
        dict_sol_all[bv]["perc_std_gap"+sol].append("-")
        dict_sol_all[bv]["stat_gap"+sol].append("-")
        dict_sol_all[bv]["obj_follower"+sol].append("-")
        if near_opt:
            dict_sol_all[bv]["obj_near_opt_follower"+sol].append("-")
            dict_sol_all[bv]["obj_near_pes_follower"+sol].append("-")
            dict_sol_all[bv]["obj_near_mean_follower"+sol].append("-")
            dict_sol_all[bv]["obj_near_std_follower"+sol].append("-")
        dict_sol_all[bv]["nb_bnb"+sol].append("-")
        dict_sol_all[bv]["var_nb_bnb"+sol].append("-")
        dict_sol_all[bv]["nb_pr_not_opt"+sol].append("-")
        dict_sol_all[bv]["time"+sol].append(float(time))
        dict_sol_all[bv]["time_lb"+sol].append(float(time))
        dict_sol_all[bv]["time_eval"+sol].append("-")
        return

    dict_sol_all[bv]["lb_1"+sol].append(float(df["LB_1"][0]))
    dict_sol_all[bv]["ub_1"+sol].append(float(df["UB_1"][0]))
    dict_sol_all[bv]["gap_1"+sol].append(min(100.0,100*float(df["GAP_1"][0])))

    dict_sol_all[bv]["status"+sol].append(int(df["STATUS"][0]))

    dict_sol_all[bv]["lb"+sol].append(float(df["LB"][0]))
    dict_sol_all[bv]["std_lb"+sol].append(math.sqrt(float(df["VAR_LB"][0])))
    dict_sol_all[bv]["perc_std_lb"+sol].append(math.sqrt(float(df["VAR_LB"][0]))/math.fabs(float(df["LB"][0])))
    dict_sol_all[bv]["stat_lb"+sol].append(float(df["STAT_LB"][0]))

    dict_sol_all[bv]["ub"+sol].append(float(df["UB"][0]))
    dict_sol_all[bv]["std_ub"+sol].append(math.sqrt(float(df["VAR_UB"][0])))
    dict_sol_all[bv]["perc_std_ub"+sol].append(math.sqrt(float(df["VAR_UB"][0]))/math.fabs(float(df["UB"][0])))
    dict_sol_all[bv]["stat_ub"+sol].append(float(df["STAT_UB"][0]))

    # dict_sol_all[bv]["gap"+sol].append(100*(float(df["GAP"][0])*float(df["UB"][0]))/math.fabs(float(df["UB"][0])))
    dict_sol_all[bv]["gap"+sol].append(min(100.0,100*(float(df["GAP"][0]))))
    dict_sol_all[bv]["abs_gap"+sol].append(float(df["GAP"][0])*float(df["UB"][0]))
    dict_sol_all[bv]["std_gap"+sol].append(math.sqrt(float(df["VAR_GAP"][0])))
    dict_sol_all[bv]["perc_std_gap"+sol].append(math.sqrt(float(df["VAR_GAP"][0]))/max(0.001,(float(df["GAP"][0])*float(df["UB"][0]))))
    dict_sol_all[bv]["stat_gap"+sol].append(float(df["STAT_GAP"][0]))

    dict_sol_all[bv]["obj_follower"+sol].append(float(df["OBJ_FOLLOWER"][0]))

    if near_opt:
        dict_sol_all[bv]["obj_near_opt_follower"+sol].append(float(df["OBJ_OPT_FOLLOWER"][0]))
        dict_sol_all[bv]["obj_near_pes_follower"+sol].append(float(df["OBJ_PES_FOLLOWER"][0]))
        dict_sol_all[bv]["obj_near_mean_follower"+sol].append(float(df["OBJ_MEAN_FOLLOWER"][0]))
        dict_sol_all[bv]["obj_near_std_follower"+sol].append(math.sqrt(float(df["OBJ_VAR_FOLLOWER"][0])))

    dict_sol_all[bv]["nb_bnb"+sol].append(int(df["NB_BNB"][0]))
    dict_sol_all[bv]["var_nb_bnb"+sol].append(float(df["VAR_NB_BNB"][0]))
    dict_sol_all[bv]["nb_pr_not_opt"+sol].append(int(df["NB_NOT_OPT"][0]))

    dict_sol_all[bv]["time"+sol].append(float(df["TIME"][0]))
    dict_sol_all[bv]["time_lb"+sol].append(min(float(time),float(df["TIME_LB"][0])))
    dict_sol_all[bv]["time_eval"+sol].append(float(df["TIME_EVAL"][0]))
    
def add_dict_sol_all_comp(bv):
    if bv != 1:
        raise ValueError("This function computes information available only for behavior 1 (Decision-dependent strong-weak).") 

    if dict_sol_all[bv]["lb_dep"][-1] != "-" and dict_sol_all[bv]["lb_saa"][-1] != "-":
        dict_sol_all[bv]["gap_lb"].append(np.clip(100*(dict_sol_all[bv]["lb_dep"][-1] - dict_sol_all[bv]["lb_saa"][-1])/math.fabs(dict_sol_all[bv]["lb_dep"][-1]),-100,100))
    else:
        dict_sol_all[bv]["gap_lb"].append("-")
    
    if dict_sol_all[bv]["eval_dep"][-1] != "-" and dict_sol_all[bv]["ub_saa"][-1] != "-":
        if dict_sol_all[bv]["ub_dep"][-1] >= 1e100:
            dict_sol_all[bv]["gap_ub"].append(100.0)
        else:
            dict_sol_all[bv]["gap_ub"].append(np.clip(100*(dict_sol_all[bv]["eval_dep"][-1] - dict_sol_all[bv]["ub_saa"][-1])/math.fabs(dict_sol_all[bv]["eval_dep"][-1]),-100.0,100.0))
    else:
        dict_sol_all[bv]["gap_ub"].append(100.0)

dict_sol_mean = dict()
for bv in behaviors:
    dict_sol_mean[bv] = {"instance":[], "type":[], "params":[]}
    if bv != 2:
        for opt in ["gap_dep","nb_bnb_dep","time_dep"]:
            dict_sol_mean[bv][opt] = []
            dict_sol_mean[bv]["int"+opt] = []
    if bv != 0:
        for opt in ["gap_1_saa","gap_saa","nb_bnb_saa","time_saa","time_lb_saa","time_eval_saa","nb_pr_not_opt_saa","perc_std_lb_saa","perc_std_ub_saa","perc_std_gap_saa"]:
            dict_sol_mean[bv][opt] = []
            dict_sol_mean[bv]["int"+opt] = []
        for opt in ["perc_std_lb_saa","perc_std_ub_saa","perc_std_gap_saa"]:
            dict_sol_mean[bv]["max"+opt] = []
    if bv == 1:
        for opt in ["gap_lb","gap_ub"]:
            dict_sol_mean[bv][opt] = []
            dict_sol_mean[bv]["int"+opt] = []

def add_dict_sol_mean(bv,tp,params,instance):
    dict_sol_mean[bv]["instance"].append(instance) 
    dict_sol_mean[bv]["type"].append(name_spec_types[bv][tp])
    dict_sol_mean[bv]["params"].append(params)
    if bv != 2: 
        for opt in ["gap_dep","nb_bnb_dep","time_dep"]:
            values = dict_sol_all[bv][opt][-len(seeds):]
            values = [float(v) for v in values if v != "-"]
            if len(values) == 0:
                dict_sol_mean[bv][opt].append("-")
                dict_sol_mean[bv]["int"+opt].append("-")
            else:
                dict_sol_mean[bv][opt].append(np.mean(values))
                dict_sol_mean[bv]["int"+opt].append(str(round(dict_sol_mean[bv][opt][-1],2)) + "+/-" + str(round((gettstudent()*np.std(values,ddof=1))/math.sqrt(len(seeds)),2)))
    if bv != 0:
        for opt in ["gap_1_saa","gap_saa","nb_bnb_saa","time_saa","time_lb_saa","time_eval_saa","nb_pr_not_opt_saa","perc_std_lb_saa","perc_std_ub_saa","perc_std_gap_saa"]:
            values = dict_sol_all[bv][opt][-len(seeds):]
            values = [float(v) for v in values if v != "-"]
            dict_sol_mean[bv][opt].append(np.mean(values))
            dict_sol_mean[bv]["int"+opt].append(str(round(dict_sol_mean[bv][opt][-1],2)) + "+/-" + str(round((gettstudent()*np.std(values,ddof=1))/math.sqrt(len(seeds)),2)))

        for opt in ["perc_std_lb_saa","perc_std_ub_saa","perc_std_gap_saa"]:
            values = dict_sol_all[bv][opt][-len(seeds):]
            values = [float(v) for v in values if v != "-"]
            dict_sol_mean[bv]["max"+opt].append(np.max(values))
    if bv == 1:
        for opt in ["gap_lb","gap_ub"]:
            values = dict_sol_all[bv][opt][-len(seeds):]
            values = [float(v) for v in values if v != "-"]
            if len(values) == 0:
                dict_sol_mean[bv][opt].append("-")
                dict_sol_mean[bv]["int"+opt].append("-")
            else:
                dict_sol_mean[bv][opt].append(np.mean(values))
                dict_sol_mean[bv]["int"+opt].append(str(round(dict_sol_mean[bv][opt][-1],2)) + "+/-" + str(round((gettstudent()*np.std(values,ddof=1))/math.sqrt(len(seeds)),2)))

def get_interval(df):
    mean = df.mean()
    margin = stats.sem(df)*stats.t.ppf(0.975, df=df.shape[0]-1)
    return str(round(mean,1)) + '+/-' + str(round(margin,1))

def get_dict_sol_mean_per_type():
    bv = 1

    dict_sol_mean_per_type = {"type":[],"params":[],"mean_gap_ub":[]}

    # Convert to DataFrame
    df = pd.DataFrame(dict_sol_all[bv]) 

    # Iterate per behavior type-param
    for tp in spec_types[bv]:
        for pr in behavior_params[bv][tp]:
            dict_sol_mean_per_type["type"].append(name_spec_types[bv][tp])
            dict_sol_mean_per_type["params"].append(pr)
            filt_df = df[(df["type"] == name_spec_types[bv][tp]) & (df["params"] == pr)]
            filt_df = filt_df[~(filt_df == "-").any(axis=1)]
            filt_df["gap_ub"] = pd.to_numeric(filt_df["gap_ub"])
            if filt_df["gap_ub"].dropna().empty:
                dict_sol_mean_per_type["mean_gap_ub"].append("-")
            else:
                dict_sol_mean_per_type["mean_gap_ub"].append(get_interval(filt_df["gap_ub"]))
    return dict_sol_mean_per_type

def get_dict_sol_mean_final(bv):
    dict_sol_mean_final = {"nx":[],"nx_bin":[],"ny":[]}
    if bv == 0:
        for pr in fix_strwk_coop:
            dict_sol_mean_final["type" + str(pr) + "_gap"] = []
            dict_sol_mean_final["type" + str(pr) + "_time"] = []
    
    if bv == 1:
        for tp in spec_types[bv]:
            dict_sol_mean_final["type" + str(int(tp)) + "_saa_gap"] = []
            dict_sol_mean_final["type" + str(int(tp)) + "_saa_time"] = []

            dict_sol_mean_final["type" + str(int(tp)) + "_dep_gap"] = []
            dict_sol_mean_final["type" + str(int(tp)) + "_dep_time"] = []
    
    if bv == 2:
        for tp in spec_types[bv]:
            dict_sol_mean_final["type" + str(int(tp)) + "_gap"] = []
            dict_sol_mean_final["type" + str(int(tp)) + "_time"] = []
        
        for tp in spec_types[bv]:
            dict_sol_mean_final["type" + str(int(tp)) + "_time_eval"] = []

    # Convert to DataFrame
    df = pd.DataFrame(dict_sol_all[bv]) 

    for instance in instance_configs:
        key = instance[0]
        value = instance[1]
        dict_sol_mean_final["nx"].append(key[0])
        dict_sol_mean_final["nx_bin"].append(key[1])
        dict_sol_mean_final["ny"].append(key[2])
        if bv == 0:
            for pr in fix_strwk_coop:
                filt_df = df[(df["params"] == pr) & (df["instance"].isin(value))]
                filt_df = filt_df.replace("-", np.nan).infer_objects(copy=False)

                filt_df["gap_dep"] = filt_df["gap_dep"].apply(pd.to_numeric)
                filt_df["time_dep"] = filt_df["time_dep"].apply(pd.to_numeric)

                dict_sol_mean_final["type" + str(pr) + "_gap"].append(get_interval(filt_df["gap_dep"]))
                dict_sol_mean_final["type" + str(pr) + "_time"].append(get_interval(filt_df["time_dep"]))

        if bv == 1:
            for tp in spec_types[bv]:
                filt_df = df[(df["type"] == name_spec_types[bv][tp]) & (df["instance"].isin(value))]
                
                filt_df = filt_df.replace("-", np.nan)
                filt_df = filt_df.infer_objects(copy=False)

                filt_df["gap_1_saa"] = filt_df["gap_1_saa"].apply(pd.to_numeric, errors='coerce')
                filt_df["time_lb_saa"] = filt_df["time_lb_saa"].apply(pd.to_numeric)
                filt_df["gap_dep"] = filt_df["gap_dep"].apply(pd.to_numeric)
                filt_df["time_dep"] = filt_df["time_dep"].apply(pd.to_numeric)

                dict_sol_mean_final["type" + str(tp) + "_saa_gap"].append(get_interval(filt_df["gap_1_saa"]))
                dict_sol_mean_final["type" + str(tp) + "_saa_time"].append(get_interval(filt_df["time_lb_saa"]))
                dict_sol_mean_final["type" + str(tp) + "_dep_gap"].append(get_interval(filt_df["gap_dep"]))
                dict_sol_mean_final["type" + str(tp) + "_dep_time"].append(get_interval(filt_df["time_dep"]))

        if bv == 2:
            for tp in spec_types[bv]:
                filt_df = df[(df["type"] == name_spec_types[bv][tp]) & (df["instance"].isin(value))]

                filt_df = filt_df.replace("-", np.nan).infer_objects(copy=False)
                filt_df["gap_1_saa"] = filt_df["gap_1_saa"].apply(pd.to_numeric, errors='coerce')
                filt_df["time_lb_saa"] = filt_df["time_lb_saa"].apply(pd.to_numeric, errors='coerce')
                filt_df["time_eval_saa"] = filt_df["time_eval_saa"].apply(pd.to_numeric, errors='coerce')

                dict_sol_mean_final["type" + str(tp) + "_gap"].append(get_interval(filt_df["gap_1_saa"]))
                dict_sol_mean_final["type" + str(tp) + "_time"].append(get_interval(filt_df["time_lb_saa"]))
                dict_sol_mean_final["type" + str(tp) + "_time_eval"].append(get_interval(filt_df["time_eval_saa"]))

    return dict_sol_mean_final

dict_sol_values_all = {"instance":[], "seed":[], "align":[]}
for cf in range(len(behaviors_config)):
    name = get_name_behavior_config(behaviors_config[cf])
    dict_sol_values_all["leader_obj_"+name] = []
for cf in range(len(behaviors_config)):
    name = get_name_behavior_config(behaviors_config[cf])
    dict_sol_values_all["follower_obj_"+name] = []
if near_opt:
    for cf in range(len(behaviors_config)):
        name = get_name_behavior_config(behaviors_config[cf])
        dict_sol_values_all["follower_obj_near_mean_"+name] = []

max_gap_sol_values = 5.0
def add_dict_sol_values_all(instance):
    nb_seeds = len(seeds)
    for sd in range(nb_seeds):
        #  Optimistic configuration is the first one here.
        spec_behavior0 = [(tp,pr) for tp in spec_types[0] for pr in behavior_params[0][tp]]
        nb_cfgs_bv0 = len(spec_behavior0)

        # Only verify for the instances for which the optimistic approach was solved with a gap inferior to 1%.
        if dict_sol_all[0]["gap_dep"][-(nb_cfgs_bv0*nb_seeds)+sd] <= 1.0:
            dict_sol_values_all["instance"].append(instance)
            dict_sol_values_all["seed"].append(seeds[sd])
            dict_sol_values_all["align"].append(dict_sol_all[0]["align"][-(nb_cfgs_bv0*nb_seeds)+sd])

            # Compute percentage compared to optimistic for each behavior.
            for bv in behaviors:
                spec_behaviors = [(tp,pr) for tp in spec_types[bv] for pr in behavior_params[bv][tp]]
                nb_cfgs = len(spec_behaviors)
                
                for cf in range(nb_cfgs):
                    name = get_name_behavior_config((bv,spec_behaviors[cf][0],spec_behaviors[cf][1]))
                    if bv == 0:
                        if dict_sol_all[bv]["gap_dep"][-((nb_cfgs-cf)*nb_seeds)+sd] <= max_gap_sol_values:
                            dict_sol_values_all["leader_obj_"+name].append(min(100.0,100*(dict_sol_all[bv]["ub_dep"][-((nb_cfgs-cf)*nb_seeds)+sd]-dict_sol_all[0]["ub_dep"][-(nb_cfgs_bv0*nb_seeds)+sd])/math.fabs(dict_sol_all[0]["ub_dep"][-(nb_cfgs_bv0*nb_seeds)+sd])))
                        else:
                            dict_sol_values_all["leader_obj_"+name].append("-")
                    else:
                        if dict_sol_all[bv]["gap_1_saa"][-((nb_cfgs-cf)*nb_seeds)+sd] <= max_gap_sol_values:
                            dict_sol_values_all["leader_obj_"+name].append(min(100.0,100*(dict_sol_all[bv]["ub_saa"][-((nb_cfgs-cf)*nb_seeds)+sd]-dict_sol_all[0]["ub_dep"][-(nb_cfgs_bv0*nb_seeds)+sd])/math.fabs(dict_sol_all[0]["ub_dep"][-(nb_cfgs_bv0*nb_seeds)+sd])))
                        else:
                            dict_sol_values_all["leader_obj_"+name].append("-")

                for cf in range(nb_cfgs):
                    name = get_name_behavior_config((bv,spec_behaviors[cf][0],spec_behaviors[cf][1]))
                    if bv == 0:
                        if dict_sol_all[bv]["gap_dep"][-((nb_cfgs-cf)*nb_seeds)+sd] <= max_gap_sol_values:
                            if (math.fabs(dict_sol_all[0]["obj_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]) >= -1e-12) and (math.fabs(dict_sol_all[0]["obj_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]) <= 1e-12):
                                dict_sol_values_all["follower_obj_"+name].append(np.clip((dict_sol_all[bv]["obj_follower_dep"][-((nb_cfgs-cf)*nb_seeds)+sd]-dict_sol_all[0]["obj_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]),-100.0,100.0))
                            else:
                                dict_sol_values_all["follower_obj_"+name].append(np.clip(100*(dict_sol_all[bv]["obj_follower_dep"][-((nb_cfgs-cf)*nb_seeds)+sd]-dict_sol_all[0]["obj_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd])/math.fabs(dict_sol_all[0]["obj_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]),-100.0,100.0))
                        else:
                            dict_sol_values_all["follower_obj_"+name].append("-")
                    else:
                        if dict_sol_all[bv]["gap_1_saa"][-((nb_cfgs-cf)*nb_seeds)+sd] <= max_gap_sol_values:
                            if (math.fabs(dict_sol_all[0]["obj_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]) >= -1e-12) and (math.fabs(dict_sol_all[0]["obj_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]) <= 1e-12):
                                dict_sol_values_all["follower_obj_"+name].append(np.clip((dict_sol_all[bv]["obj_follower_saa"][-((nb_cfgs-cf)*nb_seeds)+sd]-dict_sol_all[0]["obj_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]),-100.0,100.0))
                            else:
                                dict_sol_values_all["follower_obj_"+name].append(np.clip(100*(dict_sol_all[bv]["obj_follower_saa"][-((nb_cfgs-cf)*nb_seeds)+sd]-dict_sol_all[0]["obj_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd])/math.fabs(dict_sol_all[0]["obj_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]),-100.0,100.0))
                        else:
                            dict_sol_values_all["follower_obj_"+name].append("-")

                if near_opt:
                    for cf in range(nb_cfgs):
                        name = get_name_behavior_config((bv,spec_behaviors[cf][0],spec_behaviors[cf][1]))
                        if bv == 0:
                            if dict_sol_all[bv]["gap_dep"][-((nb_cfgs-cf)*nb_seeds)+sd] <= max_gap_sol_values:
                                if (math.fabs(dict_sol_all[0]["obj_near_mean_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]) >= -1e-12) and (math.fabs(dict_sol_all[0]["obj_near_mean_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]) <= 1e-12):
                                    dict_sol_values_all["follower_obj_near_mean_"+name].append(np.clip(dict_sol_all[bv]["obj_near_mean_follower_dep"][-((nb_cfgs-cf)*nb_seeds)+sd]-dict_sol_all[0]["obj_near_mean_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd],-100.0,100.0))
                                else:
                                    dict_sol_values_all["follower_obj_near_mean_"+name].append(np.clip(100*(dict_sol_all[bv]["obj_near_mean_follower_dep"][-((nb_cfgs-cf)*nb_seeds)+sd]-dict_sol_all[0]["obj_near_mean_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd])/math.fabs(dict_sol_all[0]["obj_near_mean_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]),-100.0,100.0))
                            else:
                                dict_sol_values_all["follower_obj_near_mean_"+name].append("-")
                        else:
                            if dict_sol_all[bv]["gap_1_saa"][-((nb_cfgs-cf)*nb_seeds)+sd] <= max_gap_sol_values:
                                if (math.fabs(dict_sol_all[0]["obj_near_mean_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]) >= -1e-12) and (math.fabs(dict_sol_all[0]["obj_near_mean_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]) <= 1e-12):
                                    dict_sol_values_all["follower_obj_near_mean_"+name].append(np.clip((dict_sol_all[bv]["obj_near_mean_follower_saa"][-((nb_cfgs-cf)*nb_seeds)+sd]-dict_sol_all[0]["obj_near_mean_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]),-100.0,100.0))
                                else:
                                    dict_sol_values_all["follower_obj_near_mean_"+name].append(np.clip(100*(dict_sol_all[bv]["obj_near_mean_follower_saa"][-((nb_cfgs-cf)*nb_seeds)+sd]-dict_sol_all[0]["obj_near_mean_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd])/math.fabs(dict_sol_all[0]["obj_near_mean_follower_dep"][-(nb_cfgs_bv0*nb_seeds)+sd]),-100,100.0))
                            else:
                                dict_sol_values_all["follower_obj_near_mean_"+name].append("-")

dict_sol_values_mean = {"instance":[]}
for cf in range(len(behaviors_config)):
    name = get_name_behavior_config(behaviors_config[cf])
    dict_sol_values_mean["leader_obj_"+name] = []
    dict_sol_values_mean["int"+"leader_obj_"+name] = []
for cf in range(len(behaviors_config)):
    name = get_name_behavior_config(behaviors_config[cf])
    dict_sol_values_mean["follower_obj_"+name] = []
    dict_sol_values_mean["int"+"follower_obj_"+name] = []
if near_opt:
    for cf in range(len(behaviors_config)):
        name = get_name_behavior_config(behaviors_config[cf])
        dict_sol_values_mean["follower_obj_near_mean_"+name] = []
        dict_sol_values_mean["int"+"follower_obj_near_mean_"+name] = []

def add_dict_sol_values_mean(instance):
    dict_sol_values_mean["instance"].append(instance)

    nb_cfgs = len(behaviors_config)

    for cf in range(nb_cfgs):
        opt = "leader_obj_" + get_name_behavior_config(behaviors_config[cf])

        values = dict_sol_values_all[opt][-len(seeds):]
        values = [float(v) for v in values if v != "-"]

        if len(values) == 0:
            dict_sol_values_mean[opt].append("-")
            dict_sol_values_mean["int"+opt].append("-")
        elif len(values) == 1:
            dict_sol_values_mean[opt].append(np.mean(values))
            dict_sol_values_mean["int"+opt].append(str(round(dict_sol_values_mean[opt][-1],2)) + "+/- 0.00")
        else:
            dict_sol_values_mean[opt].append(np.mean(values))
            dict_sol_values_mean["int"+opt].append(str(round(dict_sol_values_mean[opt][-1],2)) + "+/-" + str(round((gettstudent()*np.std(values,ddof=1))/math.sqrt(len(seeds)),2)))

    for cf in range(nb_cfgs):
        opt = "follower_obj_" + get_name_behavior_config(behaviors_config[cf])

        values = dict_sol_values_all[opt][-len(seeds):]
        values = [float(v) for v in values if v != "-"]

        if len(values) == 0:
            dict_sol_values_mean[opt].append("-")
            dict_sol_values_mean["int"+opt].append("-")
        elif len(values) == 1:
            dict_sol_values_mean[opt].append(np.mean(values))
            dict_sol_values_mean["int"+opt].append(str(round(dict_sol_values_mean[opt][-1],2)) + "+/- 0.00")
        else:
            dict_sol_values_mean[opt].append(np.mean(values))
            dict_sol_values_mean["int"+opt].append(str(round(dict_sol_values_mean[opt][-1],2)) + "+/-" + str(round((gettstudent()*np.std(values,ddof=1))/math.sqrt(len(seeds)),2)))

    if near_opt:
        for cf in range(nb_cfgs):
            opt = "follower_obj_near_mean_" + get_name_behavior_config(behaviors_config[cf])

            values = dict_sol_values_all[opt][-len(seeds):]
            values = [float(v) for v in values if v != "-"]

            if len(values) == 0:
                dict_sol_values_mean[opt].append("-")
                dict_sol_values_mean["int"+opt].append("-")
            elif len(values) == 1:
                dict_sol_values_mean[opt].append(np.mean(values))
                dict_sol_values_mean["int"+opt].append(str(round(dict_sol_values_mean[opt][-1],2)) + "+/- 0.00")
            else:
                dict_sol_values_mean[opt].append(np.mean(values))
                dict_sol_values_mean["int"+opt].append(str(round(dict_sol_values_mean[opt][-1],2)) + "+/-" + str(round((gettstudent()*np.std(values,ddof=1))/math.sqrt(len(seeds)),2)))

dict_comp_all = {"instance":[], "seed":[], "align":[], "config":[]}
for cf in range(len(behaviors_config)):
    column = get_name_behavior_config(behaviors_config[cf])
    dict_comp_all["leader_obj_mean_"+column] = []
for cf in range(len(behaviors_config)):
    column = get_name_behavior_config(behaviors_config[cf])
    dict_comp_all["follower_obj_"+column] = []
if near_opt:
    for cf in range(len(behaviors_config)):
        column = get_name_behavior_config(behaviors_config[cf])
        dict_comp_all["follower_obj_near_mean_"+column] = []
for cf in range(len(behaviors_config)):
    column = get_name_behavior_config(behaviors_config[cf])
    dict_comp_all["leader_obj_std_"+column] = []
    dict_comp_all["leader_obj_perc_std_"+column] = []
    dict_comp_all["leader_obj_ess_"+column] = []
    dict_comp_all["leader_obj_rhat_"+column] = []
    if near_opt:
        dict_comp_all["follower_obj_near_std_"+column] = []
        dict_comp_all["follower_obj_near_perc_std_"+column] = []

def add_dict_comp_all(bv_cfg,instance,seed,align,df):
    dict_comp_all["instance"].append(instance)
    dict_comp_all["seed"].append(seed)
    dict_comp_all["align"].append(align)
    dict_comp_all["config"].append(get_name_behavior_config(bv_cfg))

    # Did not find feasible solution.
    if df.shape[0] < len(behaviors_config):
        for cf in range(len(behaviors_config)):
            column = get_name_behavior_config(behaviors_config[cf])
            dict_comp_all["leader_obj_mean_"+column].append("-")
            dict_comp_all["follower_obj_"+column].append("-")
            dict_comp_all["leader_obj_std_"+column].append("-")
            dict_comp_all["leader_obj_perc_std_"+column].append("-")
            dict_comp_all["leader_obj_ess_"+column].append("-")
            dict_comp_all["leader_obj_rhat_"+column].append("-")
            if near_opt:
                dict_comp_all["follower_obj_near_mean_"+column].append("-")
                dict_comp_all["follower_obj_near_std_"+column].append("-")
                dict_comp_all["follower_obj_near_perc_std_"+column].append("-")
        return
        
    for idx, row in df.iterrows():
        bv = int(row["FOLLOWER_BEHAVIOR"])
        tp = 0 if row["COOPERATION_TYPE"] == "-" else int(row["COOPERATION_TYPE"])
        pr = [float(p) for p in row["PARAMS"].split("/")]
        if len(pr) == 1:
            pr = pr[0]

        column = get_name_behavior_config((bv,tp,pr))

        if ("leader_obj_mean_"+column) in dict_comp_all:
            dict_comp_all["leader_obj_mean_"+column].append(float(row["LEADER_OBJ_MEAN"]))
            if bv == 2:
                dict_comp_all["leader_obj_std_"+column].append(math.sqrt(float(row["LEADER_OBJ_VARIANCE"])))
                dict_comp_all["leader_obj_perc_std_"+column].append(math.sqrt(float(row["LEADER_OBJ_VARIANCE"]))/math.fabs(float(row["LEADER_OBJ_MEAN"])))
                dict_comp_all["leader_obj_ess_"+column].append(float(row["LEADER_VARS_ESS"]))
                dict_comp_all["leader_obj_rhat_"+column].append(float(row["LEADER_VARS_RHAT"]))
            else:
                dict_comp_all["leader_obj_std_"+column].append("-")
                dict_comp_all["leader_obj_perc_std_"+column].append("-")
                dict_comp_all["leader_obj_ess_"+column].append("-")
                dict_comp_all["leader_obj_rhat_"+column].append("-")
            dict_comp_all["follower_obj_"+column].append(float(row["FOLLOWER_OBJ_OPT"]))

            if near_opt:
                dict_comp_all["follower_obj_near_mean_"+column].append(float(row["FOLLOWER_OBJ_NEAROPT_MEAN"]))
                if bv == 2:
                    dict_comp_all["follower_obj_near_std_"+column].append(math.sqrt(float(row["FOLLOWER_OBJ_NEAROPT_VARIANCE"])))
                    dict_comp_all["follower_obj_near_perc_std_"+column].append(math.sqrt(float(row["FOLLOWER_OBJ_NEAROPT_VARIANCE"]))/max(0.001,math.fabs(float(row["FOLLOWER_OBJ_NEAROPT_MEAN"]))))
                else:
                    dict_comp_all["follower_obj_near_std_"+column].append("-")
                    dict_comp_all["follower_obj_near_perc_std_"+column].append("-")

dict_comp_mean = {"instance":[], "config":[]}
for cf in range(len(behaviors_config)):
    bv = behaviors_config[cf][0]
    column = get_name_behavior_config(behaviors_config[cf])

    if bv == 2:
        for opt in ["leader_obj_ess_","leader_obj_rhat_","leader_obj_std_","leader_obj_perc_std_"]:
            dict_comp_mean[opt+column] = []
            dict_comp_mean["int"+opt+column] = []
            dict_comp_mean["min"+opt+column] = []
            dict_comp_mean["max"+opt+column] = []
        
        if near_opt:
            for opt in ["follower_obj_near_std_","follower_obj_near_perc_std_"]:
                dict_comp_mean[opt+column] = []
                dict_comp_mean["int"+opt+column] = []
                dict_comp_mean["min"+opt+column] = []
                dict_comp_mean["max"+opt+column] = []
        
def add_dict_comp_mean(bv_cfg,instance):
    dict_comp_mean["instance"].append(instance)
    dict_comp_mean["config"].append(get_name_behavior_config(bv_cfg))

    nb_cfgs = len(behaviors_config)

    for cf in range(nb_cfgs):
        bv = behaviors_config[cf][0]
        column = get_name_behavior_config(behaviors_config[cf])

        if bv == 2:
            for opt in ["leader_obj_ess_","leader_obj_rhat_","leader_obj_std_","leader_obj_perc_std_"]:

                values = dict_comp_all[opt+column][-len(seeds):]
                values = [float(v) for v in values if v != "-"]

                if len(values) == 0:
                    dict_comp_mean[opt+column].append("-")
                    dict_comp_mean["int"+opt+column].append("-")
                    dict_comp_mean["min"+opt+column].append("-")
                    dict_comp_mean["max"+opt+column].append("-")
                elif len(values) == 1:
                    dict_comp_mean[opt+column].append(np.mean(values))
                    dict_comp_mean["int"+opt+column].append(str(round(dict_comp_mean[opt+column][-1],2)) + "+/- 0.00" )
                    dict_comp_mean["min"+opt+column].append(np.min(values))
                    dict_comp_mean["max"+opt+column].append(np.max(values))
                else:
                    dict_comp_mean[opt+column].append(np.mean(values))
                    dict_comp_mean["int"+opt+column].append(str(round(dict_comp_mean[opt+column][-1],2)) + "+/-" + str(round((gettstudent()*np.std(values,ddof=1))/math.sqrt(len(seeds)),2)))
                    dict_comp_mean["min"+opt+column].append(np.min(values))
                    dict_comp_mean["max"+opt+column].append(np.max(values))

            if near_opt:
                for opt in ["follower_obj_near_std_","follower_obj_near_perc_std_"]:

                    values = dict_comp_all[opt+column][-len(seeds):]
                    values = [float(v) for v in values if v != "-"]

                    if len(values) == 0:
                        dict_comp_mean[opt+column].append("-")
                        dict_comp_mean["int"+opt+column].append("-")
                        dict_comp_mean["min"+opt+column].append("-")
                        dict_comp_mean["max"+opt+column].append("-")
                    elif len(values) == 1:
                        dict_comp_mean[opt+column].append(np.mean(values))
                        dict_comp_mean["int"+opt+column].append(str(round(dict_comp_mean[opt+column][-1],2)) + "+/- 0.00")
                        dict_comp_mean["min"+opt+column].append(np.min(values))
                        dict_comp_mean["max"+opt+column].append(np.max(values))
                    else:
                        dict_comp_mean[opt+column].append(np.mean(values))
                        dict_comp_mean["int"+opt+column].append(str(round(dict_comp_mean[opt+column][-1],2)) + "+/-" + str(round((gettstudent()*np.std(values,ddof=1))/math.sqrt(len(seeds)),2)))
                        dict_comp_mean["min"+opt+column].append(np.min(values))
                        dict_comp_mean["max"+opt+column].append(np.max(values))

dict_comp_values_all = {"instance":[],"seed":[],"align":[],"config":[]}
dict_comp_values_mean = {"instance":[],"config":[]}
for cf in range(len(behaviors_config)):
    column = get_name_behavior_config(behaviors_config[cf])
    dict_comp_values_all["gap_leader_"+column] = []
    dict_comp_values_mean["gap_leader_"+column] = []
    dict_comp_values_mean["int"+"gap_leader_"+column] = []
for cf in range(len(behaviors_config)):
    column = get_name_behavior_config(behaviors_config[cf])
    dict_comp_values_all["gap_follower_"+column] = []
    dict_comp_values_mean["gap_follower_"+column] = []
    dict_comp_values_mean["int"+"gap_follower_"+column] = []
if near_opt:
    for cf in range(len(behaviors_config)):
        column = get_name_behavior_config(behaviors_config[cf])
        dict_comp_values_all["gap_follower_near_"+column] = []
        dict_comp_values_mean["gap_follower_near_"+column] = []
        dict_comp_values_mean["int"+"gap_follower_near_"+column] = []

max_gap_comp_values = 5.0  
nb_cfgs_bv0 = len([(tp,pr) for tp in spec_types[0] for pr in behavior_params[0][tp]])
nb_cfgs_bv1 = len([(tp,pr) for tp in spec_types[1] for pr in behavior_params[1][tp]])
nb_cfgs_bv2 = len([(tp,pr) for tp in spec_types[2] for pr in behavior_params[2][tp]])     
def add_dict_comp_values_all_mean(instance):
    nb_cfgs = len(behaviors_config)
    for cf in range(nb_cfgs):
        dict_comp_values_mean["instance"].append(instance)
        dict_comp_values_mean["config"].append(get_name_behavior_config(behaviors_config[cf]))

        nb_seeds = len(seeds)
        for sd in range(nb_seeds): 
            dict_comp_values_all["instance"].append(instance)
            dict_comp_values_all["config"].append(get_name_behavior_config(behaviors_config[cf]))
            dict_comp_values_all["seed"].append(seeds[sd])

            dict_comp_values_all["align"].append(dict_comp_all["align"][-((nb_cfgs-cf)*nb_seeds)+sd])
            
            for column_cf in range(nb_cfgs):
                column = get_name_behavior_config(behaviors_config[column_cf])
                
                if dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-cf)*nb_seeds)+sd] != "-" and dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd] != "-":
                    # dict_comp_values_all["gap_leader_"+column].append(100*(dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-cf)*nb_seeds)+sd] - dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd])/math.fabs(dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd]))
                
                    if behaviors_config[column_cf][0] == 0:
                        cf0 = column_cf
                        if dict_sol_all[0]["gap_dep"][-((nb_cfgs_bv0-cf0)*nb_seeds)+sd] <= max_gap_comp_values:
                            dict_comp_values_all["gap_leader_"+column].append(min(100.0,100*(dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-cf)*nb_seeds)+sd] - dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd])/math.fabs(dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd])))
                        else:
                            dict_comp_values_all["gap_leader_"+column].append("-")
                
                    if behaviors_config[column_cf][0] == 1:
                        cf1 = column_cf - nb_cfgs_bv0
                        if dict_sol_all[1]["gap_1_saa"][-((nb_cfgs_bv1-cf1)*nb_seeds)+sd] <= max_gap_comp_values:
                            dict_comp_values_all["gap_leader_"+column].append(min(100.0,100*(dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-cf)*nb_seeds)+sd] - dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd])/math.fabs(dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd])))
                        else:
                           dict_comp_values_all["gap_leader_"+column].append("-")
                
                    if behaviors_config[column_cf][0] == 2:
                        cf2 = column_cf - nb_cfgs_bv0 - nb_cfgs_bv1
                        if dict_sol_all[2]["gap_1_saa"][-((nb_cfgs_bv2-cf2)*nb_seeds)+sd] <= max_gap_comp_values:
                            dict_comp_values_all["gap_leader_"+column].append(min(100.0,100*(dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-cf)*nb_seeds)+sd] - dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd])/math.fabs(dict_comp_all["leader_obj_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd])))
                        else:
                           dict_comp_values_all["gap_leader_"+column].append("-")
                else:
                    dict_comp_values_all["gap_leader_"+column].append("-")

            for column_cf in range(nb_cfgs):
                column = get_name_behavior_config(behaviors_config[column_cf])

                if dict_comp_all["follower_obj_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd] != "-" and dict_comp_all["follower_obj_"+column][-((nb_cfgs-cf)*nb_seeds)+sd] != "-":
                    if (math.fabs(dict_comp_all["follower_obj_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd]) >= -1e-12) and (math.fabs(dict_comp_all["follower_obj_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd]) <= 1e-12):
                        dict_comp_values_all["gap_follower_"+column].append(np.clip((dict_comp_all["follower_obj_"+column][-((nb_cfgs-cf)*nb_seeds)+sd] - dict_comp_all["follower_obj_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd]),-100.0,100.0))
                    else:
                        dict_comp_values_all["gap_follower_"+column].append(np.clip(100*(dict_comp_all["follower_obj_"+column][-((nb_cfgs-cf)*nb_seeds)+sd] - dict_comp_all["follower_obj_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd])/math.fabs(dict_comp_all["follower_obj_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd]),-100.0,100.0))
                else:
                    dict_comp_values_all["gap_follower_"+column].append("-")

            if near_opt:
                for column_cf in range(nb_cfgs):
                    column = get_name_behavior_config(behaviors_config[column_cf])

                    if dict_comp_all["follower_obj_near_mean_"+column][-((nb_cfgs-cf)*nb_seeds)+sd] != "-" and dict_comp_all["follower_obj_near_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd] != "-":
                        if (math.fabs(dict_comp_all["follower_obj_near_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd]) >= -1e-12) and (math.fabs(dict_comp_all["follower_obj_near_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd]) <= 1e-12):
                            dict_comp_values_all["gap_follower_near_"+column].append(np.clip((dict_comp_all["follower_obj_near_mean_"+column][-((nb_cfgs-cf)*nb_seeds)+sd] - dict_comp_all["follower_obj_near_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd]),-100.0,100.0))
                        else:
                            dict_comp_values_all["gap_follower_near_"+column].append(np.clip(100*(dict_comp_all["follower_obj_near_mean_"+column][-((nb_cfgs-cf)*nb_seeds)+sd] - dict_comp_all["follower_obj_near_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd])/(math.fabs(dict_comp_all["follower_obj_near_mean_"+column][-((nb_cfgs-column_cf)*nb_seeds)+sd])),-100.0,100.0))
                    else:
                        dict_comp_values_all["gap_follower_near_"+column].append("-")

        for column_cf in range(nb_cfgs):
            column = get_name_behavior_config(behaviors_config[column_cf])

            values = dict_comp_values_all["gap_leader_"+column][-len(seeds):]
            values = [float(v) for v in values if v != "-"]

            if len(values) == 0:
                dict_comp_values_mean["gap_leader_"+column].append("-")
                dict_comp_values_mean["int"+"gap_leader_"+column].append("-")
            elif len(values) == 1:
                dict_comp_values_mean["gap_leader_"+column].append(np.mean(values))
                dict_comp_values_mean["int"+"gap_leader_"+column].append(str(round(dict_comp_values_mean["gap_leader_"+column][-1],2)) + "+/- 0.00")
            else:
                dict_comp_values_mean["gap_leader_"+column].append(np.mean(values))
                dict_comp_values_mean["int"+"gap_leader_"+column].append(str(round(dict_comp_values_mean["gap_leader_"+column][-1],2)) + "+/-" + str(round((gettstudent()*np.std(values,ddof=1))/math.sqrt(len(seeds)),2)))
        
        for column_cf in range(nb_cfgs):
            column = get_name_behavior_config(behaviors_config[column_cf])

            values = dict_comp_values_all["gap_follower_"+column][-len(seeds):]
            values = [float(v) for v in values if v != "-"]

            dict_comp_values_mean["gap_follower_"+column].append(np.mean(values))
            dict_comp_values_mean["int"+"gap_follower_"+column].append(str(round(dict_comp_values_mean["gap_follower_"+column][-1],2)) + "+/-" + str(round((gettstudent()*np.std(values,ddof=1))/math.sqrt(len(seeds)),2)))

        if near_opt:
            for column_cf in range(nb_cfgs):
                column = get_name_behavior_config(behaviors_config[column_cf])

                values = dict_comp_values_all["gap_follower_near_"+column][-len(seeds):]
                values = [float(v) for v in values if v != "-"]

                dict_comp_values_mean["gap_follower_near_"+column].append(np.mean(values))
                dict_comp_values_mean["int"+"gap_follower_near_"+column].append(str(round(dict_comp_values_mean["gap_follower_near_"+column][-1],2)) + "+/-" + str(round((gettstudent()*np.std(values,ddof=1))/math.sqrt(len(seeds)),2)))

############################################################################################################
#                                                 Main
############################################################################################################

def main():
    path = name_folder + 'compiled/results_nearopt' + str(near_opt) 
    if near_opt:
        path += '_eps' + str(eps_near_opt) + "_mult" + str(mul_near_opt)
    path += '.xlsx'
    with pd.ExcelWriter(path) as writer:
        # General instances.
        for instance in instances:
            # General definition of behaviors.
            for bv in behaviors:
                for tp in spec_types[bv]:
                    for pr in behavior_params[bv][tp]:
                        # Seeds for instances.
                        for seed in seeds:
                            # Read solution and comparison files.
                            init_dict_sol_all(bv,tp,pr,instance,seed)

                            if bv == 0 or bv == 2:
                                name_sol_file = get_sol_file_name(instance,seed,(bv,tp,pr))
                                name_cmp_file = get_cmp_file_name(instance,seed,(bv,tp,pr))

                                print(name_sol_file)
                                print(name_cmp_file)

                                try:
                                    df_sol = pd.read_csv(name_sol_file,delimiter=';')
                                except pd.errors.EmptyDataError:
                                    df_sol = pd.DataFrame(columns=["empty"])
                                
                                try:
                                    df_cmp = pd.read_csv(name_cmp_file,delimiter=';')
                                except pd.errors.EmptyDataError:
                                    df_cmp = pd.DataFrame(columns=["empty"])

                                align = add_dict_sol_all_align(bv,df_sol) 
                                if bv == 0:
                                    add_dict_sol_all_dep(bv,df_sol)
                                if bv == 2:
                                    add_dict_sol_all_saa(bv,df_sol)
                                
                                add_dict_comp_all((bv,tp,pr),instance,seed,align,df_cmp)

                            if bv == 1:
                                name_sol_file_dep = get_sol_file_name(instance,seed,(bv,tp,pr),"dep")
                                name_sol_file_saa = get_sol_file_name(instance,seed,(bv,tp,pr),"saa")
                                name_cmp_file = get_cmp_file_name(instance,seed,(bv,tp,pr),"saa")

                                print(name_sol_file_dep)
                                print(name_sol_file_saa)
                                print(name_cmp_file)

                                df_sol_saa = pd.read_csv(name_sol_file_saa,delimiter=';')
                                df_sol_dep = pd.DataFrame(columns=["empty"])
                                df_sol_dep = pd.read_csv(name_sol_file_dep,delimiter=';')

                                try:
                                    df_cmp = pd.read_csv(name_cmp_file,delimiter=';')
                                except pd.errors.EmptyDataError:
                                    df_cmp = pd.DataFrame(columns=["empty"])
                                
                                align = add_dict_sol_all_align(bv,df_sol_dep)
                                add_dict_sol_all_dep(bv,df_sol_dep)
                                add_dict_sol_all_saa(bv,df_sol_saa)
                                add_dict_sol_all_comp(bv)

                                add_dict_comp_all((bv,tp,pr),instance,seed,align,df_cmp)

                        # Mean for solution files. 
                        add_dict_sol_mean(bv,tp,pr,instance)

                        # Mean for comparison files.
                        add_dict_comp_mean((bv,tp,pr),instance)

            # Definition of comparison values.
            add_dict_sol_values_all(instance)
            add_dict_sol_values_mean(instance)
            add_dict_comp_values_all_mean(instance)

        # Generate folders excel file.
        for bv in behaviors:
            excel_sol_all = pd.DataFrame(dict_sol_all[bv])
            excel_sol_all.to_excel(writer,sheet_name="sol_all_"+str(bv))

            excel_sol_mean = pd.DataFrame(dict_sol_mean[bv])
            excel_sol_mean.to_excel(writer,sheet_name="sol_mean_"+str(bv))

            if bv == 1:
                excel_sol_mean_per_type = pd.DataFrame(get_dict_sol_mean_per_type())
                excel_sol_mean_per_type.to_excel(writer,sheet_name="sol_mean_per_type_"+str(bv)) 

            excel_sol_mean_final = pd.DataFrame(get_dict_sol_mean_final(bv))  
            excel_sol_mean_final.to_excel(writer,sheet_name="sol_mean_final_"+str(bv)) 

        excel_sol_values_all = pd.DataFrame(dict_sol_values_all)
        excel_sol_values_all.to_excel(writer,sheet_name="sol_values_all")

        excel_sol_values_mean = pd.DataFrame(dict_sol_values_mean)
        excel_sol_values_mean.to_excel(writer,sheet_name="sol_values_mean") 

        excel_comp_all = pd.DataFrame(dict_comp_all)
        excel_comp_all.to_excel(writer,sheet_name="comp_all") 

        excel_comp_mean = pd.DataFrame(dict_comp_mean)
        excel_comp_mean.to_excel(writer,sheet_name="comp_mean") 

        excel_comp_values_all = pd.DataFrame(dict_comp_values_all)
        excel_comp_values_all.to_excel(writer,sheet_name="comp_values_all")

        excel_comp_values_mean = pd.DataFrame(dict_comp_values_mean)
        excel_comp_values_mean.to_excel(writer,sheet_name="comp_values_mean")

if __name__ == "__main__":
    main()