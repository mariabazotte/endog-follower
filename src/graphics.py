
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from matplotlib.ticker import MaxNLocator

############################################################################################################
#                                              Behaviors
############################################################################################################

# graphic between gap de and gap saa problem 1
# graphic time and gaps per instance for different behaviors: 3 lines for behavior general types
# and multiple lines for subypes inside general type.
# graphics per alignment value and compute number of instances per alignment value

# return with cases with val follower obj = 0
# look at removing extremes for average.

# remove boxplot 
# max gap 100%
# remove cases not solved to optimality to avoid negative values at comparison gap behavior

# Follower near/optimality
near_opt     = 1
eps_near_opt = 0.01
mul_near_opt = 0

# Alignment between leader and follower objectives information.
align_types = 4
align_ranges = [[-1.0,-0.15],[-0.15,0.0],[0.0,0.15],[0.15,1.0]]

# Follower behaviors
behaviors  = [0, 1, 2] 
spec_types = {0: [0], 1: [0,1,2,3,4,5], 2: [0,1,2,3,4]} 

name_behaviors  = ["fixstrongweak","depstrongweak","depgeneral"]
name_spec_types = { 0: [""], 
                    1: ["proportional","threshold","strong","fragile","strong_power","fragile_power"], 
                    2: ["neutral","proportional","fragile","fragile_power","strong_power"]}

# Parameters for each behavior type.
fix_strwk_coop  = [1.0, 0.7, 0.5, 0.3, 0.0]
dep_strwk_thr = [0.5, 2.0, 5.0]
dep_strwk_str = [0.5, 2.0, 5.0, 10.0]
dep_strwk_frg = [0.5, 2.0, 5.0, 10.0]
dep_strwk_str_pow = [2.0, 5.0, 10.0]
dep_strwk_frg_pow = [2.0, 5.0, 10.0]
dep_gen_nbint  = [6.0]
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

name_dep_strwk = ['Prp','Thr','Str','Frg','Str-p','Frg-p']
name_dep_gen   = ['Ntr','Prp','Frg','Frg-p','Str-p']
def get_name_behavior_config_graphic(bv_config):
    bv = bv_config[0]
    tp = bv_config[1]
    pr = bv_config[2]

    if bv == 0:
        return rf'$\beta^{{{pr}}}$'
    if bv == 1:
        if tp == 0:
            return rf'$\beta_{{\mathrm{{{name_dep_strwk[tp]}}}}}(\cdot)$'
        if tp == 1:
            return rf'$\beta_{{\mathrm{{{name_dep_strwk[tp]}}}}}^{{{pr}}}(\cdot)$' #, $\delta = {pr}$'
        if tp == 2 or tp == 3: 
            return rf'$\beta_{{\mathrm{{{name_dep_strwk[tp]}}}}}^{{{pr}}}(\cdot)$' #, $a = {pr}$'
        if tp == 4 or tp == 5: 
            return rf'$\beta_{{\mathrm{{{name_dep_strwk[tp]}}}}}^{{{pr}}}(\cdot)$' #, $p = {pr}$'
    if bv == 2:
        if tp == 0:
            return rf'$\pi_{{\mathrm{{{name_dep_gen[tp]}}}}}(\cdot)$'
        if tp == 1 or tp == 2 or tp == 3 or tp == 4:
            return rf'$\pi_{{\mathrm{{{name_dep_gen[tp]}}}}}^{{{int(pr[0])},{pr[1]}}}(\cdot)$' #, $K = {pr[0]}$, $\gamma = {pr[1]}$'

label_behaviors = [get_name_behavior_config_graphic(bv_cfg) for bv_cfg in behaviors_config]

def get_name_behavior_config_general_graphic(bv,tp):
    if bv == 0:
        return rf'$\beta^{{{tp}}}$'
    if bv == 1:
        return rf'$\beta_{{\mathrm{{{name_dep_strwk[tp]}}}}}(\cdot)$'
    if bv == 2:
        return rf'$\pi_{{\mathrm{{{name_dep_gen[tp]}}}}}(\cdot)$'

############################################################################################################
#                                              Directory
############################################################################################################

# Define file and save paths.
file_path = '../results-before/compiled/results_nearopt' + str(near_opt) 
if near_opt:
    file_path += '_eps' + str(eps_near_opt) + "_mult" + str(mul_near_opt)
file_path += '.xlsx'

save_path = '../results-before/compiled/graphics/'
sol_values_path = save_path + 'sol_values/'
cmp_values_path = save_path + 'cmp_values/'
perf_path = save_path + 'perf/'
perf_path_bv = save_path + 'perf_bv/'

############################################################################################################
#                                              Graphics
############################################################################################################

colors = {'leader_obj_': '#0072B2','follower_obj_': '#009E73','follower_obj_near_mean_': '#E69F00'}
styles = {'leader_obj_':'-','follower_obj_':'--','follower_obj_near_mean_':'-.'}
hatch_styles = {'leader_obj_':'/','follower_obj_':'x','follower_obj_near_mean_':'o'}

def plot_sol_metric(y_column, y_column_name, align):
    """
    Plots the average and quantiles of a metric (y_column) by follower behavior, and saves the figure.
    
    Parameters:
    - y_column: Name of the column on the excel file to use for the y-axis. ('leader_obj_','follower_obj_','follower_obj_near_mean_')
    - y_column_name: Name to display on the y-axis of the figure.
    - align: [-1,0,1,2,3] -1 for all instances 
    """

    if y_column not in ['leader_obj_','follower_obj_','follower_obj_near_mean_']:
        raise ValueError(f"Invalid y_column: {y_column}")

    # Columns correspnding to metric y_column
    columns = [y_column+get_name_behavior_config(bv_cfg) for bv_cfg in behaviors_config]

    # Reading sheet.
    df = pd.read_excel(file_path, sheet_name="sol_values_all")

    if align >= 0 and align <= 3:
        df = df[(df["align"] >= align_ranges[align][0]) & (df["align"] <= align_ranges[align][1] - 1e-12)]

    df = df[columns]

    df = df[~(df == "-").any(axis=1)]
    df = df.apply(pd.to_numeric)

    df = df.clip(lower=-100,upper=100)

    # print(df.mean())       # line plot values
    # print(df.describe())   # boxplot min, 25%, 50%, 75%, max

    ###############################################################

    y = df.mean()
    # y = df.median()
    y_std = df.std()
    x = range(0,10*len(y),10)

    n = df.count()

    y_se = y_std / np.sqrt(n)
    ci = 1.96 * y_se  # 95% confidence interval

    plt.figure(figsize=(40, 15))
    plt.plot(x, y, color=colors[y_column], linestyle=styles[y_column], marker='o', linewidth=5)
    plt.fill_between(x, y - ci, y + ci, color=colors[y_column], alpha=0.25)

    # plt.errorbar(x,y,yerr=y_std,fmt='o',linestyle=styles[y_column],color=colors[y_column],linewidth=5,capsize=8)

    plt.xlabel(r'Follower behavior')
    plt.ylabel(f'{y_column_name}')

    # plt.ylim(bottom=-0.5)
    plt.xticks(x, label_behaviors, rotation=45)

    # plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(sol_values_path + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + "_align" + str(align) + "_line_" + f"{y_column}.png", dpi=300)
    plt.close()

    ###############################################################

    plt.figure(figsize=(40, 15))

    sns.boxplot(data=df, showfliers=True)

    plt.xlabel(r'Follower behavior')
    plt.ylabel(f'{y_column_name}')
    
    # plt.ylim(bottom=-0.5)
    plt.xticks(ticks=range(len(columns)), labels=label_behaviors, rotation=45)

    plt.grid(True, axis="y")
    plt.tight_layout()

    # plt.savefig(sol_values_path + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + "_align" + str(align) + "_boxplot_" + f"{y_column}.png", dpi=300)
    plt.close()


def plot_cmp_metric(bv_cfg, y_column, y_column_name, align):
    """
    Plots the average and quantiles of a metric (y_column) by follower behavior, and saves the figure.
    
    Parameters:
    - y_column: Name of the column to use for the y-axis. ('gap_leader_','gap_follower_','gap_follower_near_')
    - y_column_name: Name to display on the y-axis 
    - align: [-1,0,1,2,3] -1 for all instances
    """

    # Columns corresponding to metric y_column
    columns = [y_column+get_name_behavior_config(column_bv_cfg) for column_bv_cfg in behaviors_config]
    rows = get_name_behavior_config(bv_cfg)

    # Reading sheet.
    df = pd.read_excel(file_path, sheet_name="comp_values_all")

    if align >= 0 and align <= 3:
        df = df[(df["align"] >= align_ranges[align][0]) & (df["align"] <= align_ranges[align][1] - 1e-12)]

    df = df[df["config"] == rows]
    df = df[columns]

    df = df[~(df == "-").any(axis=1)]
    df = df.apply(pd.to_numeric)

    df = df.clip(lower=0)

    ###############################################################

    y = df.mean()
    # y = df.median()
    x = range(len(y))

    plt.figure(figsize=(30, 15))
    plt.plot(x, y, marker='o', linewidth=5)

    plt.xlabel(r'Follower behavior')
    plt.ylabel(f'{y_column_name}')

    # plt.ylim(bottom=0)
    plt.xticks(x, label_behaviors, rotation=45)

    # plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(cmp_values_path + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + "_align" + str(align) + "_line_" + f"{y_column}" + "_" + str(bv_cfg) + ".png", dpi=300)
    plt.close()

    ###############################################################

    plt.figure(figsize=(30, 15))

    sns.boxplot(data=df,showfliers=True)

    plt.xlabel(r'Follower behavior')
    plt.ylabel(f'{y_column_name}')
    
    # plt.ylim(bottom=0)
    plt.xticks(ticks=range(len(columns)), labels=label_behaviors, rotation=45)

    plt.grid(True, axis="y")
    plt.tight_layout()
    # plt.savefig(cmp_values_path + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + "_align" + str(align) + "_boxplot_" + f"{y_column}" + "_" + str(bv_cfg) + ".png", dpi=300)
    plt.close()

solver_columns = ['_saa','_dep']
gap_columns = ['gap_1','gap']
solver_names = [r'$\texttt{TR-SAA}$',r'$\texttt{DE}$']
solver_colors = ['#E41A1C',  '#377EB8']
solver_styles = ['-', '--']
def plot_perf(tp):
    """
    Plots the accumulated opt. gap, time, and ub gap by percentage of instance for the different methods 
    (deterministic equivalent and transformed saa) for the different types {tp} of the decision-dependent 
    strong weak behavior, and saves the figures. 
    
    Parameters:
    - tp in [-1,0,1,2,3,4,5]: -1 i for all types
    """    

    df = pd.read_excel(file_path, sheet_name="sol_all_1")
    if tp != -1:
        df = df[df["type"] == name_spec_types[1][tp]]
    x = [i/df.shape[0]*100 for i in range(1,df.shape[0]+1)]
    
    ################################
    
    # Optimality gap.
    plt.figure(figsize=(16, 10))

    for sol in [0,1]:
        y_column = gap_columns[sol] + solver_columns[sol]
        df[y_column] = df[y_column].clip(upper=100) # max gap 100% (avoid problems as we have positive and negative values)
        y = df[y_column].sort_values()

        plt.plot(x, y, label=solver_names[sol], color=solver_colors[sol], linestyle=solver_styles[sol], linewidth=5)
    
    plt.xlabel(r'Instances (\%)')
    # plt.ylabel(r'${\textnormal{gap}}_{opt}$ (\%)')
    
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(perf_path + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + f"_opt_gap_tp" + str(tp) + ".png", dpi=300)
    plt.close() 
    
    ##############################
        
    # Time
    plt.figure(figsize=(16, 10))
    for sol in [0,1]:
        y_column = 'time' + solver_columns[sol]
        y = df[y_column].sort_values()

        plt.plot(x, y, label=solver_names[sol], color=solver_colors[sol], linestyle=solver_styles[sol], linewidth=5)
        
    plt.xlabel(r'Instances (\%)')
    # plt.ylabel(r'Time (Seconds)')
    
    plt.legend()
    plt.yticks([0, 1800, 3600, 5400, 7200, 9000, 10800])

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(perf_path + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + f"_time_" + str(tp) + ".png", dpi=300)
    plt.close()
    
    ############################ 
    
    # Oportunity gap
    plt.figure(figsize=(16, 10))
    y_column = 'gap_ub'
    df[y_column] = df[y_column].clip(upper=100) # max gap 100% (avoid problems as we have positive and negative values)
    y = df[y_column].sort_values()

    plt.plot(x, y, label = solver_names[sol], color=solver_colors[sol], linestyle=solver_styles[sol], linewidth=5)
            
    plt.xlabel(r'Instances (\%)')
    # plt.ylabel(r'${\textnormal{gap}}_{ub^*}(\%)$')

    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    
    plt.savefig(perf_path + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + "_ub_gap_" + str(tp) + ".png", dpi=300)
    plt.close() 


behavior_names = [r'$\beta$',r'$\beta(\cdot)$',r'$\pi(\cdot)$']
behavior_styles = ['-', '--', '-.']
behavior_colors = ["#882255", "#44AA99", "#DDCC77"]
def plot_perf_all_bv():
    """
    Plots the accumulated (mean per instance and behavior general type) opt. gap and time by percentage of instance for the three
    different general behavior types (fixed strong weak, decision-dependent strong weak,
    generalized decision-dependent), and saves the figures. 

    """ 

    df = list()
    df.append(pd.read_excel(file_path, sheet_name="sol_all_0"))
    df.append(pd.read_excel(file_path, sheet_name="sol_all_1"))
    df.append(pd.read_excel(file_path, sheet_name="sol_all_2"))

    x = [i/df[0][["instance", "seed"]].drop_duplicates().shape[0]*100 for i in range(1,df[0][["instance", "seed"]].drop_duplicates().shape[0]+1)]

    ################################
    
    # Optimality gap.
    plt.figure(figsize=(16, 10))

    for bv in behaviors:
        y = None
        if bv == 0:
            y = (df[bv].groupby(["instance", "seed"], as_index=False)["gap_dep"].mean())["gap_dep"].sort_values()
        if bv == 1 or bv == 2:
            y = (df[bv].groupby(["instance", "seed"], as_index=False)["gap_1_saa"].mean())["gap_1_saa"].sort_values()

        plt.plot(x, y, label=behavior_names[bv], color=behavior_colors[bv], linestyle=behavior_styles[bv], linewidth=5)
    
    plt.xlabel(r'Instances (\%)')
    # plt.ylabel(r'${\textnormal{gap}}_{opt}$ (\%)')
    
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.savefig(perf_path_bv + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + f"_opt_gap_bv_all.png", dpi=300)
    plt.close() 

    ##############################
        
    # Time
    plt.figure(figsize=(16, 10))
    for bv in behaviors:
        y = None
        if bv == 0:
            y = (df[bv].groupby(["instance", "seed"], as_index=False)["time_dep"].mean())["time_dep"].sort_values()
        if bv == 1:
            y = (df[bv].groupby(["instance", "seed"], as_index=False)["time_saa"].mean())["time_saa"].sort_values()
        if bv == 2:
            y = (df[bv].groupby(["instance", "seed"], as_index=False)["time_lb_saa"].mean())["time_lb_saa"].sort_values()
    
        plt.plot(x, y, label=behavior_names[bv], color=behavior_colors[bv], linestyle=behavior_styles[bv], linewidth=5)
        
    plt.xlabel(r'Instances (\%)')
    # plt.ylabel(r'Time (Seconds)')

    plt.legend()
    plt.yticks([0, 1800, 3600, 5400, 7200, 9000, 10800])

    plt.grid(True)
    plt.tight_layout()

    plt.savefig(perf_path_bv + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + f"_time_bv_all.png", dpi=300)
    plt.close()
    
    ############################ 


spec_types_plot = {0: fix_strwk_coop, 1: spec_types[1], 2: spec_types[2]}
behavior_spec_names = {bv: {tp: get_name_behavior_config_general_graphic(bv,tp) for tp in spec_types_plot[bv]} for bv in behaviors}
behavior_spec_colors = {0: {0.0:"#CC79A7", 0.3:"#F0E442", 0.5:"#999999", 0.7:"#56B4E9", 1.0:"#D55E00"},
                        1: {0:"#882255", 1:"#44AA99", 2:"#DDCC77", 3:"#AA4499", 4:"#117733", 5:"#332288"},
                        2: {0:"#661100", 1:"#6699CC", 2:"#AA4466", 3:"#4477AA", 4:"#228833"}}
behavior_spec_styles = {0: {0.0:"-", 0.3:"--", 0.5:"-.", 0.7:":", 1.0:(0, (5, 2))}, 
                        1: {0:"-", 1:"--", 2:"-.", 3:":", 4:(0, (5, 2)), 5:(0, (3, 1, 1, 1))},
                        2: {0:"-", 1:"--", 2:"-.", 3:":", 4:(0, (5, 2))}}

def plot_perf_spec_bv():
    """
    Plots the accumulated (mean per instance and behavior specific subtypes) opt. gap and time by percentage of instance for 
    the different subtypes of each general behavior type, and saves the figures. 

    """ 

    for bv in behaviors:
        df = pd.read_excel(file_path, sheet_name="sol_all_"+ str(bv))

        x = [i/df[["instance", "seed"]].drop_duplicates().shape[0]*100 for i in range(1,df[["instance", "seed"]].drop_duplicates().shape[0]+1)]

        ################################
    
        # Optimality gap.
        plt.figure(figsize=(16, 10))

        for tp in spec_types_plot[bv]:
            y = None
            if bv == 0:
                df_tp = df[df["params"] == tp]
                y = (df_tp.groupby(["instance","seed"], as_index=False)["gap_dep"].mean())["gap_dep"].sort_values()
            if bv == 1 or bv == 2:
                df_tp = df[df["type"] == name_spec_types[bv][tp]]
                y = (df_tp.groupby(["instance", "seed"], as_index=False)["gap_1_saa"].mean())["gap_1_saa"].sort_values()

            plt.plot(x, y, label=behavior_spec_names[bv][tp], color=behavior_spec_colors[bv][tp], linestyle=behavior_spec_styles[bv][tp], linewidth=5)
        
        plt.xlabel(r'Instances (\%)')
        # plt.ylabel(r'${\textnormal{gap}}_{opt}$ (\%)')
        
        plt.legend()
        plt.grid(True)
        plt.tight_layout()

        plt.savefig(perf_path_bv + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + f"_opt_gap_bv" + str(bv) + ".png", dpi=300)
        plt.close() 

        ##############################
        
        # Time
        plt.figure(figsize=(16, 10))
        for tp in spec_types_plot[bv]:
            y = None
            if bv == 0:
                df_tp = df[df["params"] == tp]
                y = (df_tp.groupby(["instance", "seed"], as_index=False)["time_dep"].mean())["time_dep"].sort_values()
            if bv == 1:
                df_tp = df[df["type"] == name_spec_types[bv][tp]]
                y = (df_tp.groupby(["instance", "seed"], as_index=False)["time_saa"].mean())["time_saa"].sort_values()
            if bv == 2:
                df_tp = df[df["type"] == name_spec_types[bv][tp]]
                y = (df_tp.groupby(["instance", "seed"], as_index=False)["time_lb_saa"].mean())["time_lb_saa"].sort_values()
        
            plt.plot(x, y, label=behavior_spec_names[bv][tp], color=behavior_spec_colors[bv][tp], linestyle=behavior_spec_styles[bv][tp], linewidth=5)
            
        plt.xlabel(r'Instances (\%)')
        # plt.ylabel(r'Time (Seconds)')

        plt.legend()
        plt.yticks([0, 1800, 3600, 5400, 7200, 9000, 10800])

        plt.grid(True)
        plt.tight_layout()

        plt.savefig(perf_path_bv + "nearopt" + str(near_opt) + "_eps" + str(eps_near_opt) + f"_time_bv" + str(bv) +".png", dpi=300)
        plt.close()

############################################################################################################
#                                              Main
############################################################################################################

def main():
    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    os.makedirs(os.path.dirname(sol_values_path), exist_ok=True)
    os.makedirs(os.path.dirname(cmp_values_path), exist_ok=True)
    os.makedirs(os.path.dirname(perf_path), exist_ok=True)
    os.makedirs(os.path.dirname(perf_path_bv), exist_ok=True)
    
    # Set global font sizes
    plt.rcParams.update({
        'font.family': 'serif',
        # 'legend.title_fontsize': 44,
        'text.usetex': True,
        'axes.titlesize': 48,
        'axes.labelsize': 50,
        'xtick.labelsize': 35,
        'ytick.labelsize': 35,
        'legend.fontsize': 44,
        'font.size': 44,
        'axes.linewidth': 3,
        'lines.markersize': 10
    })

    for align in range(-1,align_types):
        plot_sol_metric('leader_obj_',r'$\overline{\texttt{gap}}_{F}^{b \mid opt}$ (\%)',align)
        plot_sol_metric('follower_obj_',r'$\overline{\texttt{gap}}_{f}^{b \mid opt}$ (\%)',align)
        plot_sol_metric('follower_obj_near_mean_',r'$\overline{\texttt{gap}}_{f_\epsilon}^{b \mid opt}$ (\%)',align)

    # for bv in behaviors_config:
    for pr in behavior_params[0][0]:
        bv = (0,0,pr)
        for align in range(-1,align_types):
            plot_cmp_metric(bv,'gap_leader_',r'$\overline{\texttt{gap}}_{F}^{b \leftarrow b}$ (\%)',align)
            # plot_cmp_metric(bv,'gap_follower_',r'$\overline{\texttt{gap}}_{f}^{b \leftarrow b}$ (\%)',align)
            # plot_cmp_metric(bv,'gap_follower_near_',r'$\overline{\texttt{gap}}_{f_\epsilon}^{b \leftarrow b}$ (\%)',align)

    # Only behavior 1 solved with DE and SAA.
    for tp in spec_types[1]:
        plot_perf(tp)
    plot_perf(-1) # All types together.

    # Plot performance each type of behavior.
    plot_perf_all_bv()
    plot_perf_spec_bv()



if __name__ == "__main__":
    main()