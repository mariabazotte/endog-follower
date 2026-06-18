# endog-follower

## Description

This project contains the code for the paper "Intermediate Bilevel Optimization: Modeling Endogenous Follower Tie-Breaking Behavior" by Maria Bazotte, Margarida Carvalho, and Thibaut Vidal (https://arxiv.org/abs/2606.18475). It contains the implementation of the methods introduced in this paper to solve the intermediate bilevel program (1-BO), which contains endogenous uncertainty.

## Building

### Instances

The folder `data/` contains the data files needed for experiments in the paper. The file `read_instance.py` was used to read the instances downloaded from BOBLib (file general-bilevel.tar.gz) and generate the `.txt` instance files in folder `data/instances/` used in this work. 

### Main Code 

Folder `scr/` contains the source code. To run the source code, execute:

```
make 
./exe -parameters
```

Where we have the following parameters:

#### Mandatory parameters: 

- instancefile ../data/instances/denegre-miblp_20_15_50_0110_10_1.txt    # Path to instance file
- -follower 0   # Type of follower behavior (0: strong-weak fixed intermediate (SWI-BO), 1: strong-weak decision-dependent intermediate (SWDI-BO), 2: generalized decision-dependent intermediate)
- -fix_strwk    # Subtype of strong-weak fixed intermediate (fixed cooperation parameter, e.g., beta = 0, 0.3, 0.5, 0.7, 1) (0 is the pessimistic approach, and 1 is the optimistic approach).
- -dep_strwk 0  # Subtype of strong-weak decision-dependent intermediate (0: proportional, 1: threshold, 2: sturdy, 3: fragile, 4: sturdy-p, 5: fragile-p)
- -dep_gen 0    # Subtype of generalized decision-dependent intermediate (0: neutral, 1: proportional, 2: fragile moderate, 3: fragile pronounced, 4: sturdy)

#### Optinal parameters for instance and behavior: 
- -near_opt 0                 # 0: follower is optimal, 1: follower is near-optimal
- -eps_near_optimal 0.01      # Value of deviation for near-optimality
- -thr_param  2               # Parameter (delta) for threshold strong-weak decision-dependent intermediate
- -str_param  2               # Parameter (a) for for sturdy strong-weak decision-dependent intermediate
- -frag_param  2              # Parameter (a) for for fragile strong-weak decision-dependent intermediate
- -str_pow_param 2            # Parameter (p) for for sturdy-p strong-weak decision-dependent intermediate
- -frag_pow_param 2           # Parameter (p) for for fragile-p strong-weak decision-dependent intermediate
- -scal_param 2               # Scaling parameter (gamma) for generalized decision-dependent intermediate fragile moderate, fragile pronounced, sturdy.
- -method_dep_gen 4           # Method used to solve the generalized decision-dependent intermediate bilevel program. We set it to 4, being the integer L-shaped method presented in the paper.

#### Optinal parameters for solution procedure: 
- -sol 0             # 0: Sample average approximation (SAA) (defined for strong weak and generalized decision-dependent intermediate), 1: Deterministic equivalent program (DEP) (defined for strong weak fixed and decision-dependent intermdiate)
- -timelimit 10800   # Time limit for optimization
- -nbthreads 1       # Number of threads for the Gurobi solver
- -verbose 0         # 0 -> nothing, 1-> moderated, 2-> intensive (log files)

#### Optinal parameters for SAA and MCMC: 
 - -nbproblemsSAA 1           # Number of SAA problems
 - -nbscenariosSAA 100        # Number of scenarios for each SAA problem
 - -nbthinningSAA 1           # Thinning for SAA problem (for markov chain transformation)
 - -nbvalidateproblems 4      # Number of markov chains for validation
 - -nbvalidatescenarios 10000 # Number of scenarios for validation
 - -nbvalidatethinning 15     # Thinning for validation (for MCMC Hit-and-Run)
 - -coord_har 1               # 0: Use standard MCMC Hit-and-Run, 1: Use coordinate MCMC Hit-and-Run


#### Folders main code

The main code contains the following folders:

- -`src/input/` # This folder contains the class Input, which stores information on input parameters provided by the user to the solver during code execution.
- -`src/instance/`  # This folder contains the class Instnace, which stores the data of the current instance being evaluated. It also evaluates leader decisions by returning the leader and follower objective values they define using the corresponding evaluation process, such as MCMC Hit-and-Run. This folder contains the subfolders `instance/bilevelmodel/`, which defines a general class of instances of bilevel programs investigated in this work, and `instance/performance/`, which defined performance metrics for the MCMC Hit-and-Run convergence.
- -`src/leadersolver/` # This folder contains the implementation of the upper-level problem with classes AbstractLeaderSolver, LeaderSolver and SAALeaderSolver.
- `src/followersolver/` # This folder contains the implementation of the follower problem, according to the type of measure: folder `src/followersolver/fixstrongweak` with class FixStrongWeakFollowerSolver for strong-weak fixed cooperation measures, folder `src/followersolver/depstrongweak` with classes SAADepStrWkFollowerSolver (solving SAA program) or DEDepStrWkFollowerSolver (solving DE program) for strong-weak decision-dependent cooperation measures, and folder `src/followersolver/depgeneral` with class DepGeneralFollowerSolver for generalized decision-dependent measures.

## Results

The results are saved in a repository such as: 

- results
    - compiled 
        - graphics
    - depgeneral 
        - fragile
        - fragile_power 
        - neutral
        - proportional
        - strong-power
    - depstrongweak 
        - dep
            - fragile
            - fragile_power
            - proportional
            - strong
            - strong_power
            - threshold
        - saa
            - fragile
            - fragile_power
            - proportional
            - strong
            - strong_power
            - threshold
    - fixstrongweak:


Each test generates a separate `.csv` file in the folder corresponding to the endogenous measure and method (saa or dep). The folder `results/compiled` contains an Excel file that compiles the experimental results for the paper. These file contains the information present in all tables of the paper. Th folder `results/graphics` contains the graphics (corresponding to the experimental results) of the paper.

## Replicating

The folder `scripts/` contains the sbatch files for replicating the results of the paper. Specifically, the files `jobs_dep_strwk_dep.sh`, `jobs_dep_strwk_saa.sh`, `jobs_fix_strwk.sh`, and `jobs_gen_cb.sh` contain the commands for the paper's results. After executing these sbatch files, the resulting `.csv` files are written in the folder `results`, according to the configuration explained before.
        
### Unifying paper results

The file `unify_results.py` unifies the results contained in the `.csv` files, creating the Excel file with the compiled experimental results of the paper. After running the sbatch files, execute:

```
python unify_results.py
```

You can modify this file to see the results of other parameters. The unified results are in the folder `results/compiled`.

### Generating graphics

After unifying results with `unify_results.py` and generating the Excel file, execute: 

```
python graphics.py
```

You can modify this file to see the results of other parameters. The unified results are in the folder `results/compiled/graphics`.

