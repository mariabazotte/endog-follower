#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --time=30:00:00
#SBATCH --array=1-3230
#SBATCH --output=../slurm/depgen/cb_%A_%a.out

i=1
verbose=1
nbthreads=1
timelimit=10800

nbproblems=1
nbscenarios=1500 
nbthinning=10

nbvalidateproblems=4
nbvalidatescenarios=10000
nbvalidatethinning=50

# "denegre-miblp_20_15_50_0110_10_" "denegre-miblp_20_20_50_0110_5_" 
# "denegre-miblp_20_20_50_0110_10_" "denegre-miblp_20_20_50_0110_15_" 
# "denegre-milp_30_30_30_2310_5_" "denegre-milp_30_30_30_2310_10_"
# "denegre-miblp_20_20_50_0110_5_e1_" "denegre-milp_30_30_30_2310_5_e1_" 
# "denegre-milp_30_30_30_2310_5_e2_" "denegre-milp_30_30_30_2310_10_e1_" 
# "denegre-milp_30_30_30_2310_10_e2_"
# "xuwang-bmilplib_10_" "xuwang-bmilplib_60_" 
# "xuwang-bmilplib_110_" "xuwang-bmilplib_160_" 
# "xuwang-bmilplib_210_" "xuwang-bmilplib_260_"

for instance in "denegre-miblp_20_15_50_0110_10_" "denegre-miblp_20_20_50_0110_5_" "denegre-miblp_20_20_50_0110_10_" "denegre-miblp_20_20_50_0110_15_" "denegre-milp_30_30_30_2310_5_" "denegre-milp_30_30_30_2310_10_" "denegre-miblp_20_20_50_0110_5_e1_" "denegre-milp_30_30_30_2310_5_e1_" "denegre-milp_30_30_30_2310_5_e2_" "denegre-milp_30_30_30_2310_10_e1_" "denegre-milp_30_30_30_2310_10_e2_" "xuwang-bmilplib_10_" "xuwang-bmilplib_60_" "xuwang-bmilplib_110_" "xuwang-bmilplib_160_" "xuwang-bmilplib_210_" "xuwang-bmilplib_260_"
do
for epsnearopt in 0.01
do
for seed in 1 2 3 4 5 6 7 8 9 10
do  
    # Decision-dependent general.
    for method in 4
    do 
        # Neutral 
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            ./exe -instancefile "../data/instances/${instance}${seed}_mod.txt" -follower 2 -dep_gen 0 -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbthinningSAA $nbthinning -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning -method_dep_gen $method
        fi
        (( i = $i +1 ))

        for nbint in -1
        do 
            # Proportional
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/${instance}${seed}_mod.txt" -follower 2 -dep_gen 1 -nb_intv $nbint -scal_param 1.0 -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbthinningSAA $nbthinning  -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning -method_dep_gen $method
            fi
            (( i = $i +1 ))

            for scaling in 0.5 2.0 5.0 10.0 20.0 40.0 80.0
            do  
                # Fragile
                if [ $SLURM_ARRAY_TASK_ID -eq $i ]
                then
                    ./exe -instancefile "../data/instances/${instance}${seed}_mod.txt" -follower 2 -dep_gen 2 -nb_intv $nbint -scal_param $scaling -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbthinningSAA $nbthinning -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning -method_dep_gen $method
                fi
                (( i = $i +1 ))
            done

            for scaling in 1.0 2.0 5.0 10.0 20.0
            do  
                # Fragile power
                if [ $SLURM_ARRAY_TASK_ID -eq $i ]
                then
                    ./exe -instancefile "../data/instances/${instance}${seed}_mod.txt" -follower 2 -dep_gen 3 -nb_intv $nbint -scal_param $scaling -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbthinningSAA $nbthinning -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning -method_dep_gen $method
                fi
                (( i = $i +1 ))
            done

            for scaling in 1.0 2.0 5.0 10.0 20.0
            do  
                # Strong power
                if [ $SLURM_ARRAY_TASK_ID -eq $i ]
                then
                    ./exe -instancefile "../data/instances/${instance}${seed}_mod.txt" -follower 2 -dep_gen 4 -nb_intv $nbint -scal_param $scaling -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbthinningSAA $nbthinning -nbvalidateproblems $nbvalidateproblems -nbvalidatescenarios $nbvalidatescenarios -nbvalidatethinning $nbvalidatethinning -method_dep_gen $method
                fi
                (( i = $i +1 ))
            done
        done
    done
done
done
done