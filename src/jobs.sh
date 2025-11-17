#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=50G
#SBATCH --time=7:00:00
#SBATCH --array=1-960
#SBATCH --output=../slurm/test_%A_%a.out
#SBATCH --account=ctb-dionneg1

i=1
seed_sce=0
verbose=1
nbthreads=1
timelimit=3600

nbproblems=5
nbscenarios=100 
nbvalidatescenarios=10000

for epsnearopt in 0.05 0.1
do
for instance in 1 2 3 4 5 6 7 8 9 10
do  
    # Fixed Strong-Weak.
    for fixcoop in 0.0 0.3 0.5 0.7 1.0
    do
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            ./exe -instancefile "../data/instances/denegre-miblp_20_15_50_0110_10_${instance}.txt" -follower 0 -fix_strwk $fixcoop -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -nbvalidatescenariosSAA $nbvalidatescenarios
        fi
        (( i = $i +1 ))
    done

    # Decision-dependent Strong-Weak.
    for approach in 0 1
    do
        # Proportional
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            ./exe -instancefile "../data/instances/denegre-miblp_20_15_50_0110_10_${instance}.txt" -follower 1 -dep_strwk 0 -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -seed $seed_sce -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidatescenariosSAA $nbvalidatescenarios
        fi
        (( i = $i +1 ))

        # Threshold
        for param in 0.5 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/denegre-miblp_20_15_50_0110_10_${instance}.txt" -follower 1 -dep_strwk 1 -thr_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -seed $seed_sce -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidatescenariosSAA $nbvalidatescenarios
            fi
            (( i = $i +1 ))
        done

        # Strong/Flexible
        for param in 0.5 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/denegre-miblp_20_15_50_0110_10_${instance}.txt" -follower 1 -dep_strwk 2 -str_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -seed $seed_sce -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidatescenariosSAA $nbvalidatescenarios
            fi
            (( i = $i +1 ))
        done

        # Fragile/Inflexible
        for param in 0.5 2.0 5.0 10.0
        do
            if [ $SLURM_ARRAY_TASK_ID -eq $i ]
            then
                ./exe -instancefile "../data/instances/denegre-miblp_20_15_50_0110_10_${instance}.txt" -follower 1 -dep_strwk 3 -frag_param $param -sol $approach -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -seed $seed_sce -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidatescenariosSAA $nbvalidatescenarios
            fi
            (( i = $i +1 ))
        done
    done

    # Decision-dependent general.

    # Neutral 
    if [ $SLURM_ARRAY_TASK_ID -eq $i ]
    then
        ./exe -instancefile "../data/instances/denegre-miblp_20_15_50_0110_10_${instance}.txt" -follower 2 -dep_gen 0 -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -seed $seed_sce -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidatescenariosSAA $nbvalidatescenarios
    fi
    (( i = $i +1 ))

    for nbint in 4 5
    do 
    for scaling in 0.5 2.0 5.0 10.0
    do  
        # Proportional
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            ./exe -instancefile "../data/instances/denegre-miblp_20_15_50_0110_10_${instance}.txt" -follower 2 -dep_gen 1 -nb_intv $nbint -scal_param $scaling -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -seed $seed_sce -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidatescenariosSAA $nbvalidatescenarios
        fi
        (( i = $i +1 ))

        # Strong Fragile
        if [ $SLURM_ARRAY_TASK_ID -eq $i ]
        then
            ./exe -instancefile "../data/instances/denegre-miblp_20_15_50_0110_10_${instance}.txt" -follower 2 -dep_gen 2 -nb_intv $nbint -scal_param $scaling -near_opt 1 -eps_near_opt $epsnearopt -nbthreads $nbthreads -timelimit $timelimit -verbose $verbose -seed $seed_sce -nbproblemsSAA $nbproblems -nbscenariosSAA $nbscenarios -nbvalidatescenariosSAA $nbvalidatescenarios
        fi
        (( i = $i +1 ))
    done
    done

done
done