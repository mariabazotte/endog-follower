import math
import os

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

def compute_norm(values):
    return math.sqrt(sum(v*v for v in values))

def process_instance(file_path, instance_type, seed, output_file):
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    # --- Header ---
    line_idx = 0
    nb_vars_leader, nb_constrs_leader, _, _ = map(float, lines[line_idx].split())
    nb_vars_leader = int(nb_vars_leader)
    line_idx += 1

    nb_vars_follower, nb_constrs_follower, _, _ = map(float, lines[line_idx].split())
    nb_vars_follower = int(nb_vars_follower)
    line_idx += 1

    # --- Leader variables ---
    leader_obj = []
    for _ in range(nb_vars_leader):
        parts = lines[line_idx].split()
        obj = float(parts[3])
        leader_obj.append(obj)
        line_idx += 1

    # --- Follower variables ---
    follower_obj_leader = []
    for _ in range(nb_vars_follower):
        parts = lines[line_idx].split()
        obj_leader = float(parts[3])
        follower_obj_leader.append(obj_leader)
        line_idx += 1

    # --- Compute norms ---
    norm_leader = compute_norm(leader_obj)
    norm_follower = compute_norm(follower_obj_leader)

    total = norm_leader + norm_follower
    perc_leader = norm_leader / total if total > 0 else 0
    perc_follower = norm_follower / total if total > 0 else 0

    # --- Write output ---
    with open(output_file, 'a') as f:
        f.write(f"{instance_type};{seed};{norm_leader:.6f};{norm_follower:.6f};{perc_leader:.6f};{perc_follower:.6f};\n")


def main():
    output_file = "instances/info.txt"
    with open(output_file, 'w') as f:
        f.write("instance;seed;norm_leader;norm_follower;perc_leader;perc_follower;\n")

    for instance in instances:
        for seed in seeds:
            input_file = "instances/" + instance + str(seed) + ".txt"
            process_instance(input_file, instance, seed, output_file)

if __name__ == "__main__":
    main()