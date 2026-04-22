import os

general_folder    = "../results/"
behaviors         = [0, 1, 2] 
folder_behaviors  = ["fixstrongweak/","depstrongweak/","depgeneral/"]
folder_sol_types  = { 0: [""],
                      1: ["saa/","dep/"],
                      2: [""]}
folder_spec_types = { 0: [""], 
                      1: ["proportional/","threshold/","strong/","fragile/","strong_power/","fragile_power/"], 
                      2: ["neutral/","proportional/","fragile/","fragile_power/","strong_power/"]}

# List of folder names to create
folders = [(general_folder + folder_behaviors[bv] + sol + tp) for bv in behaviors for sol in folder_sol_types[bv] for tp in folder_spec_types[bv]]
folders.append("../results/compiled/")

# Loop through the folder names and create each one
for folder in folders:
    try:
        os.makedirs(folder, exist_ok=True)
        print(f"Folder '{folder}' created successfully.")
    except FileExistsError:
        print(f"Folder '{folder}' already exists.")