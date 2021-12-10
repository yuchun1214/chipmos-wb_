import main
import json
import os
import re
import random
from shutil import copy
import pandas as pd

def _double_the_recipe(lots, machines):
    recipes = list(set(lots.recipe))
    rows = []
    for i in range(len(recipes)):
        recipe = recipes[i]
        same_kinds_of_lots = lots[lots["recipe"] == recipe]
        same_kinds_of_lots.reset_index(drop=True, inplace=True)
        for j in range(len(same_kinds_of_lots.values)):
            same_kinds_of_lots.at[j, "recipe"] += "-" + str(j % 2)
            rows.append(same_kinds_of_lots.iloc[j])
    new_lot_data = pd.DataFrame(rows, columns=lots.columns)
    
    rows = []
    recipes = list(set(machines["BOND ID"]))
    for i in range(len(recipes)):
        recipe = recipes[i]
        same_kinds_of_machines = machines[machines["BOND ID"] == recipe]
        same_kinds_of_machines.reset_index(drop=True, inplace=True)
        for j in range(len(same_kinds_of_machines.values)):
            same_kinds_of_machines.at[j, "BOND ID"] += "-" + str(j % 2)
            rows.append(same_kinds_of_machines.iloc[j])
    new_machine_data = pd.DataFrame(rows, columns=machines.columns)
    
    return new_lot_data, new_machine_data

def double_the_recipe(directory):
    lots = pd.read_csv(os.path.join(directory_name, "lots.csv"))
    machines = pd.read_csv(os.path.join(directory_name, "machines.csv"))
    new_lot_data, new_machine_data = _double_the_recipe(lots, machines)
    new_lot_data.to_csv(os.path.join(directory_name, "nlots.csv"), index=False)
    new_machine_data.to_csv(os.path.join(directory_name, "nmachines.csv"), index=False)

def _insufficient_machines(locations, machines):
    locs = list(set(locations.Location))
    ratio = 0.2
    remaining_machines = []
    for loc in locs:
        entities = list(locations[locations["Location"] == loc]["Entity"])

        size_of_machines = len(entities)
        number_of_machines = int((1-ratio)*size_of_machines)

        remaining_machines += random.sample(entities, number_of_machines)
    rows = []
    for machine in remaining_machines:
        rows += list(machines[machines["ENTITY"] == machine].values)
    new_machines = pd.DataFrame(rows, columns=machines.columns)
    return new_machines

def insufficient_machines(directory_name:str):
    locations = pd.read_csv(os.path.join(directory_name, "locations.csv"))
    machines = pd.read_csv(os.path.join(directory_name, "machines.csv"))
    new_machines = _insufficient_machines(locations, machines)
    new_machines.to_csv(os.path.join(directory_name, "nmachines.csv"))

def double_recipe_and_insufficient_machines(directory_name):
    lots = pd.read_csv(os.path.join(directory_name, "lots.csv"))
    machines = pd.read_csv(os.path.join(directory_name, "machines.csv"))
    locations = pd.read_csv(os.path.join(directory_name, "locations.csv"))
    
    new_machines = _insufficient_machines(locations, machines)
    
    print(new_machines.columns)
    new_lot_datap, new_machine_data = _double_the_recipe(lots, new_machines)
    
    new_lot_data.to_csv(os.path.join(directory_name, "nlots.csv"), index=False)
    new_machines.to_csv(os.path.join(directory_name, "nmachines.csv"), index=False)

def normal(directory_name:str):
    copy(os.path.join(directory_name, "lots.csv"), os.path.join(directory_name, "nlots.csv"))
    copy(os.path.join(directory_name, "machines.csv"), os.path.join(directory_name, "nmachines.csv"))


def copy_data(config_df):
    for i in range(len(config_df.values)):
        entry = config_df.iloc[i]
        output_dir_name = "output_" + entry.no
        copy(entry.machines, output_dir_name)
        copy(entry.locations, output_dir_name)

if __name__ == '__main__' :
    # 1. run the executable file main with the arguments : --file=config.csv -p
    # for just preprocessing files only
    # after the step 1, the relative path will contains several output_* directories
    
    config_df = pd.read_csv("config.csv")

    # 2. copy the data
    # copy the data first
    copy_data(config_df)
    
    

    # experiments_topic_generating_functions = [normal, double_the_recipe, insufficient_machines, double_recipe_and_insufficient_machines]
    experiments_topic_generating_functions = [normal]
    for i in range(len(config_df.values)):
        no = config_df.iloc[i].no

        output_dir_name = "output_" + no
        
        # for each experiment function
        for j in range(len(experiments_topic_generating_functions)):
            experiments_topic_generating_functions[j](output_dir_name)
            
            # run the executable file main
            # ...(record the output)
            # run the executable file main_gpu
            # ...(record the output)

