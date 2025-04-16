# this script gets all available Simulations from SS path and writes the full dataset to a table for putting arguments into condor submit scripts 

import os
import glob
import re
import pandas as pd
from processor_functions import *
import argparse
import yaml

def GetParser():
    """Argument parser for reading scans script.."""
    parser = argparse.ArgumentParser(
        description="Processing simulations command line options."
    )
    parser.add_argument(
        "--userconfig",
        "-c",
        type=str,
        required=True,
        help="Specify the config for the user e.g. paths",
    )

    return parser

parser = GetParser()
args = parser.parse_args()

user_config_path = args.userconfig
user_config = load_config(user_config_path)
print("Using user config: ", user_config)
SS_path = user_config['SS_path']
llp_id = user_config['llp_id']
scanfile = user_config["scanfile"]

existing_df = pd.read_csv(scanfile, index_col=0)
print("a")
print(existing_df)
if not existing_df.empty:
    max_index = existing_df.index.max()
else:
    max_index = -1 # starting from 0 if csv is empty 

# check existing_df for the SSX_prod_dec names, list their unique instances
existing_SSXs_array = existing_df['SS'].unique()
existing_SSXs_list = existing_SSXs_array.tolist()
print("Found existing SSXs in csv: ")
print(existing_SSXs_list)

# Find all SSXs in SS_path and save to list called SSXs
SSXs = []
SSXs_already_added = []
SSXs_to_add = []
#df_allscans = pd.DataFrame
for subdirpath, subdirnames, subfilenames in os.walk(SS_path):
    for subdirname in subdirnames:
        if subdirname.startswith("SS"):
            SSXs.append(subdirname)
            if subdirname in existing_SSXs_list:
                SSXs_already_added.append(subdirname)
                print(f"SSX {subdirname} already in csv")
            else:
                SSXs_to_add.append(subdirname)
                print(f"SSX {subdirname} not in csv. Append!")
print("Summary: ")
print("Simulation Swarms found in SS_path: ")
print(SSXs)
print(type(SSXs))
print("SSs matched with contents of csv already: ")
print(SSXs_already_added)
print(" ... and new SSXs to add to csv are: ")
print(SSXs_to_add)

# Get scans saved in available SSXs and convert to a lookup table to store data for SampleY's 
#for SSX in SSXs:
new_df = scan_to_df(SS_path, SSXs_to_add, llp_id)
print("New df:")
print(new_df)

new_indices = range(max_index + 1, max_index + 1 + len(new_df))
new_df.index = new_indices

print("existing df:")
print(existing_df)

updated_df = pd.concat([existing_df, new_df])
print("Updated df:")
print(updated_df)
#temp solution
updated_df = new_df

updated_df.to_csv(scanfile, index=True)
