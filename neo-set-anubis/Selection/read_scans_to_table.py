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

# Find all SSXs in SS_path and save to list called SSXs
SSXs = []
df_allscans = pd.DataFrame
for subdirpath, subdirnames, subfilenames in os.walk(SS_path):
    for subdirname in subdirnames:
        if subdirname.startswith("SS"):
            SSXs.append(subdirname)
print("Simulation Swarms found in SS_path: ")
print(SSXs)
print(type(SSXs))
# Get scans saved in available SSXs and convert to a lookup table to store data for SampleY's 
#for SSX in SSXs:
scans_df_combined = scan_to_df(SS_path, SSXs, llp_id) 
print("Combined dataframe of scans:")
print(scans_df_combined)

scans_df_combined.to_csv("samples_to_process.csv")#, index=False)
