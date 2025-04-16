
#input arguments: {index} {config_file}
# (userconfig path, unique SampleY ID number.)

# @TODO still to add:
# --- read event weight from hepmc file (should have eta cut weights) -- DONE
# --- appropriate names for h5 files -- DONE
# --- phi folding by adding pi to the phi of all saved objects -- DONE
# --- create directories that do not already exist 
from numpy.testing import assert_array_equal
from numpy.lib.recfunctions import append_fields
#from pyjet.testdata import get_event
import argparse
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from constants import *
from helpers import *
from plot import *
from geometry import *
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from colour import Color 
#import ROOT
#from root_numpy import root2array, tree2array
import uproot
import awkward as ak
from particle import PDGID, Particle, Charge
from Common_LLP_processing import *
#import pythia8
import vector
import scipy as sp
import scipy.optimize
import pyhepmc
#import fastjet # use anti-kt algorithm 
import fastjet._pyjet
import pytest
#import pyjet 
import h5py
import os
import fnmatch
import glob
import re
from processor_functions import *
import math
import random

user_config_path = "/usera/dp728/run_dir/config_hnl_madgraph.yaml"
"""
print("Using user config: ", user_config)
h5_path = user_config['h5_path']
nevents = user_config['nevents']
SS_path = user_config['SS_path']
llp_pid = user_config['llp_id']
minJetPt = user_config['minJetPt']
minChargedPt = user_config['minChargedPt']
scanfile = user_config['scanfile'] # this is samples_to_process.csv
llp_type = user_config['llp_type']
minTheta = user_config['minTheta']
maxTheta = user_config['maxTheta']
minEta = -math.log(math.tan(maxTheta/2.))
maxEta = -math.log(math.tan(minTheta/2.))

print("LLP type: ", llp_type)
# set parameters which correspond to this sampleY, logged in the scanfile

samples_to_process_df = pd.read_csv(scanfile)
llp_mass = samples_to_process_df.iloc[SampleY]['mass#9900012']
llp_width = samples_to_process_df.iloc[SampleY]['width#9900012']
llp_Ve = samples_to_process_df.iloc[SampleY]['numixing#1']
llp_id_fromcsv = samples_to_process_df.iloc[SampleY]['llp_id'] 
cross_sec = samples_to_process_df.iloc[SampleY]['cross']
SSX = samples_to_process_df.iloc[SampleY]['SS']
run_name = samples_to_process_df.iloc[SampleY]['#run_name']
"""

print("a")

bigarrayforgraphs = []

