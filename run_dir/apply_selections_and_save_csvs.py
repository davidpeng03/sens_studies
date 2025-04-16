# acceptance

import math
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
#import ROOT
#from root_numpy import root2array, tree2array
import uproot
import awkward as ak
from particle import PDGID, Particle, Charge
from Common_LLP_processing import *
import h5py
import seaborn as sns
import pickle
from colour import Color
sns.set_style("white")
import glob
import re
from processor_functions import *
import time

start_total = time.time()

def GetParser():
    """Argument parser for reading Ntuples script."""
    parser = argparse.ArgumentParser(
        description="Processing simulations command line options."
    )
    parser.add_argument(
        "--SampleY",
        "-Y",
        type=str,
        required=True,
        help="Specify the unique ID of the sample",
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
SampleY = int(args.SampleY)
print("Processing save objects as h5 for Sample", SampleY)
user_config_path = args.userconfig
user_config = load_config(user_config_path)
print("Using user config: ", user_config)
h5_path = user_config['h5_path']
nevents = user_config['nevents']
SS_path = user_config['SS_path']
llp_pid = user_config['llp_id']
minJetPt = user_config['minJetPt']
minChargedPt = user_config['minChargedPt']
scanfile = user_config['scanfile'] # this is samples_to_process.csv
llp_type = user_config['llp_type']
print("LLP type: ", llp_type)
print("SampleY to process and save as h5: ", SampleY)
# set parameters which correspond to this sampleY, logged in the scanfile
samples_to_process_df = pd.read_csv(scanfile)
llp_mass = samples_to_process_df.iloc[SampleY]['mass#9900012']
llp_width = samples_to_process_df.iloc[SampleY]['width#9900012']
llp_Ve = samples_to_process_df.iloc[SampleY]['numixing#1']
llp_id_fromcsv = samples_to_process_df.iloc[SampleY]['llp_id'] 
cross_sec = samples_to_process_df.iloc[SampleY]['cross']
SSX = samples_to_process_df.iloc[SampleY]['SS']
run_name = samples_to_process_df.iloc[SampleY]['#run_name']

# to cross check with llp_id read from config file (if they don't match, this is an error
if (llp_id_fromcsv != llp_pid):
    print("ERROR. LLP ID mismatch between scan file and expectation set in config file.")

# Read in prod_decay from name of SS directory. Each SS will be for only one prod_decay. 
print("SSX is: ", SSX)
parts_SSX = SSX.split('_',1)
prod_decay = parts_SSX[1]
print("Production and decay mode is: ", prod_decay)

# use the above info to access the relevant hepmc file:
filedescription = f"{llp_type}_{prod_decay}_mn1_{llp_mass}_ven1_{llp_Ve}_sample{SampleY}".replace('.','pt').replace('-','min')
print("sample description: ", filedescription)

start_bigloop = 0.
start_smallloop = 0.
start_tinyloop = 0.
end_bigloop = 0.
end_smallloop = 0.
end_tinyloop = 0.

# Set random seed
np.random.seed(10)#args.r)

file_llps = h5_path+f'/{SSX}/all_llps_{filedescription}.h5'
data_llps = h5py.File(file_llps,'r')
file_jets= h5_path+f'/{SSX}/jets_{filedescription}.h5'
data_jets = h5py.File(file_jets,'r')
file_charged= h5_path+f'/{SSX}/charged_{filedescription}.h5'
data_charged = h5py.File(file_charged,'r')

ctau_dec_llps = []
llp_decay_pos = []

hbar = 6.582e-25
c_light = 3.0e+8

lifetime = hbar / llp_width
auto_ctau = c_light * lifetime

store_cutflow = True

check_ctau = True
# Switch the following on if you also want to make/add to a txt file for evaluating decay positions/ctaus
save_output_txt = True

make_ctau_plot = False

#with h5py.File(file_llps, 'r') as data_llps:
#    expo_end_vertex_llp_all = data_llps['expo_end_vertex_llp_all']
#    print("Dataset Shape:", expo_end_vertex_llp_all.shape)
#    print("Dataset Contents:", expo_end_vertex_llp_all[:])  # Check the data
#with h5py.File(file_llps, 'r') as data_llps:
#    if 'expo_end_vertex_llp_all' in data_llps:
#        print("Dataset exists!")
#        expo_end_vertex_llp_all = data_llps['expo_end_vertex_llp_all']
#    else:
#        print("Dataset 'expo_end_vertex_llp_all' does not exist.")
        
expo_end_vertex_llp_all = data_llps['expo_end_vertex_llp_all']
#expo_end_vertex_llp_all = np.array(expo_end_vertex_llp_all).squeeze(axis=-1)
x_expo_end_vertex_llps = [row[0] for row in expo_end_vertex_llp_all]
y_expo_end_vertex_llps = [row[1] for row in expo_end_vertex_llp_all]
z_expo_end_vertex_llps = [row[2] for row in expo_end_vertex_llp_all]
r_expo_end_vertex = [math.sqrt(x**2 + y**2 + z**2) for x, y, z in expo_end_vertex_llp_all]
r_end_vertex = r_expo_end_vertex
x_end_vertex_llps = x_expo_end_vertex_llps
y_end_vertex_llps = y_expo_end_vertex_llps
z_end_vertex_llps = z_expo_end_vertex_llps

make_expo_dec_plots = False

if(make_expo_dec_plots == True):

    boost_llp_all = data_llps['boost_llp_all']
    beta_llp_all = data_llps['beta_llp_all']
    
    r_llps = [math.sqrt(x**2 + y**2 + z**2) for x, y, z in expo_end_vertex_llp_all]
    print("r_llps")
    print(r_llps)
    plt.figure(figsize=(12,7))
    #plt.xlim(0,450)
    #plt.hist(x_expo_end_vertex_llps, bins=50, color='green', label='HNLs',density=True)
    sns.histplot(x_expo_end_vertex_llps, bins=60, kde=False)
    plt.xlabel('x end vertex')
    plt.ylabel('frequency')
    plt.title(f'x histogram HNLs, SampleY={SampleY}, mass={llp_mass}, coupling={llp_Ve}')

    plt.figure(figsize=(12,7))
    #plt.xlim(0,450)
    sns.histplot(y_expo_end_vertex_llps, bins=60, kde=False)
    plt.xlabel('y end vertex')
    plt.ylabel('frequency')
    plt.title(f'y histogram HNLs, SampleY={SampleY}, mass={llp_mass}, coupling={llp_Ve}')

    plt.figure(figsize=(12,7))
    #plt.xlim(0,450)
    sns.histplot(z_expo_end_vertex_llps, bins=60, kde=False)
    plt.xlabel('z end vertex')
    plt.ylabel('frequency')
    plt.title(f'z histogram HNLs, SampleY={SampleY}, mass={llp_mass}, coupling={llp_Ve}')

    plt.figure(figsize=(12,7))
    #plt.xlim(0,450)
    plt.hist(r_llps, bins=50, color='green', label='HNLs',density=True)
    #sns.histplot(r_llps, bins=60, kde=False)
    plt.xlabel('r end vertex')
    plt.ylabel('frequency')
    plt.title(f'r HNLs, SampleY={SampleY}, mass={llp_mass}, coupling={llp_Ve}')

    plt.figure(figsize=(12,7))
    #plt.xlim(0,450)
    plt.hist(boost_llp_all, bins=50, color='red', label='HNLs',density=True)
    #sns.histplot(r_llps, bins=60, kde=False)
    plt.xlabel('boost')
    plt.ylabel('frequency')
    plt.title(f'boost HNLs, SampleY={SampleY}, mass={llp_mass}, coupling={llp_Ve}')

    plt.figure(figsize=(12,7))
    #plt.xlim(0,450)
    plt.hist(beta_llp_all, bins=50, color='blue', label='HNLs',density=True)
    #sns.histplot(r_llps, bins=60, kde=False)
    plt.xlabel('beta')
    plt.ylabel('frequency')
    plt.title(f'beta HNLs, SampleY={SampleY}, mass={llp_mass}, coupling={llp_Ve}')

    plt.show()

if (check_ctau == True):
    
    ctau_dec_llps = []
    llp_decay_pos = []

    hbar = 6.582e-25
    c_light = 3.0e+8

    lifetime = hbar / llp_width
    auto_ctau = c_light * lifetime

    avg_ctau = 0.
    avg_decay_pos = 0.

    r_expo_end_vertex = [math.sqrt(x**2 + y**2 + z**2) for x, y, z in expo_end_vertex_llp_all]
    
    for j in range(len(data_llps['event_number_llp_all'])):
        r = r_expo_end_vertex[j]
        #(np.sqrt(np.square(data_llps['expo_end_vertex_llp_all'][j][0]) + np.square(data_llps['expo_end_vertex_llp_all'][j][1]) + np.square(data_llps['expo_end_vertex_llp_all'][j][2])))
        llp_decay_pos.append(r)
        beta = data_llps['beta_llp_all'][j]
        gamma = data_llps['boost_llp_all'][j]
#        gamma = 1./math.sqrt(1. - beta**2)
        ctau = r / (beta*gamma)
        ctau_dec_llps.append(ctau)

    avg_ctau = sum(ctau_dec_llps) / len(ctau_dec_llps)
    avg_decay_pos = sum(llp_decay_pos) / len(llp_decay_pos)
    
    print("average ctau of llps is ", avg_ctau)
    print("c.f. madspin auto ctau is ", auto_ctau)

    if (make_ctau_plot == True):

        plt.figure(figsize=(12, 7))
        fig_ctau_min_edge = min(ctau_dec_llps)
        fig_ctau_max_edge = max(ctau_dec_llps)
        #    fig_ctau_max_edge = max(max(ctau_expo_end_vertex_all_particles), max(r_expo_end_vertex_llps))
        #fig_r_bins = 10
        fig_ctau_bins = np.linspace(fig_ctau_min_edge, fig_ctau_max_edge, 30)
        plt.hist(ctau_dec_llps, bins=fig_ctau_bins, color='blue', label='HNLs',density=True)
        #    plt.hist(r_expo_end_vertex_llps, bins=fig_r_bins, color='orange', label='HNLs', density=True)
        plt.axvline(avg_ctau, color='r', linestyle='dashed', linewidth=1,label='hepmc mean ctau')
        plt.axvline(auto_ctau, color='g', linestyle='dashed', linewidth=1,label='MG auto ctau')
        plt.xlabel('ctau (decay) (m)')
        plt.ylabel('Density')
        plt.title(f'Decay ctau of HNLs (VeN1={llp_Ve}, mN1={llp_mass})')
        plt.legend()
        plt.show()

    if (save_output_txt == True):

        llp_Ve_label = str(llp_Ve).replace('.','p')

        output_file = f"ctau_coupling{llp_Ve_label}_{prod_decay}.txt"
        file_exists = os.path.isfile(output_file)

        with open(output_file, "a") as file:
            # If the file does not exist, write the headers first
            if not file_exists:
                file.write("mass,hepmc_avgdecpos,hepmc_avgctau,madspin_autoctau,VeN1,process,SampleY\n")
            file.write(f"{llp_mass},{avg_decay_pos},{avg_ctau},{auto_ctau},{llp_Ve},{prod_decay},{SampleY}\n")

if (store_cutflow == True):

    nominal_jet_production_vertex = [origin[0],0,0] # assuming IP is corresponding to the eta,phi assigned by fastjet

    n_llps_all = 0
    n_llps_pass_geometry = 0
    n_llps_pass_geometry_met = 0
    n_llps_pass_geometry_met_deltaRjets = 0
    n_llps_pass_geometry_met_deltaRjets_deltaRcharged = 0
    n_llps_all_weighted = 0
    n_llps_pass_geometry_weighted = 0
    n_llps_pass_geometry_met_weighted = 0
    n_llps_pass_geometry_met_deltaRjets_weighted = 0
    n_llps_pass_geometry_met_deltaRjets_deltaRcharged_weighted = 0
      
    # loop over event
    # 1. check if llp intersects ceiling (geometry cut)
    # 2. if so, check if jets from that event intersect ceiling using ceiling_stat(eta, phi, decay_pos) 
    # 3. within this jet loop, if jet intersects ceiling, calculate intersection point using ceil_intersection(eta, phi, decay_pos) 

    for j in range(len(data_llps['event_number_llp_all'])):
        print("Event ", j)
        start_bigloop = time.time()
#        print("llp end vertex position (metres) ", data_llps['expo_end_vertex_llp_all'][j][0], ", ", data_llps['expo_end_vertex_llp_all'][j][1], ", ", data_llps['expo_end_vertex_llp_all'][j][2]) #@ADDED
#        if (data_llps['expo_end_vertex_llp_all'][j][0] != -1 and data_llps['n_decay_products_llp_all'][j] > 1): # shouldn't I also make sure the llp's status is the one that's long-lived here? Although I already do have the ndecayproducts>1 requirement -- which alone should already be the 1 final state LLP
        if(data_llps['n_decay_products_llp_all'][j] > 1):
#            r = (np.sqrt(np.square(data_llps['expo_end_vertex_llp_all'][j][0]) + np.square(data_llps['expo_end_vertex_llp_all'][j][1]) + np.square(data_llps['expo_end_vertex_llp_all'][j][2])))
#            print(f"r (mass = {masslabel}: ", r)
            decay_pos = (list(to_cartesian(r_expo_end_vertex[j], data_llps['eta_llp_all'][j], data_llps['phi_llp_all'][j])))
            # Add ATLAS offset:
            decay_pos[0] += origin[0]
            print("decay_pos = ", decay_pos)
            
            #should I use the atlas offset in the in_cavern function? 
            llp_is_in_cavern = in_cavern(r_expo_end_vertex[j], data_llps['eta_llp_all'][j], data_llps['phi_llp_all'][j]) # in cavern already has atlas x axis offset included, as it's not added to the r_expo in this script 
#            llp_is_in_cavern = True # for testing 
            llp_is_in_shaft = in_shaft(r_expo_end_vertex[j], data_llps['eta_llp_all'][j], data_llps['phi_llp_all'][j]) # we don't use the shaft 
            n_llps_all+=1
            n_llps_all_weighted+=data_llps['weight_llp_all'][j]
            llp_is_in_ceiling_direction = ceiling_stat(data_llps['eta_llp_all'][j], data_llps['phi_llp_all'][j], decay_pos) # this input contains the atlas offset and the function does not 
            llp_ceil_int = ceil_intersection(data_llps['eta_llp_all'][j], data_llps['phi_llp_all'][j], decay_pos) # this input contains the atlas offset and the function does not 
            print("llp_ceil_int = ", llp_ceil_int)
            if (llp_is_in_cavern):
                print("llp is in cavern")

            if (llp_is_in_ceiling_direction):
                print("llp is in ceiling direction")
                
            if(llp_is_in_cavern and llp_is_in_ceiling_direction):# and ceiling_stat(data_llps['eta_llp_all'], data_llps['phi_llp_all'], data_llps['end_vertex_llp_all'])):
                print("in cavern llp in ceiling direction!!!")
                print("llp end vertex pos from r_end_vertex (m) ", r_end_vertex[j])
                cartesian_r_end =  list(to_cartesian(r_end_vertex[j], data_llps['eta_llp_all'][j], data_llps['phi_llp_all'][j]))
                print("   and (x,y,z) from phi,eta conversion of r_end_vertex: ", cartesian_r_end)
                print("   and (x,y,z) to compare from original end_vertex_llp_all: ", x_end_vertex_llps[j], ", ", y_end_vertex_llps[j], ", ", z_end_vertex_llps[j])
                r_to_compare = math.sqrt((x_end_vertex_llps[j])**2 + (y_end_vertex_llps[j])**2 + (z_end_vertex_llps[j])**2) 
                print("   which has an r value (to compare): ", r_to_compare)
                print("(difference from phi fold)")
#                if(ceiling_stat(data_llps['eta_llp_all'][j], data_llps['phi_llp_all'][j], cartesian_r_expo)): # not sure if this is the correct input format for this function @TODO
#                   print("and passed ceiling stat too")
                # adding check here for intersection with ceiling. This could overestimate the number of llps anubis is sensitive to, as their decay products may spray at larger angles 
                n_llps_pass_geometry+=1
                n_llps_pass_geometry_weighted+=data_llps['weight_llp_all'][j]

                # Check if llp passes met cut
                met_llp = data_llps['met_llp_all'][j]
#                print("!!!!!!!!!!!!! met of llp that passes geometry cut = ", met_llp)
                if(met_llp > min_met):
                    n_llps_pass_geometry_met+=1
                    n_llps_pass_geometry_met_weighted+=data_llps['weight_llp_all'][j]                
                    # Record jets for plotting
                    jets = []

                    llp_n_ceiling = 0
                    llp_n_px14 = np.full(len(anubis_ts_height),0)
                
                    # check if jet isolation criteria passed for this llp's event number
                    deltaR_llp_alljets = []
                    for jet_index in range(len(data_jets['event_number_jets_all'])):
                        start_smallloop = time.time()
                        if (data_jets['event_number_jets_all'][jet_index] == data_llps['event_number_llp_all'][j]):
#                            print("found jets in the same event, number ", data_llps['event_number_llp_all'][j])
                            # for those jets identified in the same event, calculate whether they intersect the ceiling and if so are in range of the LLP:
                            jet_x, jet_y, jet_z = ceil_intersection(data_jets['eta_jets_all'][jet_index], data_jets['phi_jets_all'][jet_index], nominal_jet_production_vertex)

                            # Toby+Jon's functions often assume that the jets are coming from LLPs?

                            # Calculate jet momentum
                            jet_p = data_jets['pt_jets_all'][jet_index] / np.sin(to_theta(data_jets['eta_jets_all'][jet_index]))
                            if (jet_p > jet_min_p):
                                print("jet momentum passes. jet_p = ", jet_p) # already accounted for? 
                                                
                            # deltaR between LLP and any jet in its event cannot be less than min_dRJet
                            # regenerate deltaR_llp_thisjet for each jet, calculate. if the llp+jet are far enough away, append a 0. Appending a 1 means *did not pass*
                            deltaeta_llp_thisjet = abs(data_llps['eta_llp_all'][j] - data_jets['eta_jets_all'][jet_index])
                            deltaphi_llp_thisjet = abs(data_llps['phi_llp_all'][j] - data_jets['phi_jets_all'][jet_index])
                            deltaR_llp_thisjet = np.sqrt(np.square(deltaeta_llp_thisjet) + np.square(deltaphi_llp_thisjet))
                            if(deltaR_llp_thisjet > min_dRJet):
                                deltaR_llp_alljets.append(0)
                            else:
                                deltaR_llp_alljets.append(1)
                            # <= max_nCharged
                        end_smallloop = time.time()

                    if (sum(deltaR_llp_alljets) == 0):
                        start_tinyloop = time.time()
#                        print("Passed geometry+MET+deltaRjets")
                        n_llps_pass_geometry_met_deltaRjets+=1
                        n_llps_pass_geometry_met_deltaRjets_weighted+=data_llps['weight_llp_all'][j]
                        # now final cut: charged particle track isolation
                        deltaR_llp_allcharged = [] # iterate through all charged particles in the event and append a 1 if required condition NOT passed 
                        for charged_index in range(len(data_charged['event_number_charged_all'])):
                            if (data_charged['event_number_charged_all'][charged_index] == data_llps['event_number_llp_all'][j]):
#                                print("found charged energetic particles in the same event, number ", data_llps['event_number_llp_all'][j])

                                deltaeta_llp_thischarged = abs(data_llps['eta_llp_all'][j] - data_charged['eta_charged_all'][charged_index])
                                deltaphi_llp_thischarged = abs(data_llps['phi_llp_all'][j] - data_charged['phi_charged_all'][charged_index])
                                deltaR_llp_thischarged = np.sqrt(np.square(deltaeta_llp_thischarged) + np.square(deltaphi_llp_thischarged))
                                min_dRCharged = min_dRJet # Both should be 0.5
                                if(deltaR_llp_thischarged > min_dRCharged):
                                    deltaR_llp_allcharged.append(0)
                                else:
                                    deltaR_llp_allcharged.append(1)
                                # <= max_nCharged  # this was a different criterion in br_sensitivity.py that I can't find in the literature???
                        end_tinyloop = time.time()
                        if (sum(deltaR_llp_allcharged) == 0): # the llp is far enough away from all charged particles
                            n_llps_pass_geometry_met_deltaRjets_deltaRcharged+=1
                            n_llps_pass_geometry_met_deltaRjets_deltaRcharged_weighted+=data_llps['weight_llp_all'][j]


    end_bigloop = time.time()              
                                     
    # collect only those LLPs which decay within the cavern and intersect the ceiling. Then collect those jets from the same events that also intersect the ceiling. Attribute each to its event number. Store all in one large list with llp or jet label. Different list for each mass point.
    print("n_llps_all = ", n_llps_all)
    print("n_llps_pass_geometry = ", n_llps_pass_geometry)
    print("n_llps_pass_geometry_met = ", n_llps_pass_geometry_met)
    print("n_llps_pass_geometry_met_deltaRjets = ", n_llps_pass_geometry_met_deltaRjets)
    print("n_llps_pass_geometry_met_deltaRjets_deltaRcharged = ", n_llps_pass_geometry_met_deltaRjets_deltaRcharged)
    print("n_llps_all_weighted = ", n_llps_all_weighted)
    print("n_llps_pass_geometry_weighted = ", n_llps_pass_geometry_weighted)
    print("n_llps_pass_geometry_met_weighted = ", n_llps_pass_geometry_met_weighted)
    print("n_llps_pass_geometry_met_deltaRjets_weighted = ", n_llps_pass_geometry_met_deltaRjets_weighted)
    print("n_llps_pass_geometry_met_deltaRjets_deltaRcharged_weighted = ", n_llps_pass_geometry_met_deltaRjets_deltaRcharged_weighted)

    # For cutflows: need 2 lists: cutnames (label for x axis of each cut) and passedcuts (number of events passing each cut cumulatively). Can do this from pickled dictionary.
    cuts_dict = {'All': n_llps_all, 'Cavern decay': n_llps_pass_geometry, 'MET': n_llps_pass_geometry_met, 'DeltaR(LLP,jets) > 0.5': n_llps_pass_geometry_met_deltaRjets, 'DeltaR(LLP,charged) > 0.5': n_llps_pass_geometry_met_deltaRjets_deltaRcharged, 'All weighted': n_llps_all_weighted, 'Cavern decay weighted': n_llps_pass_geometry_weighted, 'MET weighted': n_llps_pass_geometry_met_weighted, 'DeltaR(LLP,jets) > 0.5 weighted': n_llps_pass_geometry_met_deltaRjets_weighted, 'DeltaR(LLP,charged) > 0.5 weighted': n_llps_pass_geometry_met_deltaRjets_deltaRcharged_weighted}
    print("cuts dictionary: ", cuts_dict)
    cuts_dict_names = ['All', 'Cavern decay', 'MET', 'DeltaR(LLP,jets) > 0.5', 'DeltaR(LLP,charged) > 0.5', 'All weighted', 'Cavern decay weighted', 'MET weighted', 'DeltaR(LLP,jets) > 0.5 weighted', 'DeltaR(LLP,charged) > 0.5 weighted']

    csv_out_dir = f"output/{SSX}"
    if not os.path.exists(csv_out_dir):
        os.makedirs(csv_out_dir)
        print(f"Directory {csv_out_dir} created.")
        
    with open(f'{csv_out_dir}/cuts_dict_{filedescription}.pkl', 'wb') as fp:
        pickle.dump(cuts_dict, fp)
        print(f'dictionary saved to pkl file for {filedescription}')

    with open(f'{csv_out_dir}/cuts_dict_{filedescription}.csv', 'w') as csvfile:
        csvwriter = csv.DictWriter(csvfile, fieldnames=cuts_dict_names)#, delimiter=' ')
        csvwriter.writeheader()
        csvwriter.writerow(cuts_dict)
        print(f'dictionary saved to csv file for {filedescription}')

end_total = time.time()

print(f"Script total time: {measure_time(start_total, end_total)} seconds")
print(f"       Big loop took: {measure_time(start_bigloop, end_bigloop)} seconds")
print(f"       Small loop took: {measure_time(start_smallloop, end_smallloop)} seconds")
print(f"       Tiny loop took: {measure_time(start_tinyloop, end_tinyloop)} seconds")

