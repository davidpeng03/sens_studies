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
#import ROOT
#from root_numpy import root2array, tree2array
import uproot
import awkward as ak
from particle import PDGID, Particle, Charge
from Common_LLP_processing import *
#import pythia8
import vector
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
minTheta = user_config['minTheta']
maxTheta = user_config['maxTheta']
minEta = -math.log(math.tan(maxTheta/2.))
maxEta = -math.log(math.tan(minTheta/2.))

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
llp_file = SS_path+f"/{SSX}/Events/{run_name}_decayed_1/tag_1_pythia8_events.hepmc"

# Check if unzipped hepmc file already exists. If not, then check if zipped hepmc file exists. 
if os.path.exists(llp_file):
    print("sample hepmc file exists.")
elif os.path.exists(llp_file+".gz"):
    print("sample hepmc file exists. unzip now.")
    file_to_unzip = llp_file+".gz"
    unzip_file(file_to_unzip, llp_file)
else:
    print("sample hepmc NOT FOUND. simulation missing.")  

do_llp_decay = True
N_llp = 0
countevts = 0
# Always put data in lists first to grow data, then convert the lists to dataframe. Never grow a df row-wise. Too intensive 
pt_neutrinos = [] 

px_llp_all = [] # element 0 in momentum 4vec in particles vector 
py_llp_all = [] # el 1 in momentum 4vec in particles vector 
pz_llp_all = [] # el 2 in momentum 4vec in particles vector 
E_llp_all = [] # el 3 in probably he 4th mom 4vec element 
#Mass_llp_all = [] # mass double in particles vector 
m_reco_llp_all = []
prod_vertex_pos = [] # trying to find production vertex location (coordinate system?) for all final state particles (required for jets and charged particles to give trajectories along with phi and thereby calculate ceiling intersection)
phi_llp_all = []
eta_llp_all = []
met_llp_all = []
m_llp_all = []
production_vertex_llp_all = []
end_vertex_llp_all = []
status_llp_all = []
boost_llp_all = []
theta_llp_all = []
pt_llp_all = []
beta_llp_all = []
ctau_llp_all = []
in_cavern_llp_all = [] # true or false for decay inside cavern volume 
n_decay_products_llp_all = []
event_number_llp_all = []
children_summass_llp_all = []
weight_llp_all = []

phi_jets_all = []
eta_jets_all = []
theta_jets_all = []
pt_jets_all = []
m_jets_all = []
px_jets_all = []
py_jets_all = []
pz_jets_all = []
E_jets_all = []
production_vertex_jets_all = []
event_number_jets_all = []
weight_jets_all = []

phi_charged_all = []
eta_charged_all = []
theta_charged_all = []
pt_charged_all = []
m_charged_all = []
px_charged_all = []
py_charged_all = []
pz_charged_all = []
E_charged_all = []
production_vertex_charged_all = []
event_number_charged_all = []
weight_charged_all = []

px_j = []
py_j = []
pz_j = []
E_j = []
Mass_j = []
prod_vertex_pos_j = []
countLLP1 = 0

jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
    
with pyhepmc.open(llp_file) as f:
    for event in f:
        countnocharge = 0
        countevts+=1
        Njetparticlesperevt = 0
        Npartperevt = 0
               
        print("event number: ", countevts)

        # append info directly from hepmc file into a df
        evtMET = 0. # initialise MET of the event (for the LLP)
        evtMETpx = 0.
        evtMETpy = 0.
        evt_px_j = [] # initialising jet momenta container for this event
        evt_py_j = []
        evt_pz_j = []
        evt_E_j = []
        evt_pt_j = []
        evt_prod_vert_pos = []
        n_llps_per_event = 0 

        for p in event.particles:
            ###### first make sure this particle has correct eta direction 
#            this_p_eta = calculate_eta(p.momentum[0],p.momentum[1],p.momentum[2])
#            if (this_p_eta <= minEta or this_p_eta >= maxEta):
#                break # discount particle 
            #####################
            Npartperevt+=1
            if(len(p.children) == 0):
                if (p.pid == 5214 or p.pid == 5212 or p.pid == 9940003 or p.pid == 9940011): # don't actually know what particle the 994... are cause it's not on the mc pid list but assuming it's not charged????? # 52xx are neutral bottom baryons.  
                    countnocharge+=1
                    # also check if not produced by LLP (check that mother1 and mother2 are not llp_pid
                    NLLPsinparents = calculate_NLLPsinparents(p,llp_pid)                                        
                    if(NLLPsinparents > 0):
                        print("id'd particle, NLLPsinparents = ", NLLPsinparents)
                    if(NLLPsinparents == 0): # checking that part is not produced by LLP
                        # add neutral particles to jets if they are also prompt:
                        evtMETpx+=p.momentum[0]
#                        print("contribution to MET by particle ", p.pid, ". px contribution = ", p.momentum[0])
                        evtMETpy+=p.momentum[1]
                        if (np.linalg.norm(p.production_vertex.position)<10.):
                            # then add to jets
                            evt_px_j.append(p.momentum[0])
                            evt_py_j.append(p.momentum[1])
                            evt_pz_j.append(p.momentum[2])
                            evt_E_j.append(p.momentum[3])
                            evt_pt_j.append(calculate_pt(p.momentum[0],p.momentum[1]))
                            evt_prod_vert_pos.append(p.production_vertex.position) # add vertex position for each particle's production 
                            Njetparticlesperevt += 1                    
                if (p.pid == 543 or p.pid == -20413): # unrecognised pdgid (in Theo's truc file), but should be a charged B meson (B_c^*+)
                    # also check if not produced by LLP (check that mother1 and mother2 are not llp_pid
                    NLLPsinparents = calculate_NLLPsinparents(p,llp_pid)                                      
                    if(NLLPsinparents > 0):
                        print("id'd particle, NLLPsinparents = ", NLLPsinparents)
                    if(NLLPsinparents == 0): # checking that part is not produced by LLP
                        evtMETpx+=p.momentum[0]
#                        print("contribution to MET by particle ", p.pid, ". px contribution = ", p.momentum[0])
                        evtMETpy+=p.momentum[1]
                        # add neutral particles to jets if they are also prompt:
                        if (np.linalg.norm(p.production_vertex.position)<10.): # <10 counts as "prompt"
                            # then add to jets
                            evt_px_j.append(p.momentum[0])
                            evt_py_j.append(p.momentum[1])
                            evt_pz_j.append(p.momentum[2])
                            evt_E_j.append(p.momentum[3])
                            evt_pt_j.append(calculate_pt(p.momentum[0],p.momentum[1]))
                            evt_prod_vert_pos.append(p.production_vertex.position) # add vertex position for each particle's production 
                            Njetparticlesperevt += 1

                    # final state, charged - now check if also prompt and energetic enough with pt>10GeV, then if so add to charged particle datasets:
                    p_mom_charged = calculate_pt(p.momentum[0], p.momentum[1])
                    if (np.linalg.norm(p.production_vertex.position)<10. and p_mom_charged > minChargedPt):
                        #                    if (p_mom_charged > minChargedPt):
#                        print("Charged particle momentum p_mom_charged = ", p_mom_charged)
                        phi_charged_all.append(calculate_phi(p.momentum[0], p.momentum[1]))
                        eta_charged_all.append(calculate_eta(p.momentum[0],p.momentum[1],p.momentum[2]))
                        weight_charged_all.append(event.weight())

                        theta_charged_all.append(to_theta(calculate_eta(p.momentum[0],p.momentum[1],p.momentum[2])))
                        pt_charged_all.append(p_mom_charged)
                        m_charged_all.append(p.generated_mass)
                        production_vertex_charged_all.append(p.production_vertex.position) # this should be (0,0,0,0) according to criteria above 
                        event_number_charged_all.append(countevts)


                if (p.pid != llp_pid and p.pid != llp_pid * (-1) and p.pid != 543 and p.pid != 9940003 and p.pid != 9940011 and p.pid != -20413 and p.pid != 5212 and p.pid != 5214):
                    #            if (p.pid != llp_pid and p.pid != 5214): ########### Do something with the 5214 to check if it's a jet candidate too. It's a bottom baryon (Sigma_b^{*0}), neutral particle.
                    if (p.pid == 12 or p.pid == 14 or p.pid == 16 or p.pid == 18): # only neutrinos 
                        pt_neutrinos.append(calculate_pt(p.momentum[0],p.momentum[1]))

                    if (p.pid != 12 and p.pid != 14 and p.pid != 16 and p.pid != 18): # from 'everything else': if not neutrinos then add to MET 
                        NLLPsinparents = calculate_NLLPsinparents(p,llp_pid)                                    
                        if(NLLPsinparents > 0):
                            print("not neutrinos or any other id'd particle, NLLPsinparents = ", NLLPsinparents)
                        if(NLLPsinparents == 0): # checking that part is not produced by LLP
                            evtMETpx+=p.momentum[0] 
#                            print("contribution to MET by particle ", p.pid, ". px contribution = ", p.momentum[0])
                            evtMETpy+=p.momentum[1]

                    # add any neutral or charged particles to jets if they are also prompt:
                    if (np.linalg.norm(p.production_vertex.position)<10.):
                        # also check if not produced by LLP (check that parents are not llp_pid
                        NLLPsinparents = calculate_NLLPsinparents(p,llp_pid)                                    
                        if(NLLPsinparents > 0):
                            print("prompt particle, NLLPsinparents = ", NLLPsinparents)
                        if(NLLPsinparents == 0): # checking that part is not produced by LLP
                            # then add to jets
                            evt_px_j.append(p.momentum[0])
                            evt_py_j.append(p.momentum[1])
                            evt_pz_j.append(p.momentum[2])
                            evt_E_j.append(p.momentum[3])
                            evt_pt_j.append(calculate_pt(p.momentum[0],p.momentum[1]))
                            evt_prod_vert_pos.append(p.production_vertex.position) # add vertex position for each particle's production 
                            Njetparticlesperevt += 1
                        
                    part = Particle.from_pdgid(p.pid)
                    if (part.charge != 0): # these are charged and final state and not llps                             
                        # final state, charged - now check if also prompt and energetic enough with pt>10GeV, then if so add to charged particle datasets:
                        p_mom_charged = calculate_pt(p.momentum[0], p.momentum[1])
                        if (np.linalg.norm(p.production_vertex.position)<10. and p_mom_charged > minChargedPt): # 10 should be in mm for the vertex position 
#                            print("!!!!!!!!!@%%%%%%%%%%%%%%% charged particle passed pt and promptness checks")
                            #                            if (p_mom_charged > minChargedPt):
#                            print("Charged particle momentum p_mom_charged = ", p_mom_charged)
                            phi_charged_all.append(calculate_phi(p.momentum[0], p.momentum[1]))
                            eta_charged_all.append(calculate_eta(p.momentum[0],p.momentum[1],p.momentum[2]))
                            weight_charged_all.append(event.weight())
                            theta_charged_all.append(to_theta(calculate_eta(p.momentum[0],p.momentum[1],p.momentum[2])))
                            pt_charged_all.append(p_mom_charged)
                            m_charged_all.append(p.generated_mass)
                            production_vertex_charged_all.append(p.production_vertex.position) # this should be (0,0,0,0) according to criteria above 
                            event_number_charged_all.append(countevts)


            if (p.pid == llp_pid): # here exited the final state particle condition 
                N_llp+=1
                n_llps_per_event+=1
                px_llp_all.append(p.momentum[0])
                py_llp_all.append(p.momentum[1])
                pz_llp_all.append(p.momentum[2])
                E_llp_all.append(p.momentum[3])
                m_llp_all.append(p.generated_mass) ##### help cause this is maybe not the same as mass calc from mom 4vec due to precision issues
                status_llp_all.append(p.status)
                production_vertex_llp_all.append(p.production_vertex.position)
                if (do_llp_decay == False):
                    end_vertex_llp_all.append((-1,-1,-1,-1))
                else: #if (do_llp_decay == True):
                    end_vertex_llp_all.append(p.end_vertex.position)
                phi_llp_all.append(calculate_phi(p.momentum[0], p.momentum[1]))
                eta_llp_all.append(calculate_eta(p.momentum[0],p.momentum[1],p.momentum[2]))
                theta_llp_all.append(to_theta(calculate_eta(p.momentum[0],p.momentum[1],p.momentum[2])))
                boost_llp_all.append(calculate_boost(calculate_pt(p.momentum[0],p.momentum[1]), calculate_eta(p.momentum[0],p.momentum[1],p.momentum[2]), p.generated_mass))
                beta_llp_all.append(calculate_beta(calculate_boost(calculate_pt(p.momentum[0],p.momentum[1]), calculate_eta(p.momentum[0],p.momentum[1],p.momentum[2]), p.generated_mass)))
                pt_llp_all.append(calculate_pt(p.momentum[0],p.momentum[1]))
                n_decay_products_llp_all.append(len(p.children))
                event_number_llp_all.append(countevts)
                weight_llp_all.append(event.weight())
                print("event weight ", event.weight())
#                if(len(p.children) == 0):
                    # add any final state LLPs to the MET for this event 
#                    evtMETpx+=p.momentum[0] 
#                    evtMETpy+=p.momentum[1]
#                    print("LLP is contributing to MET!")
                children_summass = 0
                for child in range(len(p.children)):
                    children_summass += p.children[child].generated_mass
                children_summass_llp_all.append(children_summass)
#                print("children: ", p.children) # why are the llp's children just a bunch of the same llp... ok actually pythia seems to do this so they propagate, changes their status 

# calculate MET for each llp by summing the pt of all the visible particles in each llp's event. visible means charged (electric or strong) 
        print("countnocharge (per event)= ", countnocharge)
        print("MET (per event) mpx and mpy = ", evtMETpx, " and ", evtMETpy)
#        print("evt_px_j length = ", len(evt_px_j))
        print("!!!!!!!!!Njetparticlesperevt = ", Njetparticlesperevt)
        print("!!!!!!!!!Npartperevt = ", Npartperevt)
        evtMET = calculate_pt(evtMETpx, evtMETpy) 
#        met_llp.append(evtMET)
        print("n_llps_per_event = ", n_llps_per_event)
        for llpinevent in range(n_llps_per_event):
            met_llp_all.append(evtMET)
        evt_prod_vert_pos_x, evt_prod_vert_pos_y, evt_prod_vert_pos_z, evt_prod_vert_pos_t = zip(*evt_prod_vert_pos)
        evt_jetsarray = np.array(list(zip(evt_px_j, evt_py_j, evt_pz_j, evt_E_j, evt_prod_vert_pos_t, evt_prod_vert_pos_x, evt_prod_vert_pos_y, evt_prod_vert_pos_z,evt_pt_j)),
                        dtype=[('px', float), ('py', float), ('pz', float), ('E', float), ('decay_pos_jet_t', float), ('decay_pos_jet_x', float), ('decay_pos_jet_y', float), ('decay_pos_jet_z', float),("pt",float)])
        evt_jetsarray = ak.from_numpy(evt_jetsarray)        
#        sequence = fastjet.ClusterSequence(evt_jetsarray, jetdef)# R=0.4, p=-1)#R=0.6, p=-1)
        sequence = fastjet._pyjet.AwkwardClusterSequence(evt_jetsarray, jetdef)# R=0.4, p=-1)#R=0.6, p=-1)
#        jet_cluster = fastjet.ClusterSequence(evt_jetsarray, jetdef)

        jets = sequence.inclusive_jets()
 
        print("length of jets = ", len(jets))
        print("***************************** loop end")

        for i, jet in enumerate(jets):
            thisjet_pt = np.sqrt((jet.px)**2 + (jet.py)**2)
            thisjet_spherical = list(to_spherical(jet.px, jet.py, jet.pz)) # returns r, eta, phi
            thisjet_eta = thisjet_spherical[1]
            thisjet_phi = thisjet_spherical[2]
            
#            jet_p = jet.pt / np.sin(to_theta(jet.eta))
            jet_p = thisjet_pt / np.sin(to_theta(thisjet_eta))
            if (jet_p > jet_min_p and thisjet_pt > minJetPt): ### @ TODO check which of these cuts to use - p or pt               
                phi_jets_all.append(thisjet_phi)
                event_number_jets_all.append(countevts)
                eta_jets_all.append(thisjet_eta)
                weight_jets_all.append(event.weight())
                pt_jets_all.append(thisjet_pt)
#                m_jets_all.append(jet.mass)
#                theta_jets_all.append(to_theta(jet.eta))
                px_jets_all.append(jet.px)
                py_jets_all.append(jet.py)
                pz_jets_all.append(jet.pz)
                E_jets_all.append(jet.E)
            
print("EXITING EVENT LOOP")

h5_out_dir = h5_path+f"/{SSX}"
if not os.path.exists(h5_out_dir):
    os.makedirs(h5_out_dir)
    print(f"Directory {h5_out_dir} created.")

print("eta")
print(eta_llp_all)
print("phi")
print(phi_llp_all)

phi_llp_all = [this_phi + math.pi if this_phi < 0 else this_phi for this_phi in phi_llp_all]
print("new phi:")
print(phi_llp_all)

print("eta values: ")
print(eta_llp_all)

#print("initial weights from mg before phi fold:")
#print("           llps: ", weight_llp_all)
#print("           jets: ", weight_jets_all)
#print("           charged: ", weight_charged_all)

weight_llp_all = [this_weight*0.5 for this_weight in weight_llp_all]
weight_jets_all = [this_weight*0.5 for this_weight in weight_jets_all]
weight_charged_all = [this_weight*0.5 for this_weight in weight_charged_all]
#print("new weights after phi fold:")
#print("           llps: ", weight_llp_all)
#print("           jets: ", weight_jets_all)
#print("           charged: ", weight_charged_all)
print("min eta: ", minEta)
print("max eta: ", maxEta)
count_pass_eta_llp = sum(1 for eta in eta_llp_all if minEta <= eta <= maxEta)
print("number of llp entries passing eta cut... ", count_pass_eta_llp)

count_pass_eta_jets = sum(1 for eta in eta_jets_all if minEta <= eta <= maxEta)
print("number of jets entries passing eta cut... ", count_pass_eta_jets)

count_pass_eta_charged = sum(1 for eta in eta_charged_all if minEta <= eta <= maxEta)
print("number of charged entries passing eta cut... ", count_pass_eta_charged)

pass_eta_charged_events=[]
pass_eta_jets_events=[]
pass_eta_llp_events=[]

for i, eta in enumerate(eta_llp_all):
    if minEta <= eta <= maxEta:
        this_event = event_number_llp_all[i]
        pass_eta_llp_events.append(this_event)

for i, eta in enumerate(eta_jets_all):
    if minEta <= eta <= maxEta:
        this_event = event_number_jets_all[i]
        pass_eta_jets_events.append(this_event)

for i, eta in enumerate(eta_charged_all):
    if minEta <= eta <= maxEta:
        this_event = event_number_charged_all[i]
        pass_eta_charged_events.append(this_event)

print("Unique instances:")
print("charged: ",len(set(pass_eta_charged_events)))
print("llp: ",len(set(pass_eta_llp_events)))
print("jets: ",len(set(pass_eta_jets_events)))
unique_charged = list(set(pass_eta_charged_events))
unique_jets = list(set(pass_eta_jets_events))
unique_llp = list(set(pass_eta_llp_events))
print("sets: ")
print("charged: ",unique_charged)
print("llp: ",unique_llp)
print("jets: ",unique_jets)

print("list of ALL unique instances: ")
print(len(unique_charged + unique_jets + unique_llp))
print("number of unique instances in total with no repetitions from overlap between llp/jets/charged: ")
print(len(set(unique_charged + unique_jets + unique_llp)))
print(set(unique_charged + unique_jets + unique_llp))

with pyhepmc.open(llp_file) as f:
    event_count_tot = sum(1 for _ in f)

print(f"Number of events in {llp_file} is {event_count_tot}")

with h5py.File(h5_out_dir+f'/all_llps_{filedescription}.h5', 'w') as h5f:
    h5f.create_dataset("E_llp_all", data=np.asarray(E_llp_all))
    h5f.create_dataset("production_vertex_llp_all", data=np.asarray(production_vertex_llp_all))
    h5f.create_dataset("end_vertex_llp_all", data=np.asarray(end_vertex_llp_all))
    h5f.create_dataset("px_llp_all", data=np.asarray(px_llp_all))
    h5f.create_dataset("py_llp_all", data=np.asarray(py_llp_all))
    h5f.create_dataset("pz_llp_all", data=np.asarray(pz_llp_all))
    h5f.create_dataset("m_llp_all", data=np.asarray(m_llp_all))
    h5f.create_dataset("eta_llp_all", data=np.asarray(eta_llp_all))
    h5f.create_dataset("phi_llp_all", data=np.asarray(phi_llp_all))
    h5f.create_dataset("theta_llp_all", data=np.asarray(theta_llp_all))
    h5f.create_dataset("boost_llp_all", data=np.asarray(boost_llp_all))
    h5f.create_dataset("beta_llp_all", data=np.asarray(beta_llp_all))
    h5f.create_dataset("pt_llp_all", data=np.asarray(pt_llp_all))
    h5f.create_dataset("met_llp_all", data=np.asarray(met_llp_all))
    h5f.create_dataset("status_llp_all", data=np.asarray(status_llp_all))
    h5f.create_dataset("n_decay_products_llp_all", data=np.asarray(n_decay_products_llp_all))
    h5f.create_dataset("event_number_llp_all", data=np.asarray(event_number_llp_all))
    h5f.create_dataset("children_summass_llp_all", data=np.asarray(children_summass_llp_all))
    h5f.create_dataset("weight_llp_all", data=np.asarray(weight_llp_all))    

with h5py.File(h5_out_dir+f'/jets_{filedescription}.h5', 'w') as h5f:
    # require list of jet objects, listing each jet's event number, angular information, momentum, energy, and list of constituents' angles
    h5f.create_dataset("phi_jets_all", data=np.asarray(phi_jets_all))
    h5f.create_dataset("event_number_jets_all", data=np.asarray(event_number_jets_all))
    h5f.create_dataset("eta_jets_all", data=np.asarray(eta_jets_all))
#    h5f.create_dataset("theta_jets_all", data=np.asarray(theta_jets_all))
    h5f.create_dataset("pt_jets_all", data=np.asarray(pt_jets_all))
#    h5f.create_dataset("m_jets_all", data=np.asarray(m_jets_all))
    h5f.create_dataset("px_jets_all", data=np.asarray(px_jets_all))
    h5f.create_dataset("py_jets_all", data=np.asarray(py_jets_all))
    h5f.create_dataset("pz_jets_all", data=np.asarray(pz_jets_all))
    h5f.create_dataset("E_jets_all", data=np.asarray(E_jets_all))
    #h5f.create_dataset("production_vertex_jets_all", data=np.asarray())
#    h5f.create_dataset('jet_number_jets_all', data=np.asarray(jet_number_jets_all)  
    h5f.create_dataset("weight_jets_all", data=np.asarray(weight_jets_all))
    
with h5py.File(h5_out_dir+f'/charged_{filedescription}.h5', 'w') as h5f:
    # charged final state prompt particles for isolation requirement 
    h5f.create_dataset("phi_charged_all", data=np.asarray(phi_charged_all))
    h5f.create_dataset("event_number_charged_all", data=np.asarray(event_number_charged_all))
    h5f.create_dataset("eta_charged_all", data=np.asarray(eta_charged_all))
    h5f.create_dataset("theta_charged_all", data=np.asarray(theta_charged_all))
    h5f.create_dataset("pt_charged_all", data=np.asarray(pt_charged_all))
    h5f.create_dataset("m_charged_all", data=np.asarray(m_charged_all))
    h5f.create_dataset("production_vertex_charged_all", data=np.asarray(production_vertex_charged_all))
#    h5f.create_dataset("px_charged_all", data=np.asarray(px_charged_all))
#    h5f.create_dataset("py_charged_all", data=np.asarray(py_charged_all))
#    h5f.create_dataset("pz_charged_all", data=np.asarray(pz_charged_all))
#    h5f.create_dataset("E_charged_all", data=np.asarray(E_charged_all))
    h5f.create_dataset("weight_charged_all", data=np.asarray(weight_charged_all))
