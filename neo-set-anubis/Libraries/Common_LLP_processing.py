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

# Convert r, eta, phi to Cartesian coordinates          
def to_cartesian(r, eta, phi):
    x = r * np.sin(to_theta(eta)) * np.cos(phi)
    y = r * np.sin(to_theta(eta)) * np.sin(phi)
    z = r * np.cos(to_theta(eta))
    return [x, y, z]

def Average_list(lst):
    return sum(lst) / len(lst)

def calculate_deltaR(deltaeta, deltaphi): ####### What type of delta R are we cutting on?? dRClosestJet is already in the dataframe in Toby's code
    # dR between which jets? 
    deltaR = np.sqrt(np.square(deltaeta) + np.square(deltaphi))
    return deltaR

# Calculate boost (gamma)
def calculate_boost(pt, eta, m):
    p = pt * np.cosh(eta)
    e = np.sqrt(p * p + m * m)
    return e / m

# Calculate beta
def calculate_beta(gamma):
    return np.sqrt(1 - (1 / np.square(gamma)))

# Convert theta to eta
def to_eta(theta):
    return -1 * np.log(np.tan(theta / 2))

#calculate phi (azimuthal angle, i.e. angle from x axis)
def calculate_phi(px, py):
#    return np.arctan(py/px) #### probably wrong -- update 
    if px > 0:
        phi = np.arctan(py / px)
    elif px < 0 and py >= 0:
        phi = np.arctan(py / px) + np.pi
    elif px < 0 and py < 0:
        phi = np.arctan(py / px) - np.pi
    elif px == 0 and py > 0:
        phi = np.pi / 2
    elif px == 0 and py < 0:
        phi = -np.pi / 2
    else:
        phi = np.nan
    return phi

# Convert eta to theta
def to_theta(eta):
    return 2 * np.arctan(np.exp(-eta))

def calculate_eta(px, py, pz):
    if ( ((1/2)*np.log( (np.sqrt(np.square(px) + np.square(py) + np.square(pz)) + pz)/(np.sqrt(np.square(px) + np.square(py) + np.square(pz)) - pz) )) > 1000000):
         return 1000000
    elif ( ((1/2)*np.log( (np.sqrt(np.square(px) + np.square(py) + np.square(pz)) + pz)/(np.sqrt(np.square(px) + np.square(py) + np.square(pz)) - pz) )) < -1000000):
         return -1000000
    else:
         return (1/2)*np.log( (np.sqrt(np.square(px) + np.square(py) + np.square(pz)) + pz)/(np.sqrt(np.square(px) + np.square(py) + np.square(pz)) - pz) )

def calculate_pt(px, py):
    return np.sqrt(np.square(px) + np.square(py))

# Exponential decay of particle with certain boost and ctau
def expo_decay(boost, beta, ctau):
    return np.random.exponential(scale=(boost * beta * ctau))


# Exponential decay of particle with certain boost and ctau
def decay_llp_to_r(m, pt, eta, phi, ctau):
    # Calculate boost from pt, eta, mass:
    boost = calculate_boost(pt, eta, m)
    # Calculate beta from boost: 
    beta = calculate_beta(boost)
    # Calculate decay position in spherical coords: 
#    r = np.random.exponential(scale=(boost * beta * ctau))
    r = boost*beta*ctau # in metres (depending on ctau units...) 
    # Convert to cartesian decay position: 
#    decay_pos_cartesian = to_cartesian(r, eta, phi)
    return r #decay_pos_cartesian 

    
