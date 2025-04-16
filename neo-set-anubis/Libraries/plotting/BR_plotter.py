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
import ROOT
#from root_numpy import root2array, tree2array
import uproot
import awkward as ak
from particle import PDGID, Particle, Charge
from Common_LLP_processing import *
import h5py
import seaborn as sns
import pickle
from colour import Color
from scipy.interpolate import UnivariateSpline
import scipy
sns.set_style("white")

# default = plot only one production+decay process (i.e. only one set of samples) 

# this script is for only making the BR plot(s) (not cutflows)

# Need to store xsec, neventspassed, neventstot, br

masslabels = ["0p5","0p6","0p7","0p8","0p9","1","1p1","1p2","1p3","1p4","1p5"]
massnumbers = [0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]

#elif decaymode == "qqell":
#    masslabels = ["0p5","0p6","0p7","0p8","0p9","1","1p1","1p2","1p3","1p4","1p5"]
#    massnumbers = [0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5]

br_product_w_elljj = []
br_product_w_ellellvv = []
br_product_w_qqell = []
br_product_w_qqvv = []
br_product_w_combined = []
br_product_z_combined = []
br_product_z_combined_2 = []
br_product_z_averaged = []
br_product_w_qqell_tot = []
br_product_CCDY_incl = []

num_colours = len(masslabels)
red = Color("red")
green = Color("green")
colourlist = list(red.range_to(Color("yellow"),num_colours))
colourlist2 = list(green.range_to(Color("blue"),num_colours))
## Define command line arguments
parser = argparse.ArgumentParser(description='Make cutflows.')
parser.add_argument('-pdf', default='sensitivity_plots.pdf',help='PDF file with cutflow plots.')
args = parser.parse_args()

# File for saving plots:
sensitivity_plots = PdfPages(args.pdf)

# For cutflows: need 2 lists: cutnames (label for x axis of each cut) and passedcuts (number of events passing each cut cumulatively). Can do this from pickled dictionary.

figsens, ax = plt.subplots(1, 1, figsize=(8, 6))

lumi = 3 / (10 ** -10) #iab --> HL-LHC
n_llp_target = 4

xsec_w_hnl_elljj = 4600. # pb # probably should update this with production and decay 
xsec_w_hnl_ellellvv = 460. # pb
xsec_w_hnl_qqell = 4640. # pb
xsec_w_hnl_qqvv = 1513. # pb
xsec_w_hnl_combined = 6570. # pb
xsec_z_hnl_combined = 4040. # pb
xsec_CCDY_incl = 6667. # pb

# xsec, n_hnl, nevents_input, point_label, 
# Make these into dictionaries instead of lists, so that the code is more robust 
w_elljj = [xsec_w_hnl_elljj, xsec_w_hnl_elljj * lumi, "2000", "pp_w_n1ell_elljj"]
w_ellellvv = [xsec_w_hnl_ellellvv, xsec_w_hnl_ellellvv * lumi, "2000", "pp_w_n1ell_ellellvv" ]
w_qqell = [xsec_w_hnl_qqell, xsec_w_hnl_qqell * lumi, "2000", "pp_w_n1ell_qqell"]
w_qqvv = [xsec_w_hnl_qqvv, xsec_w_hnl_qqvv * lumi, "2000", "pp_w_n1ell_qqvv"]
w_combined = [xsec_w_hnl_combined, xsec_w_hnl_combined * lumi, "5000", "pp_w_n1ell_combined"]
z_combined = [xsec_z_hnl_combined, xsec_z_hnl_combined * lumi, "1000", "pp_z_combined"]
z_combined_2 = [xsec_z_hnl_combined, xsec_z_hnl_combined * lumi, "500", "pp_z_combined"]
CCDY_incl = [xsec_CCDY_incl, xsec_CCDY_incl * lumi, "1000", "CCDY_incl"]

mass_to_plot = []
massindex = 0

for masslabel in masslabels:
    masspoint = massnumbers[massindex]
    mass_to_plot.append(masspoint)
    print("computing mass ", masspoint)
    filedescription_w_elljj = "cuts_dict_" + str(w_elljj[3]) + f"_nm1_{masslabel}_wn1_auto_qcdnonzero_" + str(w_elljj[2]) + "evt.pkl"
    print(filedescription_w_elljj)

    filedescription_w_ellellvv = "cuts_dict_" + str(w_ellellvv[3]) + f"_nm1_{masslabel}_wn1_auto_qcdnonzero_" + str(w_ellellvv[2]) + "evt.pkl"
    print(filedescription_w_ellellvv)
    
    filedescription_w_qqell = "cuts_dict_" + str(w_qqell[3]) + f"_nm1_{masslabel}_wn1_auto_qcdnonzero_" + str(w_qqell[2]) + "evt.pkl"
    print(filedescription_w_qqell)

    filedescription_w_qqvv = "cuts_dict_" + str(w_qqvv[3]) + f"_nm1_{masslabel}_wn1_auto_qcdnonzero_" + str(w_qqvv[2]) + "evt.pkl"
    print(filedescription_w_qqvv)

    filedescription_w_combined = "cuts_dict_" + str(w_combined[3]) + f"_nm1_{masslabel}_wn1_auto_qcdnonzero_" + str(w_combined[2]) + "evt.pkl"
    print(filedescription_w_combined)

    filedescription_z_combined = "cuts_dict_" + str(z_combined[3]) + f"_nm1_{masslabel}_wn1_auto_qcdnonzero_" + str(z_combined[2]) + "evt.pkl"
    print(filedescription_z_combined)

    filedescription_z_combined_2 = "cuts_dict_" + str(z_combined_2[3]) + f"_nm1_{masslabel}_wn1_auto_qcdnonzero_" + str(z_combined_2[2]) + "evt.pkl"
    print(filedescription_z_combined_2)

    filedescription_CCDY_incl = "cuts_dict_HeavyN_Maj_" + str(CCDY_incl[2]) + "evt_" + str(CCDY_incl[3]) + f"_mn1_{masslabel}.pkl"
    print(filedescription_CCDY_incl)

    with open(f'{filedescription_w_elljj}', 'rb') as fp_w_elljj, open(f'{filedescription_w_ellellvv}', 'rb') as fp_w_ellellvv, open(f'{filedescription_w_qqell}', 'rb') as fp_w_qqell, open(f'{filedescription_w_qqvv}', 'rb') as fp_w_qqvv, open(f'{filedescription_w_combined}', 'rb') as fp_w_combined, open(f'{filedescription_z_combined}', 'rb') as fp_z_combined, open(f'{filedescription_z_combined}', 'rb') as fp_z_combined_2, open(f'{filedescription_CCDY_incl}', 'rb') as fp_CCDY_incl:
        cuts_dict_w_elljj = pickle.load(fp_w_elljj)
        cuts_dict_w_ellellvv = pickle.load(fp_w_ellellvv)
        cuts_dict_w_qqell = pickle.load(fp_w_qqell)
        cuts_dict_w_qqvv = pickle.load(fp_w_qqvv)
        cuts_dict_w_combined = pickle.load(fp_w_combined)
        cuts_dict_z_combined = pickle.load(fp_z_combined)
        cuts_dict_z_combined_2 = pickle.load(fp_z_combined_2)
        cuts_dict_CCDY_incl = pickle.load(fp_CCDY_incl)

        keys_w_elljj, values_w_elljj = zip(*cuts_dict_w_elljj.items())
        keys_w_ellellvv, values_w_ellellvv = zip(*cuts_dict_w_ellellvv.items())
        keys_w_qqell, values_w_qqell = zip(*cuts_dict_w_qqell.items())
        keys_w_qqvv, values_w_qqvv = zip(*cuts_dict_w_qqvv.items())
        keys_w_combined, values_w_combined = zip(*cuts_dict_w_combined.items())
        keys_z_combined, values_z_combined = zip(*cuts_dict_z_combined.items())
        keys_z_combined_2, values_z_combined_2 = zip(*cuts_dict_z_combined_2.items())
        keys_CCDY_incl, values_CCDY_incl = zip(*cuts_dict_CCDY_incl.items())

        br_product_w_elljj.append( (n_llp_target/w_elljj[1]) * (values_w_elljj[0]/values_w_elljj[4]) )
        br_product_w_ellellvv.append( (n_llp_target/w_ellellvv[1] * (values_w_ellellvv[0]/values_w_ellellvv[4]) )) 
        br_product_w_qqell.append( (n_llp_target/w_qqell[1]) * (values_w_qqell[0]/values_w_qqell[4]) ) 
        br_product_w_qqvv.append( (n_llp_target/w_qqvv[1]) * (values_w_qqvv[0]/values_w_qqvv[4]) )
        br_product_w_qqell_tot.append( 0.5 * ( ( (n_llp_target/w_elljj[1]) * (values_w_elljj[0]/values_w_elljj[4]) ) + ( (n_llp_target/w_qqell[1]) * (values_w_qqell[0]/values_w_qqell[4])) ) )
        br_product_w_combined.append( (n_llp_target/w_combined[1]) * (values_w_combined[0]/values_w_combined[4]) )
        br_product_z_combined.append( (n_llp_target/z_combined[1]) * (values_z_combined[0]/values_z_combined[4]) )
        br_product_z_combined_2.append( (n_llp_target/z_combined_2[1]) * (values_z_combined_2[0]/values_z_combined_2[4]) )
        br_product_z_averaged.append( 0.5 * ( ( (n_llp_target/z_combined[1]) * (values_z_combined[0]/values_z_combined[4]) ) +  ( (n_llp_target/z_combined_2[1]) * (values_z_combined_2[0]/values_z_combined_2[4]) )) )
        br_product_CCDY_incl.append( 0.5 * ( ( (n_llp_target/CCDY_incl[1]) * (values_CCDY_incl[0]/values_CCDY_incl[4]) ) +  ( (n_llp_target/CCDY_incl[1]) * (values_CCDY_incl[0]/values_CCDY_incl[4]) )) )
    massindex+=1
    #m0p1, color="dodgerblue", label="mass 0.1 GeV", **kwargs)
    #m0p5, color="orange", label="mass 0.5 GeV", **kwargs)
    #m1, color="deeppink", label="mass 1 GeV", **kwargs)
    #m2p5, color="thistle", label="mass 2.5 GeV", **kwargs)
    #m5, color="palegoldenrod", label="mass 5 GeV", **kwargs)
    #m7p5, color="palegreen", label="mass 7.5 GeV", **kwargs)
    #m10, color="seagreen", label="mass 10 GeV", **kwargs)
#br_product_z_averaged[1] = 1.4e-11
#br_product_z_averaged[2] = 7.e-12
#br_product_z_averaged[5] = 1.17e-12
#print(br_product_z_averaged)
    
br_plot = plt.figure()
#plt.plot(mass_to_plot, br_product_w_elljj, linestyle='-', color='deeppink',label='W production, final state = jj + electron')
plt.plot(mass_to_plot, br_product_w_ellellvv, linestyle='-', color='thistle',label='W production, final state = electrons + neutrino')
#plt.plot(mass_to_plot, br_product_w_qqell, linestyle='-', color='blue',label='W production, final state = jj + electron')
plt.plot(mass_to_plot, br_product_w_qqell_tot, linestyle='-', color='seagreen',label='W production, final state = jets + electron')
plt.plot(mass_to_plot, br_product_w_qqvv, linestyle='-', color='palegreen',label='W production, final state = jets + neutrino')
plt.plot(mass_to_plot, br_product_w_combined, linestyle='-', color='deeppink',label='W production (inclusive)')
#plt.plot(mass_to_plot, br_product_z_combined, linestyle='-', color='black',label='Z combined')
#plt.plot(mass_to_plot, br_product_z_combined_2, linestyle='-', color='red',label='Z combined 2')
plt.plot(mass_to_plot, br_product_z_averaged, linestyle='-', color='orange',label='Z production (inclusive)')
plt.plot(mass_to_plot, br_product_CCDY_incl, linestyle='-', color='red',label='CCDY incl')
plt.title("Branching ratio sensitivity")
plt.ylabel("BR")
plt.xlabel("LLP mass [GeV]")
#plt.xscale("log")
plt.legend()
plt.yscale("log")
sensitivity_plots.savefig(br_plot,bbox_inches="tight")

smoothbr = plt.figure()
# need to calculate the coupling dependence of the cross-section (and BR)
spl = UnivariateSpline(mass_to_plot,br_product_w_elljj)
xspl = np.linspace(min(mass_to_plot),max(mass_to_plot),1000)
spl.set_smoothing_factor(0.2)
plt.plot(xspl,spl(xspl),'b',lw=3,label="W elljj")

#w_elljj_spline = scipy.interpolate.make_interp_spline(mass_to_plot, br_product_w_elljj)
#plt.plot(mass_to_plot,w_elljj_spline,'y',lw=3,label="W elljj")

spl2 = UnivariateSpline(mass_to_plot,br_product_w_ellellvv)
spl2.set_smoothing_factor(0.2)
plt.plot(xspl,spl2(xspl),'g',lw=3,label="W ellellvv")

spl3 = UnivariateSpline(mass_to_plot,br_product_w_combined)
spl3.set_smoothing_factor(0.2)
plt.plot(xspl,spl3(xspl),'r',lw=3, label="W combined")
plt.yscale("log")
#spl = UnivariateSpline(mass_to_plot,br_product4)
#plt.plot(xspl,spl(xspl),'b',lw=3)
plt.legend()

#import numpy as np
#from scipy.interpolate import interp1d
#import matplotlib.pyplot as plt

# Assuming your data is in the variables 'x' and 'y'
#x = np.array(x)
#y = np.array(y)

# Create the spline interpolation
#f = interp1d(x, y, kind='cubic')

# Create a new, smoother x-axis
#x_new = np.linspace(x.min(), x.max(), 100)

# Evaluate the spline on the new x-axis
#y_new = f(x_new)

# Plot the original and smoothed data
#plt.figure(figsize=(10, 6))
#plt.plot(x, y, 'bo-', label='Original Data')
#plt.plot(x_new, y_new, 'r-', label='Spline Interpolation')
#plt.xlabel('x')
#plt.ylabel('y')
#plt.legend()
#plt.show()

plt.show()

sensitivity_plots.savefig(smoothbr,bbox_inches="tight")



    #m0p1, color="dodgerblue", label="mass 0.1 GeV", **kwargs)
    #m0p5, color="orange", label="mass 0.5 GeV", **kwargs)
    #m1, color="deeppink", label="mass 1 GeV", **kwargs)
    #m2p5, color="thistle", label="mass 2.5 GeV", **kwargs)
    #m5, color="palegoldenrod", label="mass 5 GeV", **kwargs)
    #m7p5, color="palegreen", label="mass 7.5 GeV", **kwargs)
    #m10, color="seagreen", label="mass 10 GeV", **kwargs)


# make VeN1 scaling script and calculate sensitivity (without met cut for now, so that numbers are non-zero) across 2D scan in mass vs VeN1. Need to calculate the minimum NeN1 at which we have sensitivity for each mass point. Could make a formula which spits out VeN1 corresponding to the scale factor that's required to be multiplied to the number of output events at VeN1=1 to reach the necessary number of output events for detection. First need to find what that necessary #events is. E.g. if no backgrounds then it's 4 events? (Or maybe 90 if backgrounds..? Try 4 for now)

# Scale by VeN1:
# check PlotBRs.ipynb



# calculate branching ratio sensitivity


# check plot_br function


# compute yield
# number of accepted events


    
cutflow_plots.close()
print("Plots output to %s." % args.pdf)
