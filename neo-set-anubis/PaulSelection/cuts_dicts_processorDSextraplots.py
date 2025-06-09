import math
import argparse
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import uproot
import awkward as ak
from particle import PDGID, Particle, Charge
import seaborn as sns
import pickle
#from scipy.interpolate import UnivariateSpline
#import scipy
import os
import re
from common_funcs import *
from const_params import *
from Pheno_Eqs import * 
sns.set_style("white")

def extract_inner_keys_values(d):
    result = set()
    
    def recurse(current_dict):
        for key, value in current_dict.items():
            if isinstance(value, dict):
                recurse(value)
            else:
                result.add((key, value))
    
    recurse(d)
    return result

linestyle_str = [
     ('solid', 'solid'),      # Same as (0, ()) or '-'
     ('dotted', 'dotted'),    # Same as (0, (1, 1)) or ':'
     ('dashed', 'dashed'),    # Same as '--'
     ('dashdot', 'dashdot')]  # Same as '-.'

linestyle_tuple = [
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),
     ('long dash with offset', (5, (10, 3))),
     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]

cross_sec_dict = {"CCDY_qqe": 5026., "CCDY_eev": 494., "CCDY_qqv": 1656.,
                  "NCDY_qqe": 3145., "NCDY_eev": 308.8, "NCDY_qqv": 1032.,
                  "Wa_qqe": 3.431, "Wa_eev": 0.3365, "Wa_qqv": 1.125,
                  "ggF_qqe": 0.1992, "ggF_eev": 0.01958, "ggF_qqv": 0.06554,
                  "DS":43000000} # in ab


lumi = 3. / (10 ** -10) # iab, HL-LHC
hbar = 6.582e-25 # GeV.s
n_llp_target_nobkg = 4
n_llp_target_bkg = 90
c_const = 3.e+8
# Loop through mass points for every prod+decay process
# Loop through available run numbers for each mass point

#CCDY_qqv_path = "/r04/atlas/amullin/ANUBIS/SET_ANUBIS/output/4modes/CCDY_qqv"
#NCDY_qqv_path = "/r04/atlas/amullin/ANUBIS/SET_ANUBIS/output/4modes/NCDY_qqv"
#Wa_qqv_path = "/r04/atlas/amullin/ANUBIS/SET_ANUBIS/output/4modes/Wa_qqv"
#ggF_qqv_path = "/r04/atlas/amullin/ANUBIS/SET_ANUBIS/output/4modes/ggF_qqv"

CCDY_eev_path = "/usera/dp728/run_dir/output/SS1_CCDY_eev"
DS_path = "/usera/dp728/run_dir/CutDicts"

#NCDY_eev_path = "/r04/atlas/amullin/ANUBIS/SET_ANUBIS/output/4modes/NCDY_eev"
#Wa_eev_path = "/r04/atlas/amullin/ANUBIS/SET_ANUBIS/output/4modes/Wa_eev"
#ggF_eev_path = "/r04/atlas/amullin/ANUBIS/SET_ANUBIS/output/4modes/ggF_eev"

#CCDY_qqe_path = "/r04/atlas/amullin/ANUBIS/SET_ANUBIS/output/4modes/CCDY_qqe"
#NCDY_qqe_path = "/r04/atlas/amullin/ANUBIS/SET_ANUBIS/output/4modes/NCDY_qqe"
#Wa_qqe_path = "/r04/atlas/amullin/ANUBIS/SET_ANUBIS/output/4modes/Wa_qqe"
#ggF_qqe_path = "/r04/atlas/amullin/ANUBIS/SET_ANUBIS/output/4modes/ggF_qqe"

#masslabels = ['0p5', '0p6', '0p7', '0p8', '0p9', '1', '1p1', '1p2', '1p3', '1p4', '1p5', '1p6', '1p8', '2']
masslabels = ['0pt1', '0pt2','0pt3','0pt4','0pt5','0pt6','0pt7','0pt8','0pt9','1pt0','1pt2','1pt4', '1pt6',]
masslabels = ['0pt1', '0pt3','0pt5','0pt6','0pt7','0pt9','1pt1','1pt3','1pt5','1pt7','1pt9','2pt1','2pt3']
masslabels = ['0pt5','0pt6','0pt7','0pt9','1pt1','1pt3','1pt5','1pt7','1pt9','2pt1','2pt3']
masslabels = ['10p0','20p0','30p0','40p0','50p0','60p0']

# Speed of light in meters per second
c = 3e8

# Create log-spaced array from 0.1 to 10^8
arr = np.logspace(-1, 8, num=10)

# Divide each value by the speed of light
lifetimes = arr / c
print(lifetimes)
lifetimelabels = lifetimes

#processnames = ["CCDY_qqe", "CCDY_eev", "CCDY_qqv", 'Wa_qqe', 'ggF_qqe', 'ggF_qqv', 'ggF_eev','NCDY_qqe', 'NCDY_eev', 'NCDY_qqv']
#processnames = ["NCDY_eev", "NCDY_qqe", "NCDY_qqv", "Wa_eev", "Wa_qqe", "Wa_qqv", "CCDY_eev", "CCDY_qqe", "CCDY_qqv", "ggF_eev", "ggF_qqe", "ggF_qqv"]
processnames = ["CCDY_eev"]
processnames = {"DS"}
#processnames = ["NCDY_qqv"]
#processnames_incl = ["NCDY", "CCDY", "Wa", "ggF"] #["CCDY", "ggF", "Wa", 'NCDY']
processnames_incl = ["CCDY"]
processnames_incl = {"DS"}
#processnames_incl = ["NCDY"]

big_cuts_dict = {}
big_cuts_dict_new = {}

CCDY_qqe_br = {}
CCDY_qqv_br = {}
CCDY_eev_br = {}

NCDY_qqe_br = {}
NCDY_qqv_br = {}
NCDY_eev_br = {}

ggF_qqe_br = {}
ggF_qqv_br = {}
ggF_eev_br = {}

Wa_qqe_br = {}
Wa_qqv_br = {}
Wa_eev_br = {}

CCDY_incl_br = {}
NCDY_incl_br = {}
NCDY_incl_br = {}
Wa_incl_br = {}
import glob

#Find all the matching files
for processname in processnames:
    for masslabel in masslabels:
        for lifetimelabel in lifetimelabels:
            lifetimelabel = str(lifetimelabel).replace('.', 'p')
            pattern1 = re.compile(rf'{lifetimelabel}.*{masslabel}.*\.pickle$')
            
            print("pattern1 is: ", pattern1)
            #all_files = os.listdir(eval(f"{processname}_path")) # path to the directory for this process 
            base_path = eval(f"{processname}_path")
            all_files = glob.glob(f"{base_path}/**/*", recursive=True)

            print(f"{processname}_path")
            print(masslabel)
            print(lifetimelabel)
            # read in all files for CCDY_qqe that have this masslabel (i.e. all different runs) 
            
            matching_files1 = [f for f in all_files if re.search(pattern1,f)]
        
            #print(matching_files1)

            for file_path in matching_files1:
                with open(file_path, 'rb') as file:
                    # Load the pickle file for cuts dict
                    cuts_dict = pickle.load(file)
                    if isinstance(cuts_dict, dict): # ensuring data is a dictionary
                        file_name = os.path.basename(file_path)
                        big_cuts_dict_new[file_name] = cuts_dict
                    #print(cuts_dict)
                

#print("big_cuts_dict_new: ")
#print(big_cuts_dict_new)

            
#        big_cuts_dict[f"{processname}_{masslabel}_run{i}"] = proc_cuts

#print("big_cuts_dict: ")
#print(big_cuts_dict)


# stores info as dict of {filename: {cuts_dict}, ...} for all files (i.e. all process+mass combinations) where cuts_dict = {cut: Nevents, ... }
# Now access cuts for plotting

# need 16 different sets (dict or list)  of BR values (listed over the mass range) (12 normal plus 4 inclusive which are just the 3 decays summed I guess

br_sens_dict = dict.fromkeys(processnames)
br_sens_dict_new = dict.fromkeys(processnames)
br_sens_dict_incl = dict.fromkeys(processnames_incl)

br_sens_dict_uncertainty = dict.fromkeys(processnames)
br_sens_dict_new_uncertainty = dict.fromkeys(processnames)

# loop through all found files and calculate the BR sensitivity for each process+mass+decay combination
for processname in processnames:
    br_sens_dict[processname] = dict.fromkeys(masslabels)
    
    br_sens_dict_uncertainty[processname] = dict.fromkeys(masslabels)
    for masslabel in masslabels:
        br_sens_dict[processname][masslabel] = dict.fromkeys(lifetimelabels)
        br_sens_dict_uncertainty[processname][masslabel] = dict.fromkeys(lifetimelabels)
        for lifetimelabel in lifetimelabels:

            # cuts_dict_HeavyN_Maj_2000evt_CCDY_qqv_mn1_1_madspin_run1.pkl 
            print("doing ", processname, " for mass ", masslabel, " and lifetime ", lifetimelabel)
            lifetimelabel = str(lifetimelabel).replace('.', 'p')
            pattern = rf'.*_{processname}_.*_{lifetimelabel}_.*_{masslabel}_.*'

            regex = re.compile(pattern)
            matching_keys = [key for key in big_cuts_dict if regex.search(key)] # this should store all the keys (i.e. runs) relevant to this mass for this process
            #print(matching_keys)
            Ntot_list = []
            Nobs_list = []
            for thiskey in matching_keys:
                #print(big_cuts_dict[thiskey])
                Ntot_list.append(big_cuts_dict[thiskey]['All'])
                #final cut equivalent to all of em 
                Nobs_list.append(big_cuts_dict[thiskey]['nLLP_IsoAll'])
                if (big_cuts_dict[thiskey]['nLLP_IsoAll'] == 0):
                    print("NO EVENTS OBS ", processname, " at mass ", masslabel, " and lifetime ", lifetimelabel)
                    #print(big_cuts_dict[thiskey])
            #print(Ntot_list)
            #print(Nobs_list)
            Ntot = sum(Ntot_list)
            Nobs = sum(Nobs_list)
    #        print(Ntot)
    #        print(Nobs)
            # Now calculate BR sensitivity 
            if Nobs != 0:
                print("an event")
                br_sens_dict[processname][masslabel][lifetimelabel] = n_llp_target_bkg / ( cross_sec_dict[processname] * lumi ) * (Ntot / Nobs)# doubled from 2 LLP but halfed cuz of phi folding
                br_sens_dict_uncertainty[processname][masslabel][lifetimelabel] = br_sens_dict[processname][masslabel][lifetimelabel] *np.sqrt((np.sqrt(Nobs)/Nobs)^2+np.sqrt((Ntot)/Ntot)^2) 
            else:
                br_sens_dict[processname][masslabel][lifetimelabel] = 1. # is 0 but change to one as cant resolve any branching ratio
                print("no events observed in ", processname, " for mass ", masslabel, " and lifetime ", lifetimelabel)
            Plot = False
  



for processname in processnames:
    br_sens_dict_new[processname] = dict.fromkeys(masslabels)

    br_sens_dict_new_uncertainty[processname] = dict.fromkeys(masslabels)

    for masslabel in masslabels:
        br_sens_dict_new[processname][masslabel] = dict.fromkeys(lifetimelabels)
        br_sens_dict_new_uncertainty[processname][masslabel] = dict.fromkeys(lifetimelabels)
        for lifetimelabel in lifetimelabels:
            # cuts_dict_HeavyN_Maj_2000evt_CCDY_qqv_mn1_1_madspin_run1.pkl 
    #        print("doing ", processname, " for mass ", masslabel)
            lifetimelabel = str(lifetimelabel).replace('.', 'p')
            #pattern = rf'.*_{lifetimelabel}_.*_{masslabel}_.*'
            pattern = rf'.*{lifetimelabel}.*{masslabel}.*'

            print("lifetimelabel is: ", lifetimelabel)
            print("masslabel is: ", masslabel)
            regex = re.compile(pattern)
            matching_keys = [key for key in big_cuts_dict_new if regex.search(key)] # this should store all the keys (i.e. runs) relevant to this mass for this process
            #print("keys in big_cuts_dict_new: ", big_cuts_dict_new.keys())
            #print(matching_keys)
            Ntot_list = []
            Nobs_list = []
            for thiskey in matching_keys:
                #print("keys in big_cuts_dict_new[thiskey]: ", big_cuts_dict_new[thiskey])
                Ntot_list.append(big_cuts_dict_new[thiskey]['cutFlow']['nLLP_original'])
                #final cut equivalent to all of em 
                Nobs_list.append(big_cuts_dict_new[thiskey]['cutFlow']['nLLP_IsoAll']) #U didnt have this as NEW
                if (big_cuts_dict_new[thiskey]['cutFlow']['nLLP_IsoAll'] == 0):
                    print("NO EVENTS OBS ", processname, " at mass ", masslabel, " and lifetime ", lifetimelabel)
                    #print(big_cuts_dict_new[thiskey])
            #print(Ntot_list)
            #print(Nobs_list)
            Ntot = sum(Ntot_list)
            Nobs = sum(Nobs_list)
            print(Ntot)
            print(Nobs)
            # Now calculate BR sensitivity 
            if Nobs != 0:
                #2 for higgs
                br_sens_dict_new[processname][masslabel][lifetimelabel] = n_llp_target_bkg / (2* cross_sec_dict[processname] * lumi ) * (Ntot / Nobs)
                br_sens_dict_new_uncertainty[processname][masslabel][lifetimelabel] = br_sens_dict[processname][masslabel][lifetimelabel] *(np.sqrt(Nobs)/Nobs)

            else:
                br_sens_dict_new[processname][masslabel][lifetimelabel] = 1. # is 0 but change to one as cant resolve any branching ratio
                print("no events observed in ", processname, " for mass ", masslabel)
#print("br_sens_dict: ")
#print(br_sens_dict)
#print("br_sens_dict_new: ")
#print(br_sens_dict_new)
                
#not built for inclusive atm                
"""
for processname_incl in processnames_incl:
    br_sens_dict_incl[processname_incl] = dict.fromkeys(masslabels)
    for masslabel in masslabels:
        br_sens_dict_incl[processname_incl][masslabel] = dict.fromkeys(lifetimelabels)
        for lifetimelabel in lifetimelabels:
            br_sens_dict_incl[processname_incl][masslabel][lifetimelabel] = 0
    #        print("doing ", processname_incl, " for mass ", masslabel)
            # Require a list of summed br sens values for each 4 prod modes
            pattern_incl = rf'.*{processname_incl}.*'
            regex_incl = re.compile(pattern_incl)
            matching_keys_incl = [key for key in br_sens_dict if regex_incl.search(key)]
    #        print("matching_keys_incl: ", matching_keys_incl)
            for thiskey in matching_keys_incl:
    #            print(br_sens_dict[thiskey][masslabel])
                br_sens_dict_incl[processname_incl][masslabel][lifetimelabel] += br_sens_dict[thiskey][masslabel][lifetimelabel]
    #            print("adding ", br_sens_dict[thiskey][masslabel], " to process ", processname_incl, " for mass ", masslabel)
"""

#print(br_sens_dict_incl)        
#        CCDY_qqe_Ntot = big_cuts_dict 
#        CCDY_qqe_Nobs = 
#        CCDY_qqe_br[str(masslabel)] = n_llp_target_bkg / ( cross_sec_dict['CCDY_qqe'] * lumi ) * 

# Now make plots from the big_cuts_dict
# Select cuts values to add together by selecting keys (file names) which contain the same process 

#br_plot = plt.figure()

plotstyle_dict = {"CCDY_qqe": {'line':'dotted','colour':'orange'}, "CCDY_eev": {'line':'solid','colour':'orange'}, "CCDY_qqv": {'line':'dashed', 'colour':'orange'},
                  "NCDY_qqe": {'line':'dotted','colour':'blue'}, "NCDY_eev": {'line':'solid','colour':'blue'}, "NCDY_qqv": {'line':'dashed', 'colour':'blue'},
                  "Wa_qqe": {'line':'dotted','colour':'green'}, "Wa_eev": {'line':'solid','colour':'green'}, "Wa_qqv": {'line':'dashed', 'colour':'green'},
                  "ggF_qqe": {'line':'dotted','colour':'red'}, "ggF_eev": {'line':'solid','colour':'red'}, "ggF_qqv": {'line':'dashed', 'colour':'red'},
                  #Added one for DS
                  "DS": {'line':'dotted','colour':'blue'}
                  }


#print("!!!!!!!! LINESTYLE: linestyle_tuple[0[1]] = ", linestyle_tuple[0[1]])

#for processname in processnames:
#    masstoplot_nan = list(br_sens_dict[processname].keys())
#    masstoplot = [float(item.replace('p','.')) for item in masstoplot_nan]
##    print(masstoplot)
#    brtoplot = list(br_sens_dict[processname].values())
#    plt.plot(masstoplot, brtoplot, linestyle='-',label=processname)

#plt.title("BR sensitivity")
#plt.legend()
#plt.xlabel("Mass [Gev]")
#plt.ylabel("BR")
#plt.yscale('log')

br_plot_new, ax = plt.subplots()
for masslabel in masslabels:
    plt.clf()
    for processname in processnames:
        lifetimetoplot_nan = list(br_sens_dict_new[processname][masslabel].keys())
        #lifetimetoplot = [float(item.replace('pt','.')) for item in lifetimetoplot_nan]
        lifetimetoplot = [float(str(item).replace('p','.')) for item in lifetimetoplot_nan]

    #    print(masstoplot)
        brtoplot = list(br_sens_dict_new[processname][masslabel].values())
        brtoplotUncertainty = list(br_sens_dict_new_uncertainty[processname][masslabel].values())
        #incase of Nones
        brtoplotUncertainty_clean = [np.nan if val is None else val for val in brtoplotUncertainty]

        plt.plot(c * np.array(lifetimetoplot), brtoplot, linestyle=plotstyle_dict[processname]['line'],label=processname,color=plotstyle_dict[processname]['colour'])
        #plt.errorbar(c * np.array(lifetimetoplot), brtoplot, yerr = brtoplotUncertainty_clean, color=plotstyle_dict[processname]['colour'])


    plt.title("BR sensitivity")
    plt.legend()
    plt.xlabel(r"$c\tau$ [m]",fontsize=10)#check as there mm scaling at one point
    plt.ylabel("BR",fontsize=10)
    plt.yscale('log')
    plt.xscale('log')

    plt.grid(True, linestyle='--')
    ax.minorticks_on()
    #ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
    plt.savefig(f"BR sensitivity for Mass {masslabel} GeV ")

# convert p in masslabels into . --> plot as numbers
 #   for masslabel in masslabels:
#        print(br_sens_dict[processname][masslabel])
#        masstoplot_nan = (br_sens_dict[processname][masslabel]).replace('p', '.')
#        masstoplot.append(float(masstoplot_nan))
        #        for cut in range(masslabels):
#        brtoplot.append(br_sens_dict[processname][masslabel][cut])
#        print(br_sens_dict[processname][masslabel][cut])
        
    #plt.plot(br_sens_dict[processname][

br_plot_new_lifetime, ax = plt.subplots()

suffix = ""
makePlots = True
verbose = True
import plotting as myplot
import datetime
outputDir = "/usera/dp728/run_dir/DSplots"

cutDict = "big_cuts_dict_new"
print(br_sens_dict_new)


"""
if makePlots:
    if verbose:
        plotStartTime = datetime.datetime.now()
        print(f"Plotting... ({plotStartTime})")

    # Plot the Number of events that pass each stage of the cutFlow
    myplot.plotCutFlow(cutDict["cutFlow"], logScale = True, outputDir=outputDir, suffix=f"{suffix}")
    
    # Plot each of the columns of the dataframe for each item in the sampleDFs
    myplot.plotStandardPlots(sampleDFs, outDir=outputDir, suffix=f"{suffix}")

    # MET derived by different methods -- debugging tool. 
    #   - NOTE: This could be commented out
    #myplot.plotMETcomparison(sampleDFs, outputDir=outputDir, suffix=f"{suffix}")

    # Plot the LLP decay positions compared to the ATLAS/ANUBIS geometry at different stages of selection
    #   - Split off the Geometry positions into a subfolder
    geoOutputDir = f"{outputDir}/geoPlots/"
    if not os.path.exists(geoOutputDir):
        os.makedirs(geoOutputDir)

    myplot.plotPositionsDuringSelection(sampleDFs["LLPs"], cutDict, selection, nBins=[100,100,100], decayVertex="decayVertex_weighted", 
                                        outDir=geoOutputDir, suffix=suffix)

"""



"""
for processname in processnames:
    masstoplot_nan = list(br_sens_dict_new[processname].keys())
    masstoplot = [float(item.replace('pt','.')) for item in masstoplot_nan]
    brtoplot = list(br_sens_dict_new[processname].values())
    # Now using function dec_tot(M,Ve,Vmu,Vtau) to calculate the lifetime for each mass --> swap mass range for lifetimes
    #lifetimetoplot = [c_const*(hbar/dec_tot(mass,1,0.,0.)) for mass in masstoplot] # units of metres (ctau) 
    lifetimetoplot = lifetimes
    plt.plot(lifetimetoplot, brtoplot, linestyle='-',label=processname)
                      
plt.title("BR sensitivity (new, random seed)")
plt.legend()
plt.xlabel("Lifetime ctau [m]")
plt.ylabel("BR")
plt.yscale('log')
plt.xscale('log')
plt.grid(True, linestyle='--')
ax.minorticks_on()
plt.savefig("BR sensitivity (random seed)")


ccdy_eev_wid_tot_mg = 5.543081e-13
ccdy_eev_lifetime_mg = c_const*(hbar/ccdy_eev_wid_tot_mg)
ccdy_eev_lifetime_sofie = c_const*(hbar/dec_tot(1, 1, 0, 0))
print("CHECK LIFETIME: CCDY_eev mn1=1GeV")
print("mg lifetime = ", ccdy_eev_lifetime_mg)
print("sofie calc lifetime = ", ccdy_eev_lifetime_sofie)
masstoplotbig = [mass -0.4 for mass in masstoplot]
print("masses: ")
print(masstoplotbig)
lifetimes_converted = [c_const*(hbar/dec_tot(mass, 1., 0., 0.)) for mass in masstoplotbig]

# print list of lifetimes converted from masses:
print("lifetimes ") 
print(lifetimes_converted)

    
br_plot_incl, ax = plt.subplots()

for processname_incl in processnames_incl:
    masstoplot_nan = list(br_sens_dict_incl[processname_incl].keys())
    masstoplot = [float(item.replace('pt','.')) for item in masstoplot_nan]
#    print(masstoplot)
    brtoplot = list(br_sens_dict_incl[processname_incl].values())
    plt.plot(masstoplot, brtoplot, linestyle='-',label=processname_incl)

plt.title("BR sensitivity, inclusive")
plt.legend()
plt.xlabel("Mass [Gev]",fontsize=10)
plt.ylabel("BR", fontsize=10)
plt.yscale('log')
#ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10)
plt.savefig("BR sensitivity, inclusive")
plt.show()
"""