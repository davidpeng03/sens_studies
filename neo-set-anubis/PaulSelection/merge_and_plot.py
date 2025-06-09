# merge_and_plot.py
import pandas as pd
import glob
import matplotlib.pyplot as plt
import os
import helper as h
import defineGeometry as g
import numpy as np

# Load all output files
output_dir = "/usera/dp728/run_dir/output/HNL"
# 1. Merge all filtered_df CSVs
filtered_df_files = glob.glob(f"{output_dir}/filtered_df_*.csv")
print(filtered_df_files)
filtered_df_list = [pd.read_csv(f) for f in filtered_df_files]
filtered_df_combined = pd.concat(filtered_df_list, ignore_index=True)

# 2. Concatenate all passedevents .npy arrays
passedevents_files = glob.glob(f"{output_dir}/passedevents_*.npy")

print(passedevents_files)
passedevents_list = [np.load(f, allow_pickle=True) for f in passedevents_files]
passedevents_combined = np.concatenate(passedevents_list)

# 3. Concatenate all passedhits .npy arrays
passedhits_files = glob.glob(f"{output_dir}/passedhits_*.npy")
passedhits_list = [np.load(f, allow_pickle=True) for f in passedhits_files]
passedhits_list = [arr for arr in passedhits_list if arr.size > 0]
print(passedhits_list)
passedhits_combined = np.concatenate(passedhits_list)

print("AAAAA")
# 3. Concatenate all failhits .npy arrays
failhits_files = glob.glob(f"{output_dir}/failedHits_*.npy")
failhits_list = [np.load(f, allow_pickle=True) for f in failhits_files]
print(failhits_list)
failhits_list = [arr for arr in failhits_list if arr.size > 0]
failhits_combined = np.concatenate(failhits_list)

# 3. Concatenate all allhits .npy arrays
allhits_files = glob.glob(f"{output_dir}/allhits_*.npy")
allhits_list = [np.load(f, allow_pickle=True) for f in allhits_files]
allhits_combined = np.concatenate(allhits_list)

# 4. Save merged outputs
filtered_df_combined.to_csv("/usera/dp728/run_dir/output/MergedHNL/merged_filtered_df.csv", index=False)
np.save("/usera/dp728/run_dir/output/MergedHNL/merged_passedevents.npy", passedevents_combined)
np.save("/usera/dp728/run_dir/output/MergedHNL/merged_passedhits.npy", passedhits_combined)
np.save("/usera/dp728/run_dir/output/MergedHNL/merged_failedhits.npy", failhits_combined)
np.save("/usera/dp728/run_dir/output/MergedHNL/merged_allhits.npy", allhits_combined)

#swap to plot terminology
passedHits = passedhits_combined
allhits = allhits_combined
failedHits = failhits_combined

#setup cavern for selection
import datetime
import pickle
cav =g.ATLASCavern()

pickleDir = "/usera/dp728/run_dir/pickles"
if not os.path.exists(pickleDir):
    os.makedirs(pickleDir)

print("Making RPC Layers...")
print(datetime.datetime.now())
if not os.path.exists(f"{pickleDir}/ANUBIS_RPCs_Layer0.pickle"): #or args.remake:
    #Layer of RPCs (Triplet) attached to the ceiling.
    RPC_Pos1 = cav.ANUBIS_RPC_positions(RPCx=1, RPCy=0.06, RPCz=1.8, overlapAngleXY=0, overlapZ=0, layerRadius = cav.archRadius, ID=0)
    with open(f"{pickleDir}/ANUBIS_RPCs_Layer0.pickle", "wb") as pkl:
        pickle.dump(RPC_Pos1, pkl)
else:
    with open(f"{pickleDir}/ANUBIS_RPCs_Layer0.pickle", "rb") as pkl:
        RPC_Pos1 = pickle.load(pkl)

if not os.path.exists(f"{pickleDir}/ANUBIS_RPCs_Layer1.pickle"): #or args.remake:
    #Singlet RPCs between the Triplets (40cm below Top Triplet).
    RPC_Pos2 = cav.ANUBIS_RPC_positions(RPCx=1, RPCy=0.06, RPCz=1.8, overlapAngleXY=0, overlapZ=0, layerRadius = cav.archRadius-0.40, ID=0)
    with open(f"{pickleDir}/ANUBIS_RPCs_Layer1.pickle", "wb") as pkl:
        pickle.dump(RPC_Pos2, pkl)
else:
    with open(f"{pickleDir}/ANUBIS_RPCs_Layer1.pickle", "rb") as pkl:
        RPC_Pos2 = pickle.load(pkl)

if not os.path.exists(f"{pickleDir}/ANUBIS_RPCs_Layer2.pickle"):# or args.remake:
    #Bottom Triplet RPCs (1m below Top Triplet).
    RPC_Pos3 = cav.ANUBIS_RPC_positions(RPCx=1, RPCy=0.06, RPCz=1.8, overlapAngleXY=0, overlapZ=0, layerRadius = cav.archRadius-1, ID=0)
    with open(f"{pickleDir}/ANUBIS_RPCs_Layer2.pickle", "wb") as pkl:
        pickle.dump(RPC_Pos3, pkl)
else:
    with open(f"{pickleDir}/ANUBIS_RPCs_Layer2.pickle", "rb") as pkl:
        RPC_Pos3 = pickle.load(pkl)

print(f"Cavern Bounds: {cav.CavernBounds}")
print(f"Cavern Corners: {cav.CavernCorners}")
print(f"IP Location: {cav.IP}")
print(f"Centre of Curvature: {cav.centreOfCurvature}")
print(cav.angles)

#this is to be cav.createSimpleRPCs but when movign into helper didnt have cav defined so this might be easier to just have a copy of createSimpleRPCs in helpers than calling it
ANUBISstations = cav.createSimpleRPCs([cav.archRadius, cav.archRadius-0.4, cav.archRadius-1], RPCthickness=0.06)
print(ANUBISstations)
minStationRadius = min(min(ANUBISstations["r"]))
#finding minimum point from probably ip as constraint for hit or not

ANUBISstations = RPC_Pos1
ANUBISstations.extend(RPC_Pos2)
ANUBISstations.extend(RPC_Pos3)
ANUBISstations = cav.convertRPCList(ANUBISstations)

ANUBISstations = cav.createSimpleRPCs([cav.archRadius, cav.archRadius-0.4, cav.archRadius-1], RPCthickness=0.06)
print(ANUBISstations)
minStationRadius = min(min(ANUBISstations["r"]))
print(minStationRadius)





import matplotlib.pyplot as plt

# Extract coordinates
x_vals = [pt[0] for pt in passedHits]
y_vals = [pt[1] for pt in passedHits]
z_vals = [pt[2] for pt in passedHits]

# Define bin counts
n_bins_x = 100
n_bins_y = 100
n_bins_z = 100

# Calculate min/max for each axis
#data
#min_x, max_x = min(x_vals), max(x_vals)
#min_z, max_z = min(z_vals), max(z_vals)
#min_y, max_y = min(y_vals), max(y_vals)
#anubis
min_x, max_x = cav.CavernX[0] - 5, cav.CavernX[1] + 5
min_y, max_y = cav.CavernY[0] - 5, (cav.archRadius + cav.centreOfCurvature["y"]) + 5
min_z, max_z = cav.CavernZ[0] - 5, cav.CavernZ[1] + 5

# Calculate bin widths
width_x = (max_x - min_x) / n_bins_x
width_y = (max_y - min_y) / n_bins_y
width_z = (max_z - min_z) / n_bins_z

# Pad ranges by half a bin
range_x = [min_x - width_x / 2, max_x + width_x / 2]
range_y = [min_y - width_y / 2, max_y + width_y / 2]
range_z = [min_z - width_z / 2, max_z + width_z / 2]

# Update hitBins dictionary
hitBins = {
    "x": n_bins_x,
    "y": n_bins_y,
    "z": n_bins_z,
    "widthX": width_x,
    "widthY": width_y,
    "widthZ": width_z,
    "rangeX": range_x,
    "rangeY": range_y,
    "rangeZ": range_z
}

from matplotlib.colors import LogNorm

cav.plotFullCavern(ANUBISrpcs=[RPC_Pos1, RPC_Pos2, RPC_Pos3], plotATLAS=True, plotFailed=False, 
                    hits={"passed": allhits, "failed": failedHits, "bins": hitBins}, suffix=f"precut") #{args.suffix}

cav.plotFullCavern(ANUBISrpcs=[RPC_Pos1, RPC_Pos2, RPC_Pos3], plotATLAS=True, plotFailed=False, 
                            hits={"passed": passedHits, "failed": failedHits, "bins": hitBins}, suffix=f"_WithHits") #{args.suffix}
print(datetime.datetime.now())