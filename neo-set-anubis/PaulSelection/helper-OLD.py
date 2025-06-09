import sys, os
import csv
import yaml
import numpy as np
#import pandas as pd
from scipy import constants

import zipfile
import gzip
import shutil 
#import defineGeometry as g

#====================================================#
#===  Functions to process simulation information ===#
#====================================================#

def parseConfig(configFile):
    with open(configFile, 'r') as f:
        configuration = yaml.safe_load(f)

    return configuration

def parseCSV(csvFile):
    # Assumes the csv has the following data structure
    #   - Run ID, Run Name, LLP Mass, Coupling, cross-section, LLP Decay Width, LLP ID, Simulation Swam Name

    #SampleDicts has a key that is the RunID and then a value that is a dictionary of all of its properties
    sampleDicts={}

    with open(csvFile) as f:
        reader = csv.reader(f, delimiter=",", quotechar='#')
        for line in reader: 
            if any("llp" in x for x in line):
                continue #Skip the header

            try:
                
                sampleDicts[line[0]]={"runName": line[1],
                                      "LLPmass": float(line[2]), # GeV
                                      "coupling": float(line[3]),
                                      "crossSection": float(line[4]),
                                      "LLPdecayWidth": float(line[5]), # GeV
                                      "LLPid": int(line[6]),
                                      "SimSwarmName": line[7],
                                      "sample": f"{line[7].split('_')[1]}_{line[7].split('_')[2]}"
                                    }
                """
                #DS case
                sampleDicts[line[0]]={"runName": line[1],
                                      "Mzdmass": float(line[2]), # GeV
                                      "LLPmass": float(line[3]), # GeV
                                      "epsilon": float(line[4]), #for new UFO file change back if doesnt work
                                      "coupling": float(line[5]),
                                      "crossSection": float(line[6]),
                                      "LLPid": int(line[11]),
                                      "SimSwarmName": line[12],
                                      "sample": f"{line[12].split('_')[1]}"
                                    }
                """
                
            except:
                print(f"Something went wrong when parsing csv, for this row:")
                print(line)

                raise Exception()


    return sampleDicts
"""
                sampleDicts[line[0]]={"runName": line[1],
                                      "LLPmass": float(line[2]), # GeV
                                      "coupling": float(line[3]),
                                      "crossSection": float(line[4]),
                                      "LLPdecayWidth": float(line[5]), # GeV
                                      "LLPid": int(line[6]),
                                      "SimSwarmName": line[7],
                                      "sample": f"{line[7].split('_')[1]}_{line[7].split('_')[2]}"
"""

def groupMassAndCoupling(sampleDicts):
    # Produces a nested dictionary with the form { sample: { LLP Mass: { Coupling : [List of Run IDs] } } }
    #   - Allows for the collation of simulations with the same input parameters for higher statistics runs.

    # Get sampleDicts if csv file is provided
    if ".csv" in sampleDicts:
        sampleDicts = parseCSV(sampleDicts)

    outputDict={}
    for ID, infoDict in sampleDicts.items():
        if infoDict["sample"] not in outputDict.keys():
            outputDict[infoDict["sample"]]={}

        if infoDict["LLPmass"] not in outputDict[infoDict["sample"]].keys():
            outputDict[infoDict["sample"]][infoDict["LLPmass"]]={}

        if infoDict["coupling"] not in outputDict[infoDict["sample"]][infoDict["LLPmass"]].keys():
            outputDict[infoDict["sample"]][infoDict["LLPmass"]][infoDict["coupling"]] = []

        outputDict[infoDict["sample"]][infoDict["LLPmass"]][infoDict["coupling"]].append(ID)

    return outputDict

def is_hepmc_zipped(filename):
    zip_extensions = ['.gz', '.gzip']
    ext = os.path.splitext(filename)[1]
    return ext.lower() in zip_extensions

def unzip_file(zip_file, out_file):
    with gzip.open(zip_file, 'rb') as f_in:
        with open(out_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

#====================================================#
#===  Functions to process Simulation variables   ===#
#====================================================#


#=========== Need to verify these functions -- from Toby #====================

# Convert r, eta, phi to Cartesian coordinates          
def to_cartesian(r, eta, phi):
    x = r * np.sin(to_theta(eta)) * np.cos(phi)
    y = r * np.sin(to_theta(eta)) * np.sin(phi)
    z = r * np.cos(to_theta(eta))
    return [x, y, z]

def to_spherical(x, y, z):
    # Perform conversion
    r = np.sqrt(np.square(x) + np.square(y) + np.square(z))
    theta = np.arctan2(np.sqrt(np.square(x) + np.square(y)), z)
    if x > 0:
        phi = np.arctan(y / x)
    elif x < 0 and y >= 0:
        phi = np.arctan(y / x) + np.pi
    elif x < 0 and y < 0:
        phi = np.arctan(y / x) - np.pi
    elif x == 0 and y > 0:
        phi = np.pi / 2
    elif x == 0 and y < 0:
        phi = -np.pi / 2
    else:
        phi = np.nan

    # Scale phi to [-pi, pi]
    phi = scale_phi(phi)

    return r, to_eta(theta), phi

# Scale phi to [-pi, pi]
def scale_phi(phi):
    return ((phi + np.pi) % (2 * np.pi)) - np.pi

# Calculate Delta R = sqrt( (Delta eta)^2 + (Delta phi)^2) of a jet or charged tracj
def calculate_deltaR(deltaeta, deltaphi): 
    deltaR = np.sqrt(np.square(deltaeta) + np.square(deltaphi))
    return deltaR

# Calculate boost (gamma)
#   - NOTE: Why not use |p| directly? Should be equivalent?
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

# Convert eta to theta
def to_theta(eta):
    return 2 * np.arctan(np.exp(-eta))

#calculate phi (azimuthal angle, i.e. angle from x axis): Range [-pi, pi]
def calculate_phi(px, py):
    if px > 0:
        phi = np.arctan(py / px)
    elif px < 0 and py >= 0:
        # As px<0 np.arctan(py/px) < 0 
        phi = np.arctan(py / px) + np.pi 
    elif px < 0 and py < 0:
        # -np.pi used here as we use phi in the range [-pi,pi]
        phi = np.arctan(py / px) - np.pi 
    elif px == 0 and py > 0:
        phi = np.pi / 2
    elif px == 0 and py < 0:
        phi = -np.pi / 2
    else:
        phi = np.nan
    return phi

def calculate_eta(px, py, pz): #(+ve z means +ve eta)
    if (np.sqrt(np.square(px) + np.square(py) + np.square(pz)) - pz) == 0:
        if (np.sqrt(np.square(px) + np.square(py) + np.square(pz)) + pz) == 0:
            return np.nan # In this case the pseudorapidity is undefined. 
        else:
            return 0 # For a particle perfectly transverse pseudorapidity is 0.
    
    if ( ((1/2)*np.log( (np.sqrt(np.square(px) + np.square(py) + np.square(pz)) + pz)/(np.sqrt(np.square(px) + np.square(py) + np.square(pz)) - pz) )) > 1000000):
         return 1000000
    elif ( ((1/2)*np.log( (np.sqrt(np.square(px) + np.square(py) + np.square(pz)) + pz)/(np.sqrt(np.square(px) + np.square(py) + np.square(pz)) - pz) )) < -1000000):
         return -1000000
    else:
         return (1/2)*np.log( (np.sqrt(np.square(px) + np.square(py) + np.square(pz)) + pz)/(np.sqrt(np.square(px) + np.square(py) + np.square(pz)) - pz) )

def calculate_pt(px, py):
    return np.sqrt(np.square(px) + np.square(py))

#=======================================================================

# The below 3 reweighting equations work assuming you provide it with a row from a pandas dataframe
def reweightDecayDistByLifetime(row, lifetime):
    # Reweight decay position using lifetime using Eqn 49.14 (pg 743) from the 2022 PDG https://inspirehep.net/files/850fdb9ba039c29dfbc1f071e6afd6e6
    return (np.random.exponential(scale = lifetime * row["boost"]) * row["beta"] * constants.c) #/1000 # In mm

def reweightDecayDistByPosition(row, lifetime):
    # Reweight decay position using lifetime using Eqn 49.15 (pg 743) from the 2022 PDG https://inspirehep.net/files/850fdb9ba039c29dfbc1f071e6afd6e6
    return (np.random.exponential(scale = lifetime * row["boost"] * row["beta"] * constants.c))/1000 # In mm
            
def reweightDecayDistByRestLifetime(row, lifetime):
    # Reweight decay position using lifetime using Eqn 49.14 (pg 743) from the 2022 PDG https://inspirehep.net/files/850fdb9ba039c29dfbc1f071e6afd6e6
    # But with the particle in the rest frame (boost=1) and then boosted.
    return (np.random.exponential(scale = lifetime * 1) * row["boost"] * row["beta"] * constants.c)/1000 # In mm

# This assumes you provide it with a row from a pandas dataframe, and assumes the reweightedDecayDistance is in mm.
def getReweightedDecayPosition(row, decayDistanceColumn):
    return (list(to_cartesian(row[decayDistanceColumn], row["eta"], row["phi"]))) # In mm  

# Translate the production or decay vertices of LLP decay products 
def getDecayVertexTranslation(row, LLPdf, vertexColumn):
    if( (row["nChildren"]==0) or (row["status"] == 1)): # FourVector is x,y,z,t in mm and s
        return (-1,-1,-1,-1) # Corresponds to a final state particle that does not decay

    translation = LLPdf.loc[row['LLPindex']]['decayVertex_translation']
    translatedVertex = tuple( [ a + b for a, b in zip(row[vertexColumn], translation) ] )

    return translatedVertex

def createSimpleRPCs(self, radii, RPCthickness=0.06, origin=[]):
    if len(origin)==0: # Assume origin of IP unless otherwise stated
        origin = (self.IP["x"], self.IP["y"], self.IP["z"])

    RPCs={"r": [], "theta": [], "phi": []}

    # Radii are assumed to be relative to the centre of curvature
    for r in radii:
        RPCs["r"].append( [r - RPCthickness, r] )
        RPCs["theta"].append([min(cav.angles["theta"]), max(cav.angles["theta"])])

        if r < abs(self.CavernX[0]-origin[0]):
            RPCs["phi"].append([0,2*np.pi]) # As in this case you get a circle within the ATLAS Cavern
        else:
            phiList=[]
            for i in [1,0]:
                vec1 = self.createVector([self.CavernX[i],np.sqrt((r**2) - (self.CavernX[i] **2))], [origin[0], origin[1]])
                vec2 = self.createVector([self.CavernX[1], origin[1]], [origin[0], origin[1]])
                tempPhi = np.dot(vec1, vec2) / ( np.linalg.norm(vec1) * np.linalg.norm(vec2))
                phiList.append(np.arccos(np.clip(tempPhi, -1, 1))) # In radians
            
            RPCs["phi"].append(phiList)

    return RPCs


#====================================================#
#===           Function for Selections            ===#
#====================================================#

#TODO: Outline the Geometry Cuts
def performGeometryCut(LLPsSel,cav, ANUBISstations,minStationRadius,RPC_Pos1, RPC_Pos2, RPC_Pos3, plot = False):
    print("2")
    import datetime
    import pickle
    import argparse

    parser = argparse.ArgumentParser(description='Provide an input file')
    #args = parser.parse_args()
    

    cav.createCavernVault()


    # Create ANUBIS RPC Layers
    #   - Afterwards save to a pickle file. 
    #   - Reload RPCs from pickle file if it exists.


    print("Checking hit intersections...")
    hitStart1=datetime.datetime.now()
    print(hitStart1)
    hitBins = {"x": 100, "y": 100, "z": 100}
    # Check whether a set of points intersect ANUBIS
    """
    x = np.linspace(cav.CavernX[0] - 5, cav.CavernX[1] + 5, hitBins["x"])
    y = np.linspace(cav.CavernY[0] - 5, (cav.archRadius + cav.centreOfCurvature["y"]) + 5, hitBins["y"])
    z = np.linspace(cav.CavernZ[0] - 5, cav.CavernZ[1] + 5, hitBins["z"])

    hitBins["widthX"] = abs(x[1]-x[0])
    hitBins["widthY"] = abs(y[1]-y[0])
    hitBins["widthZ"] = abs(z[1]-z[0])
    hitBins["rangeX"] = [x[0]-(hitBins["widthX"]/2), x[-1]+(hitBins["widthX"]/2)]
    hitBins["rangeY"] = [y[0]-(hitBins["widthY"]/2), y[-1]+(hitBins["widthY"]/2)]
    hitBins["rangeZ"] = [z[0]-(hitBins["widthZ"]/2), z[-1]+(hitBins["widthZ"]/2)]
    """
    #hitbins are probably just setting up the bins and binwidths for the actual hits
    #anubisstations then just makes a list of all the stations
    """
    nHits=0
    passedHits, failedHits = [], []
    for X in x:
        for Y in y:
            for Z in z:
                #print(f"Before Intersection: {datetime.datetime.now()}")
                inCavern = cav.inCavern(X, Y, Z)
                inATLAS = cav.inATLAS(X, Y, Z, trackingOnly=False)
                intANUBIS = cav.intersectANUBISstations(X, Y, Z, ANUBISstations, origin=[])

                if ( inCavern and not inATLAS and (len(intANUBIS[1]) >= 2) ):
                    passedHits.append((X,Y,Z))
                else:
                    failedHits.append((X,Y,Z))
                #print(f"After Intersection: {datetime.datetime.now()}")
                print(f"{nHits}/{len(x)*len(y)*len(z)}")
                nHits+=1
                #input("...")

    """
    hitEnd1=datetime.datetime.now()
    print(hitEnd1)
    #print(f"Took {hitEnd1 - hitStart1}")

    #this is to be cav.createSimpleRPCs but when movign into helper didnt have cav defined so this might be easier to just have a copy of createSimpleRPCs in helpers than calling it

    #finding minimum point from probably ip as constraint for hit or not

    nHits=0
    """
    passedHits, failedHits = [], []
    for X in x:
        for Y in y:
            for Z in z:
                print(f"{nHits}/{len(x)*len(y)*len(z)}", end="\r", flush=True) #progress bar?
                inCavern = cav.inCavern(X, Y, Z, maxRadius=minStationRadius - 0.20)
                inATLAS = cav.inATLAS(X, Y, Z, trackingOnly=True)
                intANUBIS = cav.intersectANUBISstationsSimple(X, Y, Z, ANUBISstations)

                if ( inCavern and not inATLAS and (len(intANUBIS[1]) >= 2) ):
                    passedHits.append((X,Y,Z))
                else:
                    failedHits.append((X,Y,Z))
                nHits+=1
    """
    
    import pandas as pd
    hits_df = LLPsSel
    print("this is hits inputed to selection")
    print(hits_df)
    import ast
    #was for reading csv file
    #hits_df['decayVertexParsed'] = hits_df['decayVertex'].apply(ast.literal_eval)
    
    allhits = []
    passedHits, failedHits = [], []
    passedevents, failedevents = [], []
    eventeta = []
    eventphi = []

    mask = []
    for i, row in hits_df.iterrows():
        #X, Y, Z = row['decayVertexParsed'][0], row['decayVertexParsed'][1], row['decayVertexParsed'][2]
        X, Y, Z = row['decayVertex_weighted'][0], row['decayVertex_weighted'][1], row['decayVertex_weighted'][2]
        #X, Y, Z = row['decayVertex_restWeighted'][0], row['decayVertex_restWeighted'][1], row['decayVertex_restWeighted'][2]
        X = X-1.7 
        Y = Y-cav.CavernYLength/2 + 11.37
        """
        print("this is x")
        print(X)
        print("this is row['decayVertexParsed']")
        print(row['decayVertexParsed'])
        print("this is row['decayVertexParsed'][2]")
        print(row['decayVertexParsed'][2])
        """
        r, theta, phi = row['decayVertexDist_weighted'], row['theta'],row['phi']


        inCavern = cav.inCavern(X, Y, Z, maxRadius=minStationRadius - 0.20)
        inATLAS = cav.inATLAS(X, Y, Z, trackingOnly=True)
        intANUBIS = cav.intersectANUBISstationsSimple2(X, Y, Z,r, theta, phi, ANUBISstations)

        if ( inCavern and not inATLAS and (len(intANUBIS[1]) >= 2) ):
            passedHits.append((X,Y,Z))
            passedevents.append(row["eventNumber"])

        else:
            failedHits.append((X,Y,Z))
        allhits.append((X,Y,Z))

        condition =inCavern and not inATLAS and (len(intANUBIS[1]) >= 2)
        mask.append(condition)
        #eventeta.append(row["eta"])
        #eventphi.append(row["phi"])

        if i % 1000 == 0:
            print(f"Checked {i}/{len(hits_df)} hits", end="\r", flush=True)
    filtered_df = hits_df[mask].reset_index(drop=True)
    hitEnd2=datetime.datetime.now()
    print(hitEnd2)
    print(f"Took {hitEnd2 - hitEnd1}")
    print(passedevents)
    print(allhits)
    print(f"Passed: {len(passedHits)} ({len(passedHits)/(len(passedHits)+len(failedHits))}%)  | Failed: {len(failedHits)} ({len(failedHits)/(len(passedHits)+len(failedHits))}%)")
    print("first check")
    if plot:
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
        """
        plt.figure(figsize=(10,6))
        cuteventeta = [e for e in eventeta if -15 <= e <= 15]
        plt.hist(cuteventeta, bins=100, color='purple', alpha=0.75)
        plt.xlabel("η (pseudorapidity)")
        plt.ylabel("Number of hits")
        plt.title("Pseudorapidity Distribution")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("pseudorapidity_distribution.pdf")
        plt.show()

        plt.figure(figsize=(10,6))
        plt.hist(eventphi, bins=90, color='orange', alpha=0.75)
        plt.xlabel("ϕ (phi) [radians]")
        plt.ylabel("Number of hits")
        plt.title("Azimuthal Angle (ϕ) Distribution")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("phi_distribution.pdf")
        plt.show()
        """

    #for batch running with condor
        
    from datetime import datetime
    output_dir = "/usera/dp728/run_dir/output/HNL"
    
    #random id
    import uuid

    job_id = uuid.uuid4().hex[:6]
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")

    filtered_df_filename  = f"{output_dir}/filtered_df_{timestamp}_{job_id}.csv"

    filtered_df.to_csv(filtered_df_filename, index=False)
    np.save(f"{output_dir}/passedhits_{timestamp}_{job_id}.npy", passedHits)
    np.save(f"{output_dir}/allhits_{timestamp}_{job_id}.npy", allhits)
    np.save(f"{output_dir}/failedHits_{timestamp}_{job_id}.npy", failedHits)
    np.save(f"{output_dir}/passedevents_{timestamp}_{job_id}.npy", passedevents)

    print("check'/asd./'sa'd'")


    return filtered_df, failedHits # return itself for now to be able to verify selection chain function

#TODO: Outline the Cuts based on the LLP children
def checkDecayHits(LLPsGeo, cav, ANUBISstations, minStationRadius, RPC_Pos1, RPC_Pos2, RPC_Pos3, nHits=2, plot=True, each=True):

    df = LLPsGeo
    from collections import defaultdict

    passing_event_numbers = set()
    intersection_counts = defaultdict(int)

    passedHits = []

    # This assumes your df is sorted by eventNumber, but doesn't have to be
    for idx, row in df.iterrows():
        event = row['eventNumber']
        
        # Skip if we've already found 2 intersections for this event
        if event in passing_event_numbers:
            continue
        
        # Compute intersection
        X, Y, Z = row['decayVertex_weighted'][0], row['decayVertex_weighted'][1], row['decayVertex_weighted'][2]
        X = X - 1.7
        Y = Y - cav.CavernYLength / 2 + 11.37
        r, theta, phi = row['decayVertexDist_weighted'], row['theta'], row['phi']
        intANUBIS = cav.intersectANUBISstationsSimple2(X, Y, Z, r, theta, phi, ANUBISstations)
        
        #one valid interaction when hits 2 anubis rpcs 
        if len(intANUBIS[1]) >= 2:
            intersection_counts[event] += 1
        
        # If we've found 2 intersecting tracks, mark event as passing
        if intersection_counts[event] >= 2:
            passedHits.append((X,Y,Z))
            passing_event_numbers.add(event)

    # Final filtering step
    filtered_df = df[df['eventNumber'].isin(passing_event_numbers)].reset_index(drop=True)


    #grouped = filtered_df.groupby('eventNumber')

            
    #print(f"Passed events : {len(filtered_df)} ({len(filtered_df)/(len(df))}%)  | Failed events: {len(df)-len(filtered_df)} ({(len(df)-len(filtered_df))/(len(df))}%)")
    if len(df) == 0:
        print("⚠️ No events in input dataframe (df). Skipping percentage breakdown.")
        print(f"Passed events: {len(filtered_df)} | Failed events: 0")
    else:
        print(f"Passed events : {len(filtered_df)} ({len(filtered_df)/len(df)*100:.2f}%)  | Failed events: {len(df)-len(filtered_df)} ({(len(df)-len(filtered_df))/len(df)*100:.2f}%)")

    if plot:
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
        failedHits = []
        cav.plotFullCavern(ANUBISrpcs=[RPC_Pos1, RPC_Pos2, RPC_Pos3], plotATLAS=True, plotFailed=False, 
                                    hits={"passed": passedHits, "failed": failedHits, "bins": hitBins}, suffix=f"_WithHits_with2interseciton") #{args.suffix}
        #print(datetime.datetime.now())
        """
        plt.figure(figsize=(10,6))
        cuteventeta = [e for e in eventeta if -15 <= e <= 15]
        plt.hist(cuteventeta, bins=100, color='purple', alpha=0.75)
        plt.xlabel("η (pseudorapidity)")
        plt.ylabel("Number of hits")
        plt.title("Pseudorapidity Distribution")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("pseudorapidity_distribution.pdf")
        plt.show()

        plt.figure(figsize=(10,6))
        plt.hist(eventphi, bins=90, color='orange', alpha=0.75)
        plt.xlabel("ϕ (phi) [radians]")
        plt.ylabel("Number of hits")
        plt.title("Azimuthal Angle (ϕ) Distribution")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("phi_distribution.pdf")
        plt.show()
        """


    return filtered_df # return itself for now to be able to verify selection chain function
