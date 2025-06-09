import sys, os
import csv
import yaml
import numpy as np
import pandas as pd
from scipy import constants

import zipfile
import gzip
import shutil 
import pickle
import defineGeometry as g
	
from random import randint

#import defineGeometry as g

#====================================================#
#===  Functions to process simulation information ===#
#====================================================#

def parseConfig(configFile):
    with open(configFile, 'r') as f:
        configuration = yaml.safe_load(f)

    return configuration

def createGeometry(geoFile, origin=[]):
    print(geoFile)
    
    if os.path.isfile(geoFile):
        with open(geoFile, "rb") as pkl:
            cavern = pickle.load(pkl)
            print("found geofile?")
    else:
        cavern = g.ATLASCavern()
        #ANUBISstations = cavern.createSimpleRPCs([cavern.archRadius-0.2, cavern.archRadius-0.6, cavern.archRadius-1.2], RPCthickness=0.06)
        # rng 0 and 1 and then making ot need 2 1s and and 1 0 
        # check names 
        # if naming isnt right
        # manually run a few of the lower lifetimes and one at a lower life time to check if there is an obvious difference 
        # double check the weighted c*tau case - to make sure its properly reweighting the lifetime.
        # take the c tau value for the given lifetime and see if it corresponds to a bit to the left for the ctau
        # warning any plotting things in hepmc are in mm - dataframes has the units at the top - 
        
        """
        #triplet
        if intersect
            RPC_number = 3
            tolerance = 2 
            pass = False
            hit = 0
            for rpc in RPC_number:
                number = randint(1, 100)
                if number < 98:
                    hit += 1

            if hit >= tolerance:
                pass = True
        """
            

        ANUBISstations = cavern.createSimpleRPCs([cavern.archRadius-0.2, cavern.archRadius-1.2], RPCthickness=0.06)
        cavern.RPCMaxRadius = cavern.archRadius-1.2-0.5
        if origin == "IP" or len(origin)==0:
            # Set the Cavern origin to the IP to align with simulation that puts it at the IP
            cavern.posOrigin = [cavern.IP["x"], cavern.IP["y"], cavern.IP["z"]]
        else:
            cavern.posOrigin = origin
        with open(geoFile, "wb") as pkl:
            pickle.dump(cavern, pkl)

    return cavern


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
                a = 0
                ds = 1
                #if "electron" in csvFile:
                if a == 1:
                    sampleDicts[line[0]]={"runName": line[1],
                                          "LLPmass": float(line[2]), # GeV
                                          "coupling": float(line[3]),
                                          "crossSection": float(line[4]),
                                          "LLPdecayWidth": float(line[5]), # GeV
                                          "LLPid": int(line[6]),
                                          "SimSwarmName": line[7],
                                          "sample": f"{line[7].split('_')[1]}_{line[7].split('_')[2]}"
                                        }
                elif ds == 1:
                    #DS case
                    sampleDicts[line[0]]={"runName": line[1],
                                        "Mzdmass": float(line[2]), # GeV
                                        "LLPmass": float(line[3]), # GeV
                                        "epsilon": float(line[4]), #for new UFO file change back if doesnt work
                                        "coupling": float(line[5]),
                                        "crossSection": float(line[6]),
                                        "LLPdecayWidth": float(line[9]), # GeV

                                        "LLPid": int(line[11]),
                                        "SimSwarmName": line[12],
                                        "sample": f"{line[12].split('_')[1]}"
                                        }

                else:
                    sampleDicts[line[0]]={"runName": line[1],
                                          "LLPmass": float(line[2]), # GeV
                                          "coupling": float(line[8]),
                                          "crossSection": float(line[4]),
                                          "LLPdecayWidth": float(line[5]), # GeV
                                          "LLPid": int(line[6]),
                                          "SimSwarmName": line[7],
                                          "sample": f"{line[7].split('_')[1]}_{line[7].split('_')[2]}"
                                        }
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
#def reweightDecayDistByLifetime(row, lifetime, seed=""):
#    if seed !="":
#        np.random.seed(seed)
    #print(lifetime)
    #print(row["boost"])
    # Reweight decay position using lifetime using Eqn 49.14 (pg 743) from the 2022 PDG https://inspirehep.net/files/850fdb9ba039c29dfbc1f071e6afd6e6
    return (np.random.exponential(scale = lifetime * row["boost"]) * row["beta"] * constants.c)*1000 # In mm

def reweightDecayDistByPosition(row, lifetime, seed=""):     
    # Reweight decay position using lifetime using Eqn 49.15 (pg 743) from the 2022 PDG https://inspirehep.net/files/850fdb9ba039c29dfbc1f071e6afd6e6
    return (np.random.exponential(scale = lifetime * row["boost"] * row["beta"] * constants.c))*1000 # In mm
            
def reweightDecayDistByRestLifetime(row, lifetime, seed=""):
 
    # Reweight decay position using lifetime using Eqn 49.14 (pg 743) from the 2022 PDG https://inspirehep.net/files/850fdb9ba039c29dfbc1f071e6afd6e6
    # But with the particle in the rest frame (boost=1) and then boosted.
    return (np.random.exponential(scale = lifetime * 1) * row["boost"] * row["beta"] * constants.c)*1000 # In mm



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
def checkLLPdecays(row, childrenOfLLPs):
    passed=False
    if row["nChildren"] == 1:
        passed = len(childrenOfLLPs[(childrenOfLLPs.eventNumber==row["eventNumber"]) & (childrenOfLLPs.particleIndex == row["childrenIndices"][0])]["PID"] == row["PID"])>=1
    elif row["nChildren"] > 1:
        passed = True
    else:
        passed = False

    return passed

def checkInCavern(row, geometry, maxRadius, decayVertex="decayVertex"):
    # Transform to the Cavern Centre Coordinates and m
    X, Y, Z = geometry.coordsToOrigin(row[decayVertex][0]*1E-3, row[decayVertex][1]*1E-3, row[decayVertex][2]*1E-3)
    return geometry.inCavern(X, Y, Z, maxRadius)

def checkInATLAS(row, geometry, trackingOnly=False, decayVertex="decayVertex"):
    # Transform to the Cavern Centre Coordinates and m
    X, Y, Z = geometry.coordsToOrigin(row[decayVertex][0]*1E-3, row[decayVertex][1]*1E-3, row[decayVertex][2]*1E-3)
    return geometry.inATLAS(X, Y, Z, trackingOnly)

def checkIntersectionsWithANUBIS(row, geometry, decayVertex="decayVertex"):
    # Transform to the Cavern Centre Coordinates and m
    X, Y, Z = geometry.coordsToOrigin(row[decayVertex][0]*1E-3, row[decayVertex][1]*1E-3, row[decayVertex][2]*1E-3)
    nInt, tempIntersections = geometry.intersectANUBISstationsSimple(X, Y, Z, geometry.ANUBIS_RPCs, verbose=False)

    intersections=[]
    for tI in tempIntersections:
        # Transforming back out of the Cavern Centre Coordinates and to mm
        intersections.append(geometry.reverseCoordsToOrigin(tI[0]*1E3, tI[1]*1E3, tI[2]*1E3))

    return intersections

def checkChildrenIntersections(row, geometry, LLPchildren, nHits=2, requireCharge=True, decayVertex="decayVertex"):
    tempChildren = LLPchildren[LLPchildren.LLPindex == row.name]

    tempChildren["intersectionsWithANUBIS"] = tempChildren.apply(checkIntersectionsWithANUBIS, args=(geometry, decayVertex), axis=1)
    if requireCharge:
        # Select Charged children
        tempChildren = tempChildren[(tempChildren.charge != 0) & (tempChildren.charge != -0.555)] 

    # Select LLP children with 2 or more intersections with ANUBIS stations i.e. 2 or more ANUBIS layers if there's no overlaps
    childrenIntersect = tempChildren[tempChildren.apply(lambda row: len(row["intersectionsWithANUBIS"]) >=2, axis=1)]

    # Return True to keep the LLP if there are more than {nHits} children that pass these requirements
    passed = len(childrenIntersect) >= nHits
    return passed 

def performGeometryCut(LLPsSel, maxRadius, trackingOnly=False, geometry="", decayVertex="decayVertex"):
    if geometry == "":
        geometry = createGeometry("")
    # Check LLP's decay vertex is in the Cavern Volume and within a maxRadius of the centre of curvature of the ceiling
    LLPsInCavern = LLPsSel[LLPsSel.apply(checkInCavern, args=(geometry, maxRadius, decayVertex), axis=1)]
    # Check that the LLP decay vertex is outside of the ATLAS detector
    LLPsNotInATLAS = LLPsInCavern[~LLPsInCavern.apply(checkInATLAS, args=(geometry, trackingOnly, decayVertex), axis=1)]
    # Check that the LLP direction would pass through at least two layers of the ANUBIS Detector
    temp = LLPsNotInATLAS.apply(checkIntersectionsWithANUBIS, args=(geometry, decayVertex), axis=1)
    if len(temp)!=0:
        LLPsNotInATLAS["intersectionsWithANUBIS"] = LLPsNotInATLAS.apply(checkIntersectionsWithANUBIS, args=(geometry, decayVertex), axis=1)
        LLPsIntersect = LLPsNotInATLAS[LLPsNotInATLAS.apply(lambda row: len(row["intersectionsWithANUBIS"]) >=2, axis=1)]
    else:
        LLPsIntersect=LLPsNotInATLAS


    #for batch running with condor
    #from the old version
    from datetime import datetime
    output_dir = "/usera/dp728/run_dir/output/HNL"
    
    #random id
    import uuid

    job_id = uuid.uuid4().hex[:6]
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")

    LLPsIntersect_filename  = f"{output_dir}/LLPsIntersect_{timestamp}_{job_id}.csv"
    
    LLPsIntersect.to_csv(LLPsIntersect_filename, index=False)
    #np.save(f"{output_dir}/passedhits_{timestamp}_{job_id}.npy", passedHits)
    #np.save(f"{output_dir}/allhits_{timestamp}_{job_id}.npy", allhits)
    #np.save(f"{output_dir}/failedHits_{timestamp}_{job_id}.npy", failedHits)
    #np.save(f"{output_dir}/passedevents_{timestamp}_{job_id}.npy", passedevents)

    print("check'/asd./'sa'd'")

    return LLPsIntersect

def checkDecayHits(LLPsGeo, LLPchildren, nHits = 2, each=True, requireCharge=True, geometry="", decayVertex="decayVertex"):
    if geometry == "":
        geometry = createGeometry("")

    LLPsDecayInANUBIS = LLPsGeo[LLPsGeo.apply(checkChildrenIntersections, args=(geometry, LLPchildren, nHits, requireCharge, decayVertex), axis=1)]

    return LLPsDecayInANUBIS 

def getMinDeltaR(row, sampleDFs, selection):
    # Get the Associated Jets/charged tracks for the event.
    dfBG = {"jet":  sampleDFs["finalStatePromptJets"][sampleDFs["finalStatePromptJets"].eventNumber == row["eventNumber"]],
            "tracks": sampleDFs["chargedFinalStates"][sampleDFs["chargedFinalStates"].eventNumber == row["eventNumber"]]
    }

    minDeltaR = []
    for BGtype in ["jet","tracks"]:
        if len(dfBG[BGtype])!=0:
            if BGtype == ["jet"]:
                # Require a minimum jet momenta and pt
                dfBG[BGtype] = dfBG[BGtype][(dfBG[BGtype].pt > selection["minPt"]["jet"]) & (dfBG[BGtype].p > selection["minP"]["jet"])]
            else:
                # Require a minimum pt for charged track
                dfBG[BGtype] = dfBG[BGtype][dfBG[BGtype].pt > selection["minPt"]["chargedTrack"]]

            dfBG[BGtype]["deltaEta"] = row.eta - dfBG[BGtype]["eta"] 
            dfBG[BGtype]["deltaPhi"] = row.phi - dfBG[BGtype]["phi"] 
            dfBG[BGtype]["deltaR"] = np.sqrt(np.power(dfBG[BGtype]["deltaEta"], 2) + np.power(dfBG[BGtype]["deltaPhi"], 2)) 
            minDeltaR.append(dfBG[BGtype]["deltaR"].min())
        else:
            # If there are no associated Jets/Charged Tracks the isolation condition is automatically met.
            #   - DeltaR can never be negative by our definition so use the -1 values to identify these cases.
            minDeltaR.append(-1) 

    return minDeltaR
"""
def findpartners(row,LLPsInATLAS,LLPsIsoAll,LLPs):


    # Check if any LLP in the same event passed the isolation cut
    same_event_passed = LLPsIsoAll[LLPsIsoAll["eventNumber"] == row["eventNumber"]]
    if  same_event_passed.empty:
        return False

    # Check if the LLP is in ATLAS and satisfies the eta cut
    in_atlas = (LLPsInATLAS["eventNumber"] == row["eventNumber"]) & (LLPsInATLAS.index == row.name)
    if in_atlas.any():
        eta = calculate_eta(row["px"], row["py"], row["pz"])
        if abs(eta) < 1.5:
            return True

    return False
    if LLPsIsoAll["eventNumber"] == row["eventNumber"]:
        if row in LLPsIsoAll:
            return True
        elif row in LLPsInATLAS:
            if calculate_eta(row["px"], row["py"], row["pz"]) < 1.5:    # pseudorapidity requirement
                return True
    return False


def includepartners(LLPsInATLAS,LLPsIsoAll,LLPs)
    LLPs=LLPs[LLPs.apply(findpartners, args=(LLPsInATLAS, LLPsIsoAll, LLPs), axis=1)]
"""
        
def includeAtlasPartners(LLPsIsoAll, LLPsInATLAS):
    # Get set of eventNumbers that hit anubis
    passed_events = set(LLPsIsoAll["eventNumber"])

    # Find partners that hit atlas
    # pseudorapidiity condition
    atlas_partners = LLPsInATLAS[
        (LLPsInATLAS["eventNumber"].isin(passed_events)) &
        (LLPsInATLAS["eta"].abs() < 1.5)
    ]

    # Combine the original passing LLPs with the newly included ATLAS partners
    LLPsCombined = pd.concat([LLPsIsoAll, atlas_partners])

    #checks event per llp
    LLPsperEvent = LLPsCombined["eventNumber"].value_counts()
    #list of all events with 2 llps
    Eventswith2LLPS = LLPsperEvent[LLPsperEvent == 2].index
    #filter original list to only keep ones with 2 llps
    LLPsFinal = LLPsCombined[LLPsCombined["eventNumber"].isin(Eventswith2LLPS)]

    return LLPsFinal


    """
    same_event_llps = LLPsIsoAll[LLPsIsoAll["eventNumber"] == row["eventNumber"]]
    
    # basically a mask for everthing u want
    for _, other in same_event_llps.iterrows():
        #
        if other.name == row.name:
            return True

        #since all rows are done each LLP in the pair will be checked itself
        #other LLP also in cut
        #if other.name in LLPsIsoAll.index:
        #    return True 

        #other one in AtlAS
        if other.name in LLPsInATLAS.index:
            if calculate_eta(other["px"],other["py"],other["pz"]) < 1.5: #pseudorapidity requirement
                return True
    """
    
    
