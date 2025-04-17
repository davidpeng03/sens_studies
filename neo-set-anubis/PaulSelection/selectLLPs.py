import sys, os
import subprocess
import numpy as np
import pandas as pd
try:
    import pyhepmc
except:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pyhepmc"])
    import pyhepmc
try:
    import fastjet._pyjet
except:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "fastjet"])
    import fastjet._pyjet
from particle import PDGID, Particle, Charge
import datetime
import argparse 
import awkward as ak
import matplotlib
import matplotlib.pyplot as plt
import pickle

import helper as h
import defineGeometry as g
import yaml
from scipy import constants
import plotting as myplot
import json

def createDFfromHEPMC(fileName):
    print(f"== Creating df from: {fileName}...")
    print("ASSDJSNADIASNIDNSIADNISAND")
    columns=["eventNumber","particleIndex","px","py","pz","pt","E","mass","prodVertex","prodVertexDist","decayVertex","decayVertexDist",
             "boost","phi","eta","METx","METy","MET","theta","beta","PID","charge","nParents","nChildren","parentIndices","childrenIndices",
             "weight","status","ctau"]

    hepMCdict={}
    for col in columns:
        hepMCdict[col]=[]

    unknownPIDs=[]
    
    start=datetime.datetime.now()
    with pyhepmc.open(fileName) as f:
        eventNumber=0
        for event in f:
            #eventNumber = event.event_number # The HEPMC Event number appears to reset after 124 events, and would no longer be a unique identifier.
            if eventNumber%100 == 0:
                print(f"Event {eventNumber}: time elapsed {datetime.datetime.now()-start}")
            for p in event.particles:
                # For each particle in each event save all relevant properties - apply cuts later 
                hepMCdict["eventNumber"].append(eventNumber)
                hepMCdict["particleIndex"].append(p.id)
                # A unique ID per event can be composed of the event number and particle index
                hepMCdict["px"].append(p.momentum[0])
                hepMCdict["py"].append(p.momentum[1])
                hepMCdict["pz"].append(p.momentum[2])
                pt = h.calculate_pt(p.momentum[0], p.momentum[1])
                hepMCdict["pt"].append(pt)
                hepMCdict["E"].append(p.momentum[3])
                hepMCdict["mass"].append(p.generated_mass)
                prodVertex = p.production_vertex.position
                
                hepMCdict["prodVertex"].append( (prodVertex.x, prodVertex.y, prodVertex.z, prodVertex.t) )
                hepMCdict["prodVertexDist"].append(p.production_vertex.position.length())
                if len(p.children) != 0: # Some particles may not decay and so put in a placeholder position
                    endVertex = p.end_vertex.position
                    hepMCdict["decayVertex"].append( (endVertex.x, endVertex.y, endVertex.z, endVertex.t) ) # FourVector is x, y, z, t - in mm and s
                    hepMCdict["decayVertexDist"].append(endVertex.length()) #Only take the position part of the FourVector
                else: 
                    hepMCdict["decayVertex"].append((-1,-1,-1,-1)) # FourVector is x, y, z, t - in mm and s
                    hepMCdict["decayVertexDist"].append(0)
                if eventNumber % 100 == 0:
                    print(f"Event {eventNumber}: decayVertex example {endVertex.length()}")
                    #temporaryvariableidk = hepMCdict["decayVertexDist"]
                    #print(f"check hepMCdict: {temporaryvariableidk}")
                hepMCdict["phi"].append(h.calculate_phi(p.momentum[0], p.momentum[1]))
                eta = h.calculate_eta(p.momentum[0], p.momentum[1],p.momentum[2])
                hepMCdict["eta"].append(eta)
                boost = h.calculate_boost(pt,eta,p.generated_mass)
                hepMCdict["boost"].append(boost)
                hepMCdict["theta"].append(h.to_theta(eta))
                beta = h.calculate_beta(boost)
                hepMCdict["beta"].append(beta)
                hepMCdict["nParents"].append(len(p.parents))
                hepMCdict["nChildren"].append(len(p.children))
                hepMCdict["parentIndices"].append([x.id for x in p.parents])
                hepMCdict["childrenIndices"].append([x.id for x in p.children])
                hepMCdict["weight"].append(event.weight())
                hepMCdict["status"].append(p.status)
                hepMCdict["ctau"].append( (hepMCdict["decayVertexDist"][-1])/(boost*beta) )

                if eventNumber % 100 == 0:
                    print("--------------------------------------")
                    print(f"Event {eventNumber}: ctau example {endVertex.length()}")
                    #nexttermpvariable = hepMCdict["decayVertexDist"][-1]
                    #print(F"check hepMCdict: {nexttermpvariable}")

                pid = p.pid
                hepMCdict["PID"].append(pid)

                # Useful list of IDs: http://home.thep.lu.se/~bierlich/misc/keyword-manual/ParticleData.html
                if abs(pid) in [5214, 5212, 9940003, 9900012, 9940011, 9942003, 9942103, 9942033]: 
                    # Special cases for known neutral particles 
                    charge=0
                elif abs(pid) in [543, 5314, 20413, 9932103]:
                    # Special cases for known charged particles 
                    charge=-1 if pid < 0 else 1
                elif abs(pid) in [35]:
                    charge=0
                else:
                    try:
                        part = Particle.from_pdgid(p.pid)
                        charge = part.charge
                        
                    except:
                        #print(f"Unknown PID: {p.pid}, cannot assign charge properly please implement manually")
                        #print(f"Particle will be assigned a charge of -0.555 in the meantime")
                        unknownPIDs.append(p.pid)
                        charge = -0.555

                hepMCdict["charge"].append(charge)

                # Placeholder - later will use (pt of neutrinos) - (pt of charged + neutral tracks) to get MET for LLP  
                hepMCdict["MET"].append(0.0)
                hepMCdict["METx"].append(0.0)
                hepMCdict["METy"].append(0.0)
            
            eventNumber+=1

    print(f"Unknown PIDs: {list(set(unknownPIDs))}")
    end = datetime.datetime.now()
    print(f"total time taken: {end-start}")

    df = pd.DataFrame.from_dict(hepMCdict)
    # Ensure the MET branches are treated as floats
    df["MET"] = df["MET"].astype(float)
    df["METx"] = df["METx"].astype(float)
    df["METy"] = df["METy"].astype(float)
    print(df)
    df.to_csv("./tempDF.csv")
    
    return df

def childrenHunt(dataFrame, eventNumber, index, childList=[]):
    # NOTE: Mixing dataframe indices and particleIndices??
    indices = dataFrame.index[(dataFrame["eventNumber"] == eventNumber) & (dataFrame["particleIndex"] == index)].tolist()
    childrenIndexList = dataFrame.loc[indices].childrenIndices.tolist()

    if len(childrenIndexList)==0:
        return childList
    else:
        for childrenIndices in childrenIndexList:
            for idx in childrenIndices:
                if idx in childList:
                    continue
                childList.append(idx)
                childrenHunt(dataFrame, eventNumber, idx, childList)
            return childList

 

def createSampleDataFrames(df, LLPid, minPt):
    print(f"Creating sample dataframes...")
    print(datetime.datetime.now())
    # Get a list of event numbers:
    eventNos = list(set(df['eventNumber']))

    # Get all the final state particles
    finalStates = df[(df.nChildren == 0) & (df.status==1)]

    print(LLPid)
    # Get the LLPs with the given ID and their indices
    LLPindices = df.index[df["PID"]==LLPid].tolist()
    LLPs = df[df.PID == LLPid]
    print("________________________________________________________________")
    print(f"Average LLP Decay Distance (mm): {LLPs['decayVertexDist'].mean()}")
    print(f"Average LLP ctau (mm): {LLPs['ctau'].mean()}")
    print("________________________________________________________________")

    # Get the decay particles of the LLP
    #   - To allow for selection based on their tracks

    # 2 llps per event so theres hould be 2*events 
    # number of cbildren associated with a particle
    # can produce child thats itself 
    #
    #take the first 2-3 events just do 2 or 3 events 
    #save to csv file and look at it 
    #
    #
    #


    LLPchildrenDict={"eventNumber": [], "childrenIndices": [], "LLPindex": []}
    
    for index in LLPindices:
        tempList=[]
        tempEvtNo = df.loc[index].eventNumber
        tempParticleIndex = df.loc[index].particleIndex
        # Children Hunt gives the list of particle IDs of the children -- resets per event.
        tempChildIndices = childrenHunt(df, tempEvtNo, tempParticleIndex, tempList) #Returns particleID list
        print(index)
        for chIdx in tempChildIndices:
            LLPchildrenDict["eventNumber"].append(tempEvtNo)
            LLPchildrenDict["childrenIndices"].append(chIdx)
            LLPchildrenDict["LLPindex"].append(index)
    #print(len(LLPindices))

    # The above gives the per event children indices, this gives a list of the dataframe indices for those
    childrenDFindices = []
    originatingLLP = []
    for i in range(len(LLPchildrenDict["eventNumber"])):
        childrenDFindices.extend(df[(df.eventNumber == LLPchildrenDict["eventNumber"][i]) &\
                                    (df.particleIndex == LLPchildrenDict["childrenIndices"][i])].index.tolist())
        originatingLLP.append(LLPchildrenDict["LLPindex"][i])

    LLPchildren = df.iloc[childrenDFindices]
    #LLPchildren = LLPchildren[~LLPchildren.index.duplicated(keep='first')]


    # Sometimes a MC generator can propagate a particle without a decay such that the 'decay chain' is of the form:
    #   X ---> X ---> X ---> decay products
    # Therefore, the "True" LLPs are the first in this chain as that represents the actual production.
    # i.e. if an LLP is also in the LLP children list, then it should be removed.
    #so makes a list of all the 
    newLLPindices = [ x for x in LLPindices if x not in childrenDFindices]

    #newLLPindices = LLPindices

    """
    newLLPindicesfordistances = [
    x for x in LLPindices
    if not any(chIdx in LLPindices for chIdx in childrenHunt(df, df.loc[x].eventNumber, df.loc[x].particleIndex))
    ]
    """
    
    print(newLLPindices)
    LLPs = df.iloc[newLLPindices]
    #LLPSfordistancecalc = df.iloc[newLLPindicesfordistances]
    # Add a column to identify which LLP the decay products came from originally.
    childrenDF = pd.DataFrame(originatingLLP, columns=["LLPindex"], index=childrenDFindices)
    childrenDF = childrenDF[childrenDF.LLPindex.isin(newLLPindices)]
    
    print(f"LLPchildren rows: {len(LLPchildren)}")
    print(f"childrenDF['LLPindex'] length: {len(childrenDF['LLPindex'])}")
    print(childrenDF.columns)
    print(LLPchildren.columns)

    LLPchildren["LLPindex"] = childrenDF["LLPindex"].values
    LLPchildren = LLPchildren[LLPchildren.PID != LLPid] # Remove the LLP from the decay products
    # i believe the way this works is that it deletes all of the LLPs in the chains but then reassigns the  




    # Select Final states that don't originate from the LLP
    finalStates_NoLLP = finalStates[~finalStates.index.isin(newLLPindices)] # Remove the LLPs 
    finalStates_NoLLP = finalStates_NoLLP[~finalStates_NoLLP.index.isin(childrenDFindices)] # Remove the LLP decay products
    finalStates_NoLLP = finalStates_NoLLP[finalStates_NoLLP.PID != LLPpid] # Sanity check - Ensure that the LLP isn't included
    
    # Select the final State neutrinos
    finalStates_Neutrinos = finalStates_NoLLP[finalStates_NoLLP.PID.isin([12,14,16,18])]
    finalStates_NoLLP = finalStates_NoLLP[~finalStates_NoLLP.PID.isin([12,14,16,18])] # Remove the neutrinos with their PIDs 

    # Select Charged Final states that don't originate from the LLP
    #   Require it to be prompt, i.e. have a production vertex <10mm from the IP
    #   Require a minimum pT given in the config file
    chargedFinalStates = finalStates_NoLLP[(finalStates_NoLLP.charge != 0) & (finalStates_NoLLP.charge != -0.555)] 
    chargedFinalStates = chargedFinalStates[chargedFinalStates.prodVertexDist < 10]
    chargedFinalStates = chargedFinalStates[chargedFinalStates.pt > minPt["chargedParticle"]]

    # Select Neutral Final states that don't originate from the LLP
    #   Require it to be prompt, i.e. have a production vertex <10mm from the IP
    neutralFinalStates = finalStates_NoLLP[finalStates_NoLLP.charge == 0] 
    neutralFinalStates = neutralFinalStates[neutralFinalStates.prodVertexDist < 10] 

    # For each Event get the associated pt of the final states and assign this to the LLP as MET
    #   - This doesn't include neutrinos, the LLPs themselves or the LLPs' decay products
    #   - We care about the total MET in an event for selection purposes, and so each LLP in an event will be assigned the total MET for ease of reference. 
    for event in eventNos:
        METx = finalStates_NoLLP[finalStates_NoLLP.eventNumber == event]["px"].sum()
        METy = finalStates_NoLLP[finalStates_NoLLP.eventNumber == event]["py"].sum()
        if len(LLPs[LLPs.eventNumber == event])!=0:
            LLPs.loc[LLPs.loc[LLPs.eventNumber == event].index, "METx"] = METx
            LLPs.loc[LLPs.loc[LLPs.eventNumber == event].index, "METy"] = METy
            LLPs.loc[LLPs.loc[LLPs.eventNumber == event].index, "MET"] = h.calculate_pt(METx, METy)

    outputDict={"finalStates": finalStates,
                "LLPs": LLPs,
                "LLPchildren": LLPchildren, 
                "finalStates_NoLLP": finalStates_NoLLP, # Prompt final states
                "finalStates_Neutrinos": finalStates_Neutrinos, # Prompt neutrino final states
                "chargedFinalStates": chargedFinalStates, # Charged Prompt Final States
                "neutralFinalStates": neutralFinalStates,
                #"LLPSfordistancecalc": LLPSfordistancecalc,                
    }
    
    return outputDict

def reweightDecayPositions(LLPs, LLPchildren, lifetime, LLPid): 
    print(f"Reweighting the LLP ({LLPid}) decay positions using lifetime {lifetime} s...")
    print(datetime.datetime.now())

    #Turn off SettingWithCopyWarning, as it is not a problem here
    pd.options.mode.chained_assignment = None

    # Reweight the Decay vertex of the LLP based on a given lifetime
    LLPs['decayVertexDist_weighted'] = LLPs.apply(h.reweightDecayDistByLifetime, args=(lifetime,), axis=1)
    LLPs['decayVertex_weighted'] = LLPs.apply(h.getReweightedDecayPosition, args=("decayVertexDist_weighted",), axis=1)

    # Get the difference between the original decay position and the weighted one 
    #   - Used to translate LLP decay products positions by the same amount
    #   - The last element in this FourVector represents time, so this is preserved
    LLPs['decayVertex_translation'] = [ (b[0]-a[0], b[1]-a[1], b[2]-a[2], 0) for a, b in LLPs[["decayVertex","decayVertex_weighted"]].to_numpy()]

    # Additional reweighting methods for debugging purposes - should be similar to above
    LLPs['decayVertexDist_posWeighted'] = LLPs.apply(h.reweightDecayDistByPosition, args=(lifetime,), axis=1)
    LLPs['decayVertex_posWeighted'] = LLPs.apply(h.getReweightedDecayPosition, args=("decayVertexDist_posWeighted",), axis=1)
    
    LLPs['decayVertexDist_restWeighted'] = LLPs.apply(h.reweightDecayDistByRestLifetime, args=(lifetime,), axis=1)
    LLPs['decayVertex_restWeighted'] = LLPs.apply(h.getReweightedDecayPosition, args=("decayVertexDist_restWeighted",), axis=1)

    # Translate the decay products of the LLPs by the same distance the LLPs were. 
    LLPchildren['prodVertex_weighted'] = LLPchildren.apply(h.getDecayVertexTranslation, args=(LLPs,"prodVertex",), axis=1)
    LLPchildren['decayVertex_weighted'] = LLPchildren.apply(h.getDecayVertexTranslation, args=(LLPs,"decayVertex",), axis=1)

    # Turning back on the SettingWithCopyWarning
    pd.options.mode.chained_assignment = 'warn'

def createJetDF(eventNumbers, chargedFinalStates, neutralFinalStates):
    print("Creating jets...")
    print(datetime.datetime.now())
    # Create a set of Jets to be added to the DataFrame based on the neutral and charged final state particles
    jetDef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
    jetDict = {"eventNumber": [], "p": [], "pt": [], "px": [], "py": [], "pz": [], "E": [], "eta": [], "phi": [], "weight": []} 
    for eventNumber in eventNumbers:
        cFS = chargedFinalStates[chargedFinalStates.eventNumber == eventNumber]
        nFS = neutralFinalStates[neutralFinalStates.eventNumber == eventNumber]

        if len(cFS["prodVertex"])!=0:
            charged_PV_x, charged_PV_y, charged_PV_z, charged_PV_t = zip(*cFS["prodVertex"])
            chargedJetList = list(zip(cFS["px"], cFS["py"], cFS["pz"], cFS["E"], charged_PV_t, charged_PV_x, charged_PV_y,charged_PV_z, cFS["pt"]))
        else:
            chargedKetList = []

        if len(nFS["prodVertex"])!=0:
            neutral_PV_x, neutral_PV_y, neutral_PV_z, neutral_PV_t = zip(*nFS["prodVertex"])
            neutralJetList = list(zip(nFS["px"], nFS["py"], nFS["pz"], nFS["E"], neutral_PV_t, neutral_PV_x, neutral_PV_y,neutral_PV_z, nFS["pt"]))
        else:
            neutralJetList = []

        jetList = chargedJetList + neutralJetList
        if len(jetList) == 0:
            continue

        jetsArray = np.array(jetList,
                    dtype=[('px', float), ('py', float), ('pz', float), ('E', float), 
                           ('decay_pos_jet_t', float), ('decay_pos_jet_x', float), ('decay_pos_jet_y', float), ('decay_pos_jet_z', float),("pt",float)])

        jetsArray = ak.from_numpy(jetsArray)        
        sequence = fastjet._pyjet.AwkwardClusterSequence(jetsArray, jetDef)
        jets = sequence.inclusive_jets()

        weight = list(set(cFS["weight"]+nFS["weight"])) # Both should have the same event weight
        for i, jet in enumerate(jets):
            jetDict["eventNumber"].append(eventNumber)
            thisjet_spherical = list(h.to_spherical(jet.px, jet.py, jet.pz)) # returns r, eta, phi
            jetDict["eta"].append(thisjet_spherical[1])
            jetDict["phi"].append(thisjet_spherical[2])
            jetDict["pt"].append(np.sqrt((jet.px)**2 + (jet.py)**2))
            jetDict["p"].append(np.sqrt((jet.px)**2 + (jet.py)**2 + (jet.pz)**2))
            jetDict["px"].append(jet.px)
            jetDict["py"].append(jet.py)
            jetDict["pz"].append(jet.pz)
            jetDict["E"].append(jet.E)
            jetDict["weight"].append(weight[0])

    finalStatePromptJets =  pd.DataFrame.from_dict(jetDict)

    return finalStatePromptJets

# Now Have:
#   LLPs dataframe - containing the info on the LLPs in the sample
#   neutralFinalStates and chargedFinalStates dataframes - containing the info on the prompt non-LLP neutral and charged final states
#   finalStatePromptJets dataframe - info on the prompt jets compiled from the above final states.

#=====================================#
#=             Selection              #
#=====================================#

def applySelection(SDFs, selection={"nStations": 2, 
                                    "eachTrack": False, 
                                    "minMET": 30, 
                                    "minPt": {'LLP': 10, 'chargedTrack': 10, 'neutralTrack': 10, "jet": 10},
                                    "minP": {'LLP': 10, 'chargedTrack': 10, 'neutralTrack': 10, "jet": 10},
                                    "minDR": {'jet': 0.4, 'chargedTrack': 0.4, 'neutralTrack': 0.4},
                                    }
                                ):
    print(f"Applying Selection...")
    print(datetime.datetime.now())
    cutFlow = {}
    cutIndices = {}

    LLPs = SDFs["LLPs"]
    LLPchildren = SDFs["LLPchildren"]
    # Record number of signal events before selection
    cutFlow["nLLP_original"] = len(LLPs.index)
    cutFlow["nLLP_original_weighted"] = LLPs["weight"].sum()

    # Require LLPs that Decay 
    LLPsSel = LLPs[LLPs.nChildren > 1]
    cutFlow["nLLP_LLPdecay"] = len(LLPsSel.index)
    cutFlow["nLLP_LLPdecay_weighted"] = LLPsSel["weight"].sum()
    cutIndices["nLLP_LLPdecay"] = LLPsSel.index

    #setup cavern for selection


    cav =g.ATLASCavern()
    
    pickleDir = "./pickles"
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
    # Geometry Selection:
    #   For the decay position to be within ANUBIS use a function in which the boundaries of ANUBIS have been defined.
    print("++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++")
    print("++++++++++++++++++++++++++++++++")
    #h.precutcavern(LLPsSel)
    print("this is llp seleciton")
    print(LLPsSel)
    LLPsGeo, boundaries = h.performGeometryCut(LLPsSel,cav, ANUBISstations,minStationRadius,RPC_Pos1, RPC_Pos2, RPC_Pos3, plot = True) 
    print(LLPsGeo)
    cutFlow["nLLP_Geometry"] = len(LLPsGeo.index)
    cutFlow["nLLP_Geometry_weighted"] = LLPsGeo["weight"].sum()

    # Require two hits in ANUBIS tracking stations from LLP decay products (where a 'hit' is that a charged track intersects the tracking stations).
    
    LLPsDecays = h.checkDecayHits(LLPsGeo, cav, ANUBISstations,minStationRadius,RPC_Pos1, RPC_Pos2, RPC_Pos3,nHits = selection["nStations"],  plot = True, each=selection["eachTrack"]) # each -> each decay product has to give nHits.
    cutFlow["nLLP_Tracker"] = len(LLPsDecays.index)
    cutFlow["nLLP_Tracker_weighted"] = LLPsDecays["weight"].sum()
    cutIndices["nLLP_Tracker"] = LLPsDecays.index

    # Perform a MET selection - LLPs need to have a minimum associated MET
    LLPsMET = LLPsDecays[LLPsDecays.MET > selection["minMET"]]
    cutFlow["nLLP_MET"] = len(LLPsMET.index)
    cutFlow["nLLP_MET_weighted"] = LLPsMET["weight"].sum()
    cutIndices["nLLP_MET"] = LLPsMET.index

    # Next isolation cuts are applied, these are done on a per event basis
    #   - Delta R requirement enforced between LLP and Jet, if there is no Jet then between LLP and charged particle tracks
    passEventIdx={"both": [], "jet": [], "charged": []}
    for eventNumber in list(set(LLPsMET["eventNumber"])):
        eventJets = SDFs["finalStatePromptJets"][SDFs["finalStatePromptJets"].eventNumber == eventNumber]
        eventCharged = SDFs["chargedFinalStates"][SDFs["chargedFinalStates"].eventNumber == eventNumber]
        eventLLPs = LLPs[LLPs.eventNumber == eventNumber] 

        #NOTE: This is assuming that there is one LLP per event
        #      With multiple LLPs in an event each LLP event would need to pass the isolation requirement
        jetPassed=False
        chargedPassed=False
        if len(eventJets)!=0:
            # Require a minimum jet momenta and pt
            eventJets = eventJets[(eventJets.pt > selection["minPt"]["jet"]) & (eventJets.p > selection["minP"]["jet"])]

            eventJets["deltaEta"] = eventLLPs["eta"] - eventJets["eta"] 
            eventJets["deltaPhi"] = eventLLPs["phi"] - eventJets["phi"] 
            eventJets["deltaR"] = np.sqrt(np.power(eventJets["deltaEta"], 2) + np.power(eventJets["deltaPhi"], 2)) 
            # If there are no jets within the minimum delta R then pass the event
            if len(eventJets[eventJets.deltaR < selection["minDR"]["jet"]] == 0): 
                jetPassed = True
                passEventIdx["jet"].append(eventNumber)
        
        if len(eventCharged) != 0:
            # Require a minimum pt for charged track
            eventCharged = eventCharged[eventCharged.pt > selection["minPt"]["chargedTrack"]]

            eventCharged["deltaEta"] = eventLLPs["eta"] - eventCharged["eta"] 
            eventCharged["deltaPhi"] = eventLLPs["phi"] - eventCharged["phi"] 
            eventCharged["deltaR"] = np.sqrt(np.power(eventCharged["deltaEta"], 2) + np.power(eventCharged["deltaPhi"], 2)) 
            # If there are no charged tracks within the minimum delta R then pass the event
            if len(eventCharged[eventCharged.deltaR < selection["minDR"]["chargedTrack"]] == 0): 
                chargedPassed = True
                passEventIdx["charged"].append(eventNumber)

        if jetPassed and chargedPassed:
            passEventIdx["both"].append(eventNumber)
    
    # Get the isolated cases for just jets and just charged particles
    LLPsIsoJets = LLPsMET[LLPsMET.eventNumber.isin(passEventIdx["jet"])] 
    cutFlow["nLLP_IsoJet"] = len(LLPsIsoJets.index)
    cutFlow["nLLP_IsoJet_weighted"] = LLPsIsoJets["weight"].sum()
    cutIndices["nLLP_isoJet"] = LLPsIsoJets.index

    LLPsIsoCharged = LLPsMET[LLPsMET.eventNumber.isin(passEventIdx["charged"])] 
    cutFlow["nLLP_IsoCharged"] = len(LLPsIsoCharged.index)
    cutFlow["nLLP_IsoCharged_weighted"] = LLPsIsoCharged["weight"].sum()
    cutIndices["nLLP_isoCharged"] = LLPsIsoCharged.index

    LLPsIsoAll = LLPsMET[LLPsMET.eventNumber.isin(passEventIdx["both"])] 
    cutFlow["nLLP_IsoAll"] = len(LLPsIsoAll.index)
    cutFlow["nLLP_IsoAll_weighted"] = LLPsIsoAll["weight"].sum()
    cutIndices["nLLP_isoCharged"] = LLPsIsoAll.index

    outputDict = {"cutFlow": cutFlow,
                  "cutIndices": cutIndices,
                  "finalDF": LLPsIsoAll,
                  }

    return outputDict 


def createCTauFile(dataFrame, runID, infoDict, suffix="", outputDir=""):
    outDir=f"{outputDir}/ctaus_HNL/{runID}/"
    #outDir=f"{outputDir}/ctaus_DS/{runID}/"
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    #with open(f"{outDir}/HEPMC_DS_ctau_{infoDict[f'{runID}']['LLPmass']}.txt","w") as f:
    with open(f"{outDir}/HEPMC_HNL_ctau_{infoDict[f'{runID}']['LLPmass']}.txt","w") as f:
        f.write("mass,width,Average DecayDist/mm,Average ctau/mm\n")
        f.write(f"{infoDict[f'{runID}']['LLPmass']},{infoDict[f'{runID}']['LLPdecayWidth']},{dataFrame['decayVertexDist'].mean()},{dataFrame['ctau'].mean()}\n")



if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--coupling", nargs="+", type=float, default = [1.0])
    parser.add_argument("-m","--mass", nargs="+", type=float, default = [1.0])
    parser.add_argument("-s","--sample", nargs="+", type=str, default = ["CCDY_eev"])
    parser.add_argument("--runID", nargs="+", type=float, default = [])
    #parser.add_argument("--csv", type=str, default="", help="xsv file that summarises information on the simulated samples")
    parser.add_argument("--outDir", type=str, default = "selectedPlots") 
    parser.add_argument("--suffix", type=str, default = "") 
    parser.add_argument("--remake", action="store_true", default = False) 
    parser.add_argument("--config", type=str, default = "configs/config_HNL_m1_c1.yaml", #should change this and config file name to DS at some point 
                        help = "Configuration file containing necessary things for the run", required=True )
    args = parser.parse_args()

    # Create outDir if it doesn't exist
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)
    outputDir = args.outDir

    # Process the Configuration file
    config = h.parseConfig(args.config)

    # Process the Simulation Information table
    sampleInfoDicts = h.parseCSV(config["scanfile"])
    runIDsByMassCoupling = h.groupMassAndCoupling(sampleInfoDicts) # var[sample][mass][coupling] = list of associated runIDs
    
    # Hard-coding some inputs for testing
    runID = '0'
    sampleInfo = sampleInfoDicts[runID]
    
    #===========================================================================================================================================================
    
    #was commented out for DS 
    # Get Expected HNL Lifetime:
    HNLCalc_file = "/usera/amullin/Documents/ANUBIS/acceptance_results/HNLCalc-main/HNLCalc_ctau_table.csv"
    HNLCalc_file = "/usera/dp728/HNLCalc_ctau_table.csv"
    HNLCalc_df = pd.read_csv(HNLCalc_file)
    HNLCalc_ctau = HNLCalc_df.loc[(HNLCalc_df["Mass"] == sampleInfo["LLPmass"]) & (HNLCalc_df["Coupling"] == sampleInfo["coupling"]), "Lifetime"].values
    if HNLCalc_ctau.size > 0:
        HNLCalc_ctauVal = HNLCalc_ctau[0]  # Extract the first value (if multiple rows match)
        HNLCalc_lifetime = float(HNLCalc_ctauVal / constants.c)
        print(f"The lifetime value is: {HNLCalc_lifetime}, ctau = {HNLCalc_ctauVal}")
    else:
        print("No matching lifetime found for the given mass and coupling.")
    
   #===========================================================================================================================================================
  
    #fileName = f"/r04/atlas/amullin/ANUBIS/Simulations/SimulationSwarms/{sampleInfo['SimSwarmName']}/Events/{sampleInfo['runName']}_decayed_1/tag_1_pythia8_events.hepmc"
    fileName = f"{config['SimSwarmDir']}/{sampleInfo['SimSwarmName']}/Events/{sampleInfo['runName']}_decayed_1/tag_1_pythia8_events.hepmc"

    if os.path.exists(fileName):
        print("sample hepmc file exists.")
    elif os.path.exists(fileName+".gz"):
        print("sample hepmc file exists. unzip now.")
        file_to_unzip = fileName+".gz"
        h.unzip_file(file_to_unzip, fileName)
    else:
        print("sample hepmc NOT FOUND. simulation missing.")  
   
    #===========================================================================================================================================================

    LLPpid = 9900012
    minPt = {"chargedParticle": 0.5, "jet": 10}
    minTheta=0.828849
    maxTheta=1.51461
    
    print(config.keys())
    selection=config["selection"]
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)

    #===========================================================================================================================================================
    #   Program Flow:
    #       + If coupling of 1:
    #           - Process the relevant HEPMC file(s) and create a pre-selection h5 file or pickle of sample dataframes
    #           - Apply the selection to the sample dataframes and save a set of indices per sample dataframe for each stage of the cutflow
    #             ( This would mean that the initial saved h5 can be reused and easily get the values at each stage of the selection )
    #           - Save a dictionary of the number of remaining events at each stage of the cutflow
    #       + If coupling is not 1:
    #           - Check if the equivalent mass h5 file with a coupling of 1 has been generated, if not then create it.
    #           - Load the coupling of 1 h5 file
    #           - Reweight the decayVertex position based on the given lifetime for that coupling and mass: save a h5 file for this new dataframe
    #              (Could potentially just save the relevant columns if these h5 files are large).
    #           - Apply the selection and save the cutflow dictionary as before
    #
    #       + If plotting is specified then create a set of distributions based on the dataframes, before and after selection.
    #===========================================================================================================================================================

    start = datetime.datetime.now()
    print(f"Processing: {fileName}")
    print(start)
    parquetName = f"./tempDF_{runID}.parquet"
    if not os.path.isfile(parquetName) or args.remake:
        df = createDFfromHEPMC(fileName)
        df.to_parquet(parquetName)
    else:
        print(f"Loading dataframe from {parquetName}...")
        df = pd.read_parquet(parquetName)

    # Get a list of event numbers:
    eventNos = list(set(df['eventNumber']))


    print(f"Average general Decay Distance (mm): {df['decayVertexDist'].mean()}")
    print(f"Average general ctau (mm): {df['ctau'].mean()}")
    print(sampleInfo["LLPid"])
    #seperates it into the outputs but crucially i think it removes the chain of LLP's but it removes the later LLps which might have the longer decay distance
    sampleDFs = createSampleDataFrames(df, sampleInfo["LLPid"], minPt)

    print(f"max LLP Decay Distance (mm): {sampleDFs['LLPs']['decayVertexDist'].max()}")
    print(f"max LLP ctau (mm): {sampleDFs['LLPs']['ctau'].max()}")
    
    #print(f"max LLP Decay Distance (mm) USING LAST LLP: {sampleDFs['LLPSfordistancecalc']['decayVertexDist'].max()}")
    #print(f"max LLP ctau (mm) USING LAST LLP: {sampleDFs['LLPSfordistancecalc']['ctau'].max()}")


    #df.to_csv("./tempDF-cut?.csv")

    reweightDecayPositions(sampleDFs["LLPs"], sampleDFs["LLPchildren"], HNLCalc_lifetime, sampleInfo["LLPid"]) 
    sampleDFs["finalStatePromptJets"] = createJetDF(eventNos, sampleDFs["chargedFinalStates"], sampleDFs["neutralFinalStates"])
    print(sampleDFs)
    cutDict =  applySelection(sampleDFs, selection=selection)
    cutDict = sampleDFs
    createCTauFile(sampleDFs["LLPs"], runID, sampleInfoDicts, suffix="", outputDir="./")
    print(f"LLP Mass: {np.mean(sampleDFs['LLPs']['mass'])}")
    print(f"Average LLP Decay Distance (mm): {sampleDFs['LLPs']['decayVertexDist'].mean()}")
    print(f"Average LLP ctau (mm): {sampleDFs['LLPs']['ctau'].mean()}")


    print(f"max LLP Decay Distance (mm): {sampleDFs['LLPs']['decayVertexDist'].max()}")
    print(f"max LLP ctau (mm): {sampleDFs['LLPs']['ctau'].max()}")

    #LLPDecayDistance = myplot.plotColumnAsHist(sampleDFs["LLPSfordistancecalc"].dropna(subset=["decayVertexDist"]), "decayVertexDist", nBins = 20, normalise=False, logScale=False, outputDir=outputDir) 
    #LLPctau = myplot.plotColumnAsHist(sampleDFs["LLPSfordistancecalc"].dropna(subset=["ctau"]), "ctau", nBins = 20, normalise=False,  expfit=True, logScale=False, outputDir=outputDir) 
    
    LLPDecayDistance = myplot.plotColumnAsHist(sampleDFs["LLPs"].dropna(subset=["decayVertexDist"]), "decayVertexDist", nBins = 20, normalise=False, logScale=False, outputDir=outputDir) 
    LLPctau = myplot.plotColumnAsHist(sampleDFs["LLPs"].dropna(subset=["ctau"]), "ctau", nBins = 20, normalise=False,  expfit=True, logScale=False, outputDir=outputDir) 
   
    
    METxplot = myplot.plotColumnAsHist(sampleDFs["LLPs"], "METx", nBins = 20, normalise=True, logScale=False, outputDir=outputDir)
    METyplot = myplot.plotColumnAsHist(sampleDFs["LLPs"], "METy", nBins = 20, normalise=True, logScale=False, outputDir=outputDir)
    METplot = myplot.plotColumnAsHist(sampleDFs["LLPs"], "MET", nBins = 20, axisRange = (0,200), metcutoff = True, normalise=True, logScale=False, outputDir=outputDir)
    #ßmyplot.plotCutFlow(cutDict["cutFlow"], logScale = True, outputDir=outputDir, suffix=f"_{runID}")

    # Debug MET by getting the pt of the neutrinos and LLPs directly
    #one apparently is sum of missing stuff like neutrinos
    #the other is all the non missing stuff minus total?
    testMET = []
    for event in list(set(sampleDFs["LLPs"]["eventNumber"])):
        METx = sampleDFs["LLPs"][(sampleDFs["LLPs"].eventNumber == event)]["px"].sum(min_count=1)
        tempMETx= sampleDFs["finalStates_Neutrinos"][(sampleDFs["finalStates_Neutrinos"].eventNumber == event)]["px"].sum(min_count=1)
        if not np.isnan(tempMETx): #There may not be any neutrinos in the final state
            METx += tempMETx
        METy = sampleDFs["LLPs"][(sampleDFs["LLPs"].eventNumber == event)]["py"].sum(min_count=1)
        tempMETy= sampleDFs["finalStates_Neutrinos"][(sampleDFs["finalStates_Neutrinos"].eventNumber == event)]["py"].sum(min_count=1)
        if not np.isnan(tempMETy): #There may not be any neutrinos in the final state
            METy += tempMETy

        if np.isnan(METx) and np.isnan(METy):
            continue

        if np.isnan(METx):
            METx = 0

        if np.isnan(METy):
            METy = 0

        MET = h.calculate_pt(METx, METy)
        testMET.append(MET)



   #================================================================-
    """
    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlabel("Decay Distance (mm)", fontsize = 20)
    ax.set_ylabel(f'Counts / Bin', fontsize=20)
    ax.hist(LLPDecayDistance, density = False, bins = 16, range=(-10,150), color = "coral", ec = "red", fill=False, label="DirectMET")
    plt.grid()              
    plt.legend(loc="best")
    plotName=outputDir+f"/LLPNG_decayVertexDist{args.suffix}.png"
    #plt.show()
    plt.savefig(plotName)
    plt.close(fig)
   #================================================================-

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlabel("ctau (mm)", fontsize = 20)
    ax.set_ylabel(f'Counts / Bin', fontsize=20)
    ax.hist(LLPctau, density = False, bins = 16, range=(-10,150), color = "coral", ec = "red", fill=False, label="DirectMET")

    import scipy as sp
    def exp_decay(x, a, b):
        return a * np.exp(-b * x)
    counts, bin_edges = np.histogram(LLPctau, bins=16, range=(-10,150), density=False)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Compute bin centers
    valid = bin_centers > 0
    popt, pcov = sp.optimize.curve_fit(exp_decay, bin_centers, counts, p0=(max(counts), 0.1))

    x_fit = np.linspace(min(bin_centers), max(bin_centers), 100)
    y_fit = exp_decay(x_fit, *popt)
    plt.plot(x_fit, y_fit, 'r-', linewidth=2, label="Exponential Fit")
 
    plt.grid()              
    plt.legend(loc="best")
    plotName=outputDir+f"/LLPNG_ctau{args.suffix}.png"
    #plt.show()
    plt.savefig(plotName)
    plt.close(fig)
    """
    #================================================================-

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlabel("Direct MET", fontsize = 20)
    ax.set_ylabel(f'Counts / Bin', fontsize=20)
    ax.hist(testMET, density = False, bins = 16, range=(-10,150), color = "coral", ec = "red", fill=False, label="DirectMET")
    ax.hist(sampleDFs["LLPs"]["MET"].to_numpy(), range=(-10,150), density = False, bins = 16, color = "turquoise", ec = "teal", fill=False, label="MET")
    """
    import scipy as sp
    def exp_decay(x, a, b):
        return a * np.exp(-b * x)
    counts, bin_edges = np.histogram(testMET, bins=16, range=(-10,150), density=False)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Compute bin centers
    valid = bin_centers > 0
    popt, pcov = sp.optimize.curve_fit(exp_decay, bin_centers, counts, p0=(max(counts), 0.1))

    x_fit = np.linspace(min(bin_centers), max(bin_centers), 100)
    y_fit = exp_decay(x_fit, *popt)
    plt.plot(x_fit, y_fit, 'r-', linewidth=2, label="Exponential Fit")
    """
    plt.grid()              
    plt.legend(loc="best")
    print(outputDir)
    plotName=outputDir+f"/LLPNG_MET_FromLLPsAndNeutrinos{args.suffix}.png"
    #plt.show()
    plt.savefig(plotName)
    plt.savefig("test")
    plt.close(fig)

    #===========================================================================================================================================================

                



    """
    import csv

    csv_out_dir = f"output/SS_DS"

    if not os.path.exists(csv_out_dir):
        os.makedirs(csv_out_dir)
        print(f"Directory {csv_out_dir} created.")

    filedescription = f"HS_mn1_c1_{sampleInfoDicts[runID]}"

      
    with open(f'{csv_out_dir}/cuts_dict_{filedescription}.pkl', 'wb') as fp:
        pickle.dump(cutDict, fp)
        print(f'dictionary saved to pkl file for {filedescription}')

    with open(f'{csv_out_dir}/cuts_dict_{filedescription}.csv', 'w') as csvfile:
        csvwriter = csv.DictWriter(csvfile, fieldnames=cutDict)#, delimiter=' ') #its definitely not cutDict
        csvwriter.writeheader()
        csvwriter.writerow(cutDict)
        print(f'dictionary saved to csv file for {filedescription}')
"""

#can save cutflows in some format like json so dont have to rerun over and over
#
#
#
#
#
#
#
#
#
#
    end = datetime.datetime.now()
    print(f"total time taken: {end-start}")

    #event number visualiser



    # Helper function to project track to the ceiling
    def projectToCeiling(start_vertex, end_vertex, y_ceiling):
        """
        Given two points (start_vertex, end_vertex) in 3D, project the track to the given z_ceiling.
        """
        x1, y1 = start_vertex
        x2, y2 = end_vertex
        
        # Assuming the track follows a straight line, we need to find the intersection of the line with z = z_ceiling
     
        if y2 != y1 and y_ceiling > y1:  # Avoid division by zero
            t = (y_ceiling - y1) / (y2 - y1)
            projected_x = x1 + t * (x2 - x1)
            return (projected_x, y_ceiling)
        else:
            return (x2, y2)  # Parallel to ceiling


        
    
    def plot_event(event_number, LLPs, LLPchildren, finalStates):
        """
        Visualize an event's LLP and its decay products' tracks in 3D.

        Parameters:
        - event_number: The event number to visualize.
        - LLPs: DataFrame containing information about the LLPs.
        - LLPchildren: DataFrame containing information about the decay products of LLPs.
        - finalStates: DataFrame of all final-state particles.

        This function assumes the following columns in the input DataFrames:
        - 'prodVertex', 'decayVertex', 'px', 'py', 'pz' (for momentum components),
        'pt', 'charge', etc.
        """
        
        # Filter data for the given event number
        event_LLPs = LLPs[LLPs['eventNumber'] == event_number]
        event_children = LLPchildren[LLPchildren['eventNumber'] == event_number]
        
        if event_LLPs.empty:
            print(f"No LLPs found for event {event_number}.")
            return
        #set up cavern
        
        # Create a 3D plot
            #setup cavern for selection


        cav =g.ATLASCavern()
        
        pickleDir = "./pickles"
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
        plotRPCs={"xy": True, "xz": False, "zy": False, "3D": False}
        ANUBISrpcs=[RPC_Pos1, RPC_Pos2, RPC_Pos3]
        suffix ="idk"
        #XY
        fig, ax = plt.subplots(1, figsize=(10, 16), dpi=100)
                # ANUBISrpcs is a list of RPClayers, each contain a list of RPCs in the format:
        #   - "corners": 8 (x,y,z) coordinates corresponding to their corners,
        #   - "midPoint": The midPoint of the RPC in (x,y,z), 
        #   - "LayerID" and "RPCid": A Layer ID and RPC ID to uniquely identify the RPC
        #   - "plane": A Sympy plane in the eta-phi plane that passes through the midpoint

        # Get the Cavern ceiling data grid
        cavernArch = cav.createCavernVault(doPlot=False)
        # Get the Access Shafts data grid 
        accessShafts = cav.createAccessShafts()
        # Cavern Boundaries
        cavernBounds = { "x": np.linspace(cav.CavernX[0], cav.CavernX[1],100),
                         "y": np.linspace(cav.CavernY[0], cav.CavernY[1],100),
                         "z": np.linspace(cav.CavernZ[0], cav.CavernZ[1],100),}
        tempSuffix= { "xy": suffix, "xz": suffix, "zy": suffix, "3D": suffix}

        if len(ANUBISrpcs)!=0:
            LayerColours = ["violet", "coral", "green", "maroon", "chartreuse", "gray"]
            
            for label, decision in plotRPCs.items():
                if decision:
                    tempSuffix[label]+=f"_withANUBIS_{len(ANUBISrpcs)}Layers"


        # Cavern Boundaries
        ax.plot([cav.CavernX[0], cav.CavernX[0]], [cav.CavernY[0], cav.CavernY[1]], c="paleturquoise")
        ax.plot([cav.CavernX[1], cav.CavernX[1]], [cav.CavernY[0], cav.CavernY[1]], c="paleturquoise")
        ax.plot([cav.CavernX[0], cav.CavernX[1]], [cav.CavernY[0], cav.CavernY[0]], c="paleturquoise")
        # Cavern Ceiling
        ax.plot(cavernArch["x"], cavernArch["y"], c="paleturquoise")
        # Access Shafts
        ax.plot([cav.PX14_Centre["x"]-cav.PX14_Radius,cav.PX14_Centre["x"]-cav.PX14_Radius],
                [cav.PX14_LowestY, cav.PX14_Centre["y"]+cav.PX14_Height], c="red", label="PX14", alpha=0.5)
        ax.plot([cav.PX14_Centre["x"]+cav.PX14_Radius,cav.PX14_Centre["x"]+cav.PX14_Radius],
                [cav.PX14_LowestY, cav.PX14_Centre["y"]+cav.PX14_Height], c="red", label="PX14", alpha=0.5)

        ax.plot([cav.PX16_Centre["x"]-cav.PX16_Radius,cav.PX16_Centre["x"]-cav.PX16_Radius],
                [cav.PX16_LowestY, cav.PX16_Centre["y"]+cav.PX16_Height], c="turquoise", label="PX16", alpha=0.5)
        ax.plot([cav.PX16_Centre["x"]+cav.PX16_Radius,cav.PX16_Centre["x"]+cav.PX16_Radius],
                [cav.PX16_LowestY, cav.PX16_Centre["y"]+cav.PX16_Height], c="turquoise", label="PX16", alpha=0.5)
        
        if len(ANUBISrpcs)!=0 and plotRPCs["xy"]: #Include ANUBIS positions
            nLayer=0
            for rpcLayer in ANUBISrpcs:
                tempRPCList = cav.convertRPCList(rpcLayer)
                initialZ = tempRPCList["corners"][0][0][2]
                # Plot midpoints
                ax.scatter([x[0] for x in tempRPCList["midPoint"]], [y[1] for y in tempRPCList["midPoint"]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
                # Plot RPC Boundaries
                nInterations=0
                for i in range(len(tempRPCList["corners"])):
                    c = tempRPCList["corners"][i]
                    
                    if c[0][2] != initialZ:
                        continue
                    nInterations+=1

                    # For Front of RPC in z, use 1--3, 1--5, 3--7, 5--7 for back of RPC.
                    ax.plot( [c[0][0], c[2][0]], [c[0][1],c[2][1]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
                    ax.plot( [c[0][0], c[4][0]], [c[0][1],c[4][1]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
                    ax.plot( [c[2][0], c[6][0]], [c[2][1],c[6][1]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
                    ax.plot( [c[4][0], c[6][0]], [c[4][1],c[6][1]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
                nLayer+=1
            
            # Naive envelopes
            ax.add_patch( matplotlib.patches.Arc((cav.centreOfCurvature["x"], cav.centreOfCurvature["y"]),
                                                width = 2*(cav.archRadius), height = 2*(cav.archRadius), 
                                                theta1=min(cav.angles["theta"])*(180/np.pi), theta2=max(cav.angles["theta"])*(180/np.pi),
                                                color=LayerColours[0], fill=False, ls="--") )
            ax.add_patch( matplotlib.patches.Arc((cav.centreOfCurvature["x"], cav.centreOfCurvature["y"]),
                                                width = 2*(cav.archRadius-0.06), height = 2*(cav.archRadius-0.06),
                                                theta1=min(cav.angles["theta"])*(180/np.pi), theta2=max(cav.angles["theta"])*(180/np.pi),
                                                color=LayerColours[0], fill=False, ls="--") )

            ax.add_patch( matplotlib.patches.Arc((cav.centreOfCurvature["x"], cav.centreOfCurvature["y"]),
                                                width = 2*(cav.archRadius-0.4), height = 2*(cav.archRadius-0.4),
                                                theta1=min(cav.angles["theta"])*(180/np.pi), theta2=max(cav.angles["theta"])*(180/np.pi),
                                                color=LayerColours[1], fill=False, ls="--") )
            
            ax.add_patch( matplotlib.patches.Arc((cav.centreOfCurvature["x"], cav.centreOfCurvature["y"]),
                                                width =2*(cav.archRadius-0.4-0.06), height = 2*(cav.archRadius-0.4-0.06),
                                                theta1=min(cav.angles["theta"])*(180/np.pi), theta2=max(cav.angles["theta"])*(180/np.pi),
                                                color=LayerColours[1], fill=False, ls="--") )

            ax.add_patch( matplotlib.patches.Arc((cav.centreOfCurvature["x"], cav.centreOfCurvature["y"]),width = 2*(cav.archRadius-1), height = 2*(cav.archRadius-1),
                                                theta1=min(cav.angles["theta"])*(180/np.pi), theta2=max(cav.angles["theta"])*(180/np.pi),
                                                color=LayerColours[2], fill=False, ls="--") )
            ax.add_patch( matplotlib.patches.Arc((cav.centreOfCurvature["x"], cav.centreOfCurvature["y"]),
                                                width = 2*(cav.archRadius-1-0.06), height = 2*(cav.archRadius-1-0.06),
                                                theta1=min(cav.angles["theta"])*(180/np.pi), theta2=max(cav.angles["theta"])*(180/np.pi),
                                                color=LayerColours[2], fill=False, ls="--") )

        # Plot the ATLAS experiment boundaries
        plotATLAS = True
        if plotATLAS:
            ax.add_patch( plt.Circle((cav.ATLAS_Centre["x"], cav.ATLAS_Centre["y"]), cav.radiusATLAS, color="blue", fill=False, ls="--") )

        # Mark the Cavern Centre, IP, and Centre of Curvature for the ceiling
        ax.scatter(0, 0, c="r", marker = "x", label="Cavern Centre")
        ax.annotate("Centre", (0,0))
        ax.scatter(cav.IP["x"], cav.IP["y"], c="g", marker = "o", label="IP")
        ax.annotate("IP", (cav.IP["x"], cav.IP["y"]))
        ax.scatter(cav.centreOfCurvature["x"], cav.centreOfCurvature["y"], c="b", marker = "D", label="Ceiling Centre of Curvature")
        ax.annotate("Centre of Curvature (Ceiling)", (cav.centreOfCurvature["x"], cav.centreOfCurvature["y"]))


        """
        plt.savefig(f"{outDir}/ATLASCavern_XY_withShafts{tempSuffix['xy']}.pdf")
        plt.gca().set_aspect('equal')
        ax.set_ylim(top=24)
        plt.savefig(f"{outDir}/ATLASCavern_XY{tempSuffix['xy']}.pdf")
        ax.set_xlim(-16,-5)
        ax.set_ylim(10,20)
        plt.savefig(f"{outDir}/ATLASCavern_XY_Zoomed{tempSuffix['xy']}.pdf")
        plt.close()
        """
        # Plot LLPs and their decay vertices
        for _, LLP in event_LLPs.iterrows():
            X, Y = LLP['prodVertex'][0], LLP['prodVertex'][1]
            X = X-1.7 
            Y = Y-cav.CavernYLength/2 + 11.37

            X2, Y2 = LLP['decayVertex'][0], LLP['decayVertex'][1]
            X2 = X2-1.7 
            Y2 = Y2-cav.CavernYLength/2 + 11.37

            # Plot the LLP production vertex (assumed to be at (prodVertex))
            ax.scatter(X,Y,  color='blue', label='LLP Production Vertex')

            # Plot the LLP decay vertex
            ax.scatter(X2,Y2, color='red', label='LLP Decay Vertex')
            production_vertex = [X,Y]
            decay_vertex= [X2,Y2]
            y_ceiling = cav.CavernY[1]  # or use the ceiling's z-coordinate if it's defined elsewhere
            projection_to_ceiling = projectToCeiling(production_vertex, decay_vertex, y_ceiling)
            
            # Plot the projected track to ceiling
            print(f"Production: {production_vertex}, Projection: {projection_to_ceiling}")

            ax.plot([production_vertex[0], projection_to_ceiling[0]], [production_vertex[1], projection_to_ceiling[1]], c="blue", linestyle="--")

        # Plot decay products' tracks
        for _, decay in event_children.iterrows():
            X, Y = decay['prodVertex'][0], decay['prodVertex'][1]
            X = X-1.7 
            Y = Y-cav.CavernYLength/2 + 11.37

            X2, Y2 = decay['decayVertex'][0], decay['decayVertex'][1]
            X2 = X2-1.7 
            Y2 = Y2-cav.CavernYLength/2 + 11.37
        
            # Plot the decay product production vertex
            ax.scatter(X, Y, color='green', label='Decay Product Vertex')

            # Project the decay product track
            child_vertex =[X, Y]
            child_decay_vertex = [X2, Y2]

            child_projection = projectToCeiling(child_vertex, child_decay_vertex, y_ceiling)
            print(f"Production: {child_vertex}, Projection: {child_decay_vertex}")


            # Project the decay product's track to the ceiling
            ax.plot([child_vertex[0], child_projection[0]], [child_vertex[1], child_projection[1]], c="green", label="Child Projected Track", linestyle="-.")

        
        # Set title
        ax.set_title(f"Event {event_number}: LLP and Decay Products")
        plt.xlabel(f"x /m")
        plt.ylabel(f"y /m")
        plt.title("ATLAS Cavern")
        # Show the plot
        plt.savefig('/usera/dp728/event_plot.png', dpi=300)
        plt.show()

    def select_event(LLPs, LLPchildren, finalStates):
        """
        Let the user select an event to visualize.
        """

        # Get a list of unique event numbers from the data
        event_numbers = sorted(set(LLPs['eventNumber']))
        
        print(f"Available event numbers: {event_numbers}")
        
        # Prompt the user to select an event number
        x = 0
        while x == 0:
        
            event_number = int(input(f"Enter an event number from the available list: "))
            
            if event_number in event_numbers:
                print(f"Visualizing event {event_number}...")
                plot_event(event_number, LLPs, LLPchildren, finalStates)

            else:
                print(f"Invalid event number. Please choose from the available events: {event_numbers}")
        


    # Example of how to call this function:
    # Assuming 'LLPs', 'LLPchildren', and 'finalStates' are your dataframes
    select_event(sampleDFs['LLPs'], sampleDFs["LLPchildren"], sampleDFs["finalStates"])
