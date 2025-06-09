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
        
            """

            event_iter = iter(f)  # Convert to iterator explicitly
            while True:
                try:
                    event = next(event_iter)
                except StopIteration:
                    break
                except Exception as e:
                    print(f"[WARN] Skipping malformed event {eventNumber}: {e}")
                    skippedEvents += 1
                    eventNumber += 1
                    continue
            """
        
            #eventNumber = event.event_number # The HEPMC Event number appears to reset after 124 events, and would no longer be a unique identifier.
            if eventNumber%100 == 0:
                print(f"Event {eventNumber}: time elapsed {datetime.datetime.now()-start}")
        
            for p in event.particles:
                #print(p)
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
                pid = p.pid
                hepMCdict["PID"].append(pid)

                # Useful list of IDs: http://home.thep.lu.se/~bierlich/misc/keyword-manual/ParticleData.html
                if abs(pid) in [35, 5214, 5212, 5322, 5324, 203122, 9940003, 9900012, 9940011, 9940103, 9940023, 9941003, 
                                9941103, 9942003, 9942103, 9942033, 9950003, 9950005, 9951003, 9951103, 9951203]: 
                    # Special cases for known neutral particles 
                    charge=0
                elif abs(pid) in [543, 4124, 5314, 20413, 9932103]:
                    # Special cases for known charged particles 
                    charge=-1 if pid < 0 else 1
                elif abs(pid) in [4424]:
                    # Special cases for known doubly charged particles 
                    charge=-2 if pid < 0 else 2
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
                """
                except Exception as e:
                    print(f"[ERROR] Failed to process event {eventNumber}: {e}")
                    skippedEvents += 1
                """

            
            eventNumber+=1

    print(f"Unknown PIDs: {list(set(unknownPIDs))}")
    end = datetime.datetime.now()
    print(f"total time taken: {end-start}")

    df = pd.DataFrame.from_dict(hepMCdict)
    # Ensure the MET branches are treated as floats
    df["MET"] = df["MET"].astype(float)
    df["METx"] = df["METx"].astype(float)
    df["METy"] = df["METy"].astype(float)

    #dfDir = "./dataframes"
    #if not os.path.exists(dfDir):
    #    os.makedirs(dfDir)

    #df.to_csv(f"{dfDir}/tempDF.csv")

    return df, unknownPIDs





def doPhiFold(row, foldEvents, useWeighted=False):
    prodVertexName = "prodVertex"
    decayVertexName = "decayVertex"
    if useWeighted:
        prodVertexName+="_weighted"
        decayVertexName+="_weighted"

    if row["eventNumber"] in foldEvents:
        prodVertex = row[prodVertexName]
        decayVertex = row[decayVertexName]

        newProdVertex = np.array((-prodVertex[0], -prodVertex[1], prodVertex[2], prodVertex[3]))
        if np.all(decayVertex == (-1,-1,-1,-1)):
            newDecayVertex = decayVertex
        else:
            newDecayVertex= np.array((-decayVertex[0], -decayVertex[1], decayVertex[2], decayVertex[3]))

        return [-row["px"], -row["py"], row["phi"]+np.pi, newProdVertex, newDecayVertex]

    return [row["px"], row["py"], row["phi"], row[prodVertexName], row[decayVertexName]]





def childrenHunt(dataFrame, eventNumber, index, childList=[], LLPids=[]):
    indices = dataFrame.index[(dataFrame["eventNumber"] == eventNumber) & (dataFrame["particleIndex"] == index)].tolist()
    childrenIndexList = dataFrame.loc[indices].childrenIndices.tolist()

    if len(indices)!=1:
        # This should not be possible functionally with how these have been constructed, but this is a sanity check
        raise Exception("This dataframe cannot uniquely identify events based on eventNumber and particleIndex")

    if len(childrenIndexList)==0:
        return childList
    else:
        for childrenIndices in childrenIndexList:
            for idx in childrenIndices:
                if idx in childList:
                    continue

                # In situations where the parent LLP has one child that is also an LLP do not make this a child
                # This is because in MadGraph you can get intermediate particles where LLP -> LLP -> LLP but these are all virtually identical.
                # This selection should ensure that only the last LLP in the chain is selected -- preserving the decay information e.g. nChildren.
                parent = dataFrame.loc[indices[0]]
                childPID = dataFrame[(dataFrame["eventNumber"]==eventNumber) & (dataFrame["particleIndex"] == idx)]["PID"].iloc[0]
                if (childPID in LLPids) and (parent["PID"] in LLPids) and\
                   (childPID == parent["PID"]) and (parent["nChildren"]==1):
                    continue

                childList.append(idx)
                childrenHunt(dataFrame, eventNumber, idx, childList)
            return childList

 

def createSampleDataFrames(df, LLPid, minPt):

    # Get all the final state particles
    finalStates = df[(df.nChildren == 0) & (df.status==1)]

    # Get the LLPs with the given ID and their indices
    LLPindices = df.index[df["PID"]==LLPid].tolist()
    LLPs = df[df.PID == LLPid]

    # Get the decay particles of the LLP
    #   - To allow for selection based on their tracks
    LLPchildrenDict={"eventNumber": [], "childrenIndices": [], "LLPindex": []}
    for index in LLPindices:
        tempList=[]
        tempEvtNo = df.loc[index].eventNumber
        tempParticleIndex = df.loc[index].particleIndex
        # Children Hunt gives the list of particle IDs of the children -- resets per event.
        tempChildIndices = childrenHunt(df, tempEvtNo, tempParticleIndex, tempList, LLPids = [LLPid]) #Returns particleID list

        for chIdx in tempChildIndices:
            LLPchildrenDict["eventNumber"].append(tempEvtNo)
            LLPchildrenDict["childrenIndices"].append(chIdx)
            LLPchildrenDict["LLPindex"].append(index)

    print("After Getting LLPchildren")
    print(datetime.datetime.now())

    # The above gives the per event children indices, this gives a list of the dataframe indices for those
    childrenDFindices = []
    originatingLLP = []
    for i in range(len(LLPchildrenDict["eventNumber"])):
        tmpChildrenDFindices = df[(df.eventNumber == LLPchildrenDict["eventNumber"][i]) &\
                                  (df.particleIndex == LLPchildrenDict["childrenIndices"][i])].index.tolist()
        childrenDFindices.extend(tmpChildrenDFindices)
        originatingLLP.extend([LLPchildrenDict["LLPindex"][i]]*len(tmpChildrenDFindices))

    LLPchildren = df.iloc[childrenDFindices]

    # Add a column to identify which LLP the decay products came from originally.
    childrenDF = pd.DataFrame(originatingLLP, columns=["LLPindex"], index=childrenDFindices)

    LLPchildren["LLPindex"] = childrenDF["LLPindex"].values
    LLPchildren = LLPchildren[LLPchildren.PID != LLPid] # Remove the LLP from the decay products
    LLPchildren = LLPchildren[~LLPchildren.index.duplicated(keep='first')]
    
    # Sometimes a MC generator can propagate a particle without a decay such that the 'decay chain' is of the form:
    #   X ---> X ---> X ---> decay products
    # NOTE: The actual values for the intermediates' columns (e.g. decayVertex or momenta) is effectively identical to the original.
    #
    # The below line allows for the final LLP in the decay chain to be stored by examining the LLP children (that don't contain LLPs)
    # However, if the LLP doesn't decay then this won't be captured, so also include final state LLPs.
    LLPs = LLPs[( LLPs.index.isin(list(set(LLPchildren["LLPindex"].tolist()))) ) | (LLPs.status == 1)]

    print("After Disambiguating LLPchildren and LLPs")
    print(datetime.datetime.now())

    # Select Final states that don't originate from the LLP
    finalStates_NoLLP = finalStates[~finalStates.index.isin(childrenDFindices)] # Remove the LLP decay products
    finalStates_NoLLP = finalStates_NoLLP[finalStates_NoLLP.PID != LLPid] # Sanity check - Ensure that the LLP isn't included
    
    # Select the final State neutrinos
    finalStates_Neutrinos = finalStates_NoLLP[finalStates_NoLLP.PID.isin([12,14,16,18])]
    finalStates_NoLLP = finalStates_NoLLP[~finalStates_NoLLP.PID.isin([12,14,16,18])] # Remove the neutrinos with their PIDs 

    # Select Charged Final states that don't originate from the LLP
    #   Require it to be prompt, i.e. have a production vertex <10mm from the IP
    #   Require a minimum pT given in the config file
    chargedFinalStates = finalStates_NoLLP[(finalStates_NoLLP.charge != 0) & (finalStates_NoLLP.charge != -0.555)] 
    chargedFinalStates = chargedFinalStates[chargedFinalStates.prodVertexDist < 10]
    chargedFinalStates = chargedFinalStates[chargedFinalStates.pt > minPt["chargedTrack"]]

    # Select Neutral Final states that don't originate from the LLP
    #   Require it to be prompt, i.e. have a production vertex <10mm from the IP
    neutralFinalStates = finalStates_NoLLP[finalStates_NoLLP.charge == 0] 
    neutralFinalStates = neutralFinalStates[neutralFinalStates.prodVertexDist < 10] 

    print("After Getting finalStates")
    print(datetime.datetime.now())

    # For each Event get the associated pt of the final states and assign this to the LLP as MET
    #   - This doesn't include neutrinos, the LLPs themselves or the LLPs' decay products
    #   - We care about the total MET in an event for selection purposes, and so each LLP in an event will be assigned the total MET for ease of reference. 
    for event in list(set(df['eventNumber'])):
        METx = finalStates_NoLLP[finalStates_NoLLP.eventNumber == event]["px"].sum()
        METy = finalStates_NoLLP[finalStates_NoLLP.eventNumber == event]["py"].sum()
        if len(LLPs[LLPs.eventNumber == event])!=0:
            LLPs.loc[LLPs.loc[LLPs.eventNumber == event].index, "METx"] = METx
            LLPs.loc[LLPs.loc[LLPs.eventNumber == event].index, "METy"] = METy
            LLPs.loc[LLPs.loc[LLPs.eventNumber == event].index, "MET"] = h.calculate_pt(METx, METy)

    print("After Getting MET")
    print(datetime.datetime.now())

    outputDict={"finalStates": finalStates,
                "LLPs": LLPs,
                "LLPchildren": LLPchildren, 
                "finalStates_NoLLP": finalStates_NoLLP, # Prompt final states
                "finalStates_Neutrinos": finalStates_Neutrinos, # Prompt neutrino final states
                "chargedFinalStates": chargedFinalStates, # Charged Prompt Final States
                "neutralFinalStates": neutralFinalStates,
    }
    
    return outputDict



def reweightDecayPositions(LLPs, LLPchildren, lifetime, LLPid, seed=""): 

    #Turn off SettingWithCopyWarning, as it is not a problem here
    pd.options.mode.chained_assignment = None

    np.random.seed(seed)  # Set the random seed for reproducibility

    # Reweight the Decay vertex of the LLP based on a given lifetime
    LLPs['decayVertexDist_weighted'] = LLPs.apply(h.reweightDecayDistByLifetime, args=(lifetime,), axis=1)# getting rid of seed the seed might persist into the fucntion although not sure 
    #LLPs['decayVertexDist_weighted'] = LLPs.apply(h.reweightDecayDistByLifetime, args=(lifetime,seed), axis=1)
    LLPs['decayVertex_weighted'] = LLPs.apply(h.getReweightedDecayPosition, args=("decayVertexDist_weighted",), axis=1)
    LLPs['ctau_weighted'] = LLPs['decayVertexDist_weighted']/(LLPs['boost']*LLPs['beta']) 

    # Get the difference between the original decay position and the weighted one 
    #   - Used to translate LLP decay products positions by the same amount
    #   - The last element in this FourVector represents time, so this is preserved
    # NOTE: DOUBLE CHECK!!!!!
    LLPs['decayVertex_translation'] = [ (b[0]-a[0], b[1]-a[1], b[2]-a[2], 0) for a, b in LLPs[["decayVertex","decayVertex_weighted"]].to_numpy()]

    # Additional reweighting methods for debugging purposes - should be similar to above
    LLPs['decayVertexDist_posWeighted'] = LLPs.apply(h.reweightDecayDistByPosition, args=(lifetime,), axis=1)
    LLPs['decayVertex_posWeighted'] = LLPs.apply(h.getReweightedDecayPosition, args=("decayVertexDist_posWeighted",), axis=1)
    LLPs['ctau_posWeighted'] = LLPs['decayVertexDist_posWeighted']/(LLPs['boost']*LLPs['beta']) 
    
    LLPs['decayVertexDist_restWeighted'] = LLPs.apply(h.reweightDecayDistByRestLifetime, args=(lifetime,), axis=1)
    LLPs['decayVertex_restWeighted'] = LLPs.apply(h.getReweightedDecayPosition, args=("decayVertexDist_restWeighted",), axis=1)
    LLPs['ctau_restWeighted'] = LLPs['decayVertexDist_restWeighted']/(LLPs['boost']*LLPs['beta']) 

    # Translate the decay products of the LLPs by the same distance the LLPs were. 
    #LLPchildren['prodVertex_weighted'] = LLPchildren.apply(h.getDecayVertexTranslation, args=(LLPs,"prodVertex",), axis=1)
    if not LLPchildren.empty:
        LLPchildren["prodVertex_weighted"] = LLPchildren.apply(
            h.getDecayVertexTranslation, args=(LLPs, "prodVertex",), axis=1
        )
    else:
        # Assign placeholder tuple to all rows (none in this case, but this keeps structure consistent)
        LLPchildren["prodVertex_weighted"] = [(-2, -2, -2, -2)] * len(LLPchildren)
    
    #LLPchildren['decayVertex_weighted'] = LLPchildren.apply(h.getDecayVertexTranslation, args=(LLPs,"decayVertex",), axis=1)
    # Assign 'decayVertex_weighted' safely
    if not LLPchildren.empty:
        LLPchildren["decayVertex_weighted"] = LLPchildren.apply(
            h.getDecayVertexTranslation, args=(LLPs, "decayVertex",), axis=1
        )
    else:
        LLPchildren["decayVertex_weighted"] = [(-2, -2, -2, -2)] * len(LLPchildren)
    # Turning back on the SettingWithCopyWarning
    pd.options.mode.chained_assignment = 'warn'



def createJetDF(eventNumbers, chargedFinalStates, neutralFinalStates):
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
            chargedJetList = []

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
                                    "geometry": "",}
                                    ):
    #Turn off SettingWithCopyWarning, as it is not a problem here
    pd.options.mode.chained_assignment = None

    cutFlow, cutIndices = {}, {}

    LLPs = SDFs["LLPs"]
    # Record number of signal events before selection
    cutFlow["nLLP_original"] = len(LLPs.index)
    cutFlow["nLLP_original_weighted"] = LLPs["weight"].sum()

    # Require LLPs that Decay
    LLPsSel = LLPs[LLPs.status != 1]
    cutFlow["nLLP_LLPdecay"] = len(LLPsSel.index)
    cutFlow["nLLP_LLPdecay_weighted"] = LLPsSel["weight"].sum()
    cutIndices["nLLP_LLPdecay"] = LLPsSel.index.tolist()

    # Geometry Selection:
    #   For the decay position to be within ANUBIS use a function in which the boundaries of ANUBIS have been defined.
    """
    LLPsGeo = h.performGeometryCut(LLPsSel, selection["geometry"].RPCMaxRadius, trackingOnly=True, geometry=selection["geometry"], decayVertex="decayVertex_weighted") 
    cutFlow["nLLP_Geometry"] = len(LLPsGeo.index)
    cutFlow["nLLP_Geometry_weighted"] = LLPsGeo["weight"].sum()
    cutIndices["nLLP_Geometry"] = LLPsSel.index

    myplot.plotDecayVertexPositions(LLPsGeo, selection["geometry"], nBins=[100,100,100], outDir=outDir, suffix="AfterGeoCut", decayVertex="decayVertex_weighted")
    """

    # Check LLP's decay vertex is in the Cavern Volume and within a maxRadius of the centre of curvature of the ceiling
    LLPsInCavern = LLPsSel[LLPsSel.apply(h.checkInCavern, args=(selection["geometry"], selection["geometry"].RPCMaxRadius, "decayVertex_weighted"), axis=1)]
    cutFlow["nLLP_InCavern"] = len(LLPsInCavern.index)
    cutFlow["nLLP_InCavern_weighted"] = LLPsInCavern["weight"].sum()
    cutIndices["nLLP_InCavern"] = LLPsInCavern.index.tolist()
    
    """
    print("in cavern cut")
    passed_count = LLPsInCavern.shape[0]
    total = LLPsSel.shape[0]
    print("total")
    print(total)
    failed_count = total - passed_count
    print("passed count")

    print(passed_count)
    if not total > 0:
        total = 1


    print(f"Passed: {passed_count} ({passed_count/total:.2%})  | Failed: {failed_count} ({failed_count/total:.2%})")
    """

    # Check that the LLP decay vertex is outside of the ATLAS detector
    LLPsNotInATLAS = LLPsInCavern[~LLPsInCavern.apply(h.checkInATLAS, args=(selection["geometry"], True, "decayVertex_weighted"), axis=1)]
    
    LLPsInATLAS = LLPsInCavern[~np.isin(LLPsInCavern, LLPsNotInATLAS)] 
    
    cutFlow["nLLP_NotInATLAS"] = len(LLPsNotInATLAS.index)
    cutFlow["nLLP_NotInATLAS_weighted"] = LLPsNotInATLAS["weight"].sum()
    cutIndices["nLLP_NotInATLAS"] = LLPsNotInATLAS.index.tolist()

    """
    print("outside detecter cut")
    passed_count = LLPsNotInATLAS.shape[0]
    failed_count = total - passed_count
    print(passed_count)
    print(f"Passed: {passed_count} ({passed_count/total:.2%})  | Failed: {failed_count} ({failed_count/total:.2%})")
    """

    # Check that the LLP direction would pass through at least two layers of the ANUBIS Detector
    temp = LLPsNotInATLAS.apply(h.checkIntersectionsWithANUBIS, args=(selection["geometry"], "decayVertex_weighted"), axis=1)
    if len(temp)!=0:
        LLPsNotInATLAS["intersectionsWithANUBIS"] = LLPsNotInATLAS.apply(h.checkIntersectionsWithANUBIS, args=(selection["geometry"], "decayVertex_weighted"), axis=1)
        LLPsIntersect = LLPsNotInATLAS[LLPsNotInATLAS.apply(lambda row: len(row["intersectionsWithANUBIS"]) >=2, axis=1)]
    else:
        LLPsIntersect=LLPsNotInATLAS
    """
    print("directionally correct cut")
    passed_count = LLPsIntersect.shape[0]
    failed_count = total - passed_count
    print(passed_count)
    print(f"Passed: {passed_count} ({passed_count/total:.2%})  | Failed: {failed_count} ({failed_count/total:.2%})")
    """

    LLPsGeo = LLPsIntersect

    cutFlow["nLLP_Geometry"] = len(LLPsGeo.index)
    cutFlow["nLLP_Geometry_weighted"] = LLPsGeo["weight"].sum()
    cutIndices["nLLP_Geometry"] = LLPsGeo.index.tolist()


    # Require two hits in ANUBIS tracking stations from LLP decay products (where a 'hit' is that a charged track intersects the tracking stations).
    LLPsDecays = h.checkDecayHits(LLPsGeo, SDFs["LLPchildren"], nHits = selection["nStations"], requireCharge=True, geometry=selection["geometry"], 
                                  decayVertex="decayVertex_weighted") 
    cutFlow["nLLP_Tracker"] = len(LLPsDecays.index)
    cutFlow["nLLP_Tracker_weighted"] = LLPsDecays["weight"].sum()
    cutIndices["nLLP_Tracker"] = LLPsDecays.index.tolist()

    """
    print("2 hit requirement cut")
    passed_count = LLPsDecays.shape[0]
    failed_count = total - passed_count
    print(passed_count)
    print(f"Passed: {passed_count} ({passed_count/total:.2%})  | Failed: {failed_count} ({failed_count/total:.2%})")
    """

    # Perform a MET selection - LLPs need to have a minimum associated MET
    LLPsMET = LLPsDecays[LLPsDecays.MET > selection["minMET"]]
    cutFlow["nLLP_MET"] = len(LLPsMET.index)
    cutFlow["nLLP_MET_weighted"] = LLPsMET["weight"].sum()
    cutIndices["nLLP_MET"] = LLPsMET.index.tolist()

    # Next isolation cuts are applied, these are done on a per event basis
    #   - Delta R requirement enforced between LLP and Jet, if there is no Jet then between LLP and charged particle tracks
    minDeltaRColumns = list(LLPsMET.apply(h.getMinDeltaR, args=(SDFs,selection), axis=1)) # Create columns of the Isolation variables
    if len(minDeltaRColumns)==0: # Due to LLPsMET being an empty dataframe
        LLPsMET[["minDeltaR_Jets","minDeltaR_Tracks"]] = None
    else:
        LLPsMET[["minDeltaR_Jets","minDeltaR_Tracks"]] = minDeltaRColumns
    
    # Get the isolated cases for just jets and just charged particles
    #   - min DeltaR of -1 is the case where there are no associated jets or charged tracks in the event
    #LLPsIsoJets = LLPsMET[LLPsMET.eventNumber.isin(passEventIdx["jet"])] 
    LLPsIsoJets = LLPsMET[(LLPsMET.minDeltaR_Jets > selection["minDR"]["jet"]) | (LLPsMET.minDeltaR_Jets==-1)] 
    cutFlow["nLLP_IsoJet"] = len(LLPsIsoJets.index)
    cutFlow["nLLP_IsoJet_weighted"] = LLPsIsoJets["weight"].sum()
    cutIndices["nLLP_isoJet"] = LLPsIsoJets.index.tolist()

    #LLPsIsoCharged = LLPsMET[LLPsMET.eventNumber.isin(passEventIdx["charged"])] 
    LLPsIsoCharged = LLPsMET[(LLPsMET.minDeltaR_Tracks > selection["minDR"]["chargedTrack"]) | (LLPsMET.minDeltaR_Tracks==-1)] 
    cutFlow["nLLP_IsoCharged"] = len(LLPsIsoCharged.index)
    cutFlow["nLLP_IsoCharged_weighted"] = LLPsIsoCharged["weight"].sum()
    cutIndices["nLLP_isoCharged"] = LLPsIsoCharged.index.tolist()

    #LLPsIsoAll = LLPsMET[LLPsMET.eventNumber.isin(passEventIdx["both"])] 
    LLPsIsoAll = LLPsIsoJets[(LLPsIsoJets.minDeltaR_Tracks > selection["minDR"]["chargedTrack"]) | (LLPsMET.minDeltaR_Tracks==-1)] 
    cutFlow["nLLP_IsoAll"] = len(LLPsIsoAll.index)
    cutFlow["nLLP_IsoAll_weighted"] = LLPsIsoAll["weight"].sum()
    cutIndices["nLLP_isoAll"] = LLPsIsoAll.index.tolist()
    
    #loop over surviving events
    #possibly log both llps and one in anubis one in atlas
    
    LLPinclusive = h.includeAtlasPartners(LLPsIsoAll, LLPsInATLAS)
    cutFlow["nLLP_partners"] = len(LLPinclusive.index)
    cutFlow["nLLP_partners_weighted"] = LLPinclusive["weight"].sum()
    cutIndices["nLLP_partners"] = LLPinclusive.index.tolist()
    

    
    outputDict = {"cutFlow": cutFlow,
                  "cutIndices": cutIndices,
                  "finalDF": LLPsIsoAll,
                  "finalfinalDF": LLPinclusive,
                  }

    # Turning back on the SettingWithCopyWarning
    pd.options.mode.chained_assignment = 'warn'

    return outputDict 



# Debugging Text file to check the default lifetime distribution
def createCTauFile(dataFrame, runID, infoDict, suffix="", outputDir=""):
    outDir=f"{outputDir}/ctaus_DS/{runID}/"
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    with open(f"{outDir}/HEPMC_DS_ctau_{infoDict[f'{runID}']['LLPmass']}.txt","w") as f:
        f.write("mass,width,Average DecayDist/mm,Average ctau/mm\n")
        f.write(f"{infoDict[f'{runID}']['LLPmass']},{infoDict[f'{runID}']['LLPdecayWidth']},{dataFrame['decayVertexDist'].mean()},{dataFrame['ctau'].mean()}\n")


def saveCutDict(cutDict, runID, truerunID, infoDict, outputDir="", suffix="", verbose=True, saveAsPickle=False):
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    #fileName=outputDir+f"/CutDict_C{infoDict[f'{runID}']['coupling']}_M{infoDict[f'{runID}']['LLPmass']}_ID{runID}{suffix}".replace('.','p')
    #Using lifetime instead as we are manually  setting the lifetime instead of using the coupling
    LLPlifetime = args.lifetime
    fileName=outputDir+f"/CutDict__lifetime{LLPlifetime}_Mass{infoDict[f'{truerunID}']['LLPmass']}_ID{runID}{suffix}_".replace('.','p')
    with open(fileName+".json", "w") as js:
        json.dump(cutDict["cutFlow"], js, indent = 2)

    if saveAsPickle:
        with open(fileName+".pickle", "wb") as pkl:
            pickle.dump(cutDict, pkl)

    if verbose:
        print(json.dumps(cutDict["cutFlow"], indent = 2))

        print(f"Saved cutdict to {fileName}.json")

def getLifetimeFromCSV(sampleInfo, config):
    if "lifetimeFile" not in config.keys():
        raise Exception(f"No lifetime file provided in config with key 'lifetimeFile' and no manual lifetime set")

    # lifetimeFile is assumed to be a csv with the format: mass, coupling, lifetime
    #   - Different models may require different formats, adjust as needed.
    lifetime_df = pd.read_csv(config["lifetimeFile"])
    ctau = lifetime_df.loc[(lifetime_df["Mass"] == sampleInfo["LLPmass"]) & (lifetime_df["Coupling"] == sampleInfo["coupling"]), "Lifetime"].values
    if ctau.size > 0:
        ctauVal = ctau[0]  # Extract the first value (if multiple rows match)
        lifetime = float(ctauVal / constants.c)
        print(f"The lifetime value is: {lifetime}, ctau = {ctauVal}")
    else:
        print("No matching lifetime found for the given mass and coupling.")
        lifetime = ""

    return lifetime


# Per file processing of dataframes.
#   TODO: Adapt to allow for input of multiple files at a time? 
def selectionProcedure(runID, truerunID, config, sampleInfo, cavern, lifetime=-99, outDir="", dataFrameDir="", cutDictDir="", suffix="", phiFold=True, 
                       remakeDFs=True, cacheDFs=False, lifetimeReweighting=False, cacheSelectedDF=False, makePlots=False, verbose=True):
    print("test2")
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    #sampleInfo = sampleInfoDicts[runID]
    selection=config["selection"]
    selection["geometry"] = cavern # Add the Geometry Class object for geometric acceptance
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
    
    unknownPIDs, failedFiles = [], ""
    print("test3")

    #TODO: In future this path should come from the database itself (i.e. sampleInfo) rather than being constructed as below
    fileName = f"{config['SimSwarmDir']}/{sampleInfo['SimSwarmName']}/Events/{sampleInfo['runName']}_decayed_1/tag_1_pythia8_events.hepmc"
    
    print("test4")

    #TODO: Add a delete to the unzipped filename at the end of a run to save memory. 
    if not os.path.exists(fileName):
        if os.path.exists(fileName+".gz"):
            print("sample hepmc file exists. unzip now.")
            file_to_unzip = fileName+".gz"
            h.unzip_file(file_to_unzip, fileName)
        else:
            print("sample hepmc NOT FOUND. simulation missing.")  
            failedFiles = fileName

            return {} # Return early an empty dict

    start = datetime.datetime.now()
    print(f"Processing: {fileName}")
    if verbose:
        print(f"Start Time: {start}")

    #================================================#
    #   Process full HEPMC file to pandas dataframe  #
    #================================================#
    if dataFrameDir=="":
        dataFrameDir = "./dataframes"

    if not os.path.exists(dataFrameDir):
        os.makedirs(dataFrameDir)

    parquetName = f"{dataFrameDir}/tempDF_{runID}.parquet"
    if not os.path.isfile(parquetName) or remakeDFs:
        df, unknownPIDs = createDFfromHEPMC(fileName)

        #batch running so want multiple diferent files
        cacheDFs = False
        if cacheDFs: # Cache pandas dataframes as parquet files to avoid having to reprocess the file
            df.to_parquet(parquetName)
    else:
        print(f"Loading dataframe from {parquetName}...")
        df = pd.read_parquet(parquetName)
 
    if os.path.exists(fileName):
        os.remove(fileName) # Remove the unzipped file to save space
    #=====================================#
    #  Perform Phi Folding if specified   #
    #=====================================#
    print(f"phiFold is set to: ", phiFold)
    print(df[(df.PID == 35)][["eventNumber", "PID", "phi"]])

    phiFold = False

    if phiFold:
        if verbose:
            print(f"Performing Phi Folding...")
            phiFoldStartTime = datetime.datetime.now()
            print(f"Phi Fold start time: {phiFoldStartTime}")
        print("Max phi before folding:", df[df.PID.isin([sampleInfo["LLPid"]])]["phi"].max())
        print("Min folding:", df[df.PID.isin([sampleInfo["LLPid"]])]["phi"].min())

        # To Maximise statistics on the LLPs we can study reflect LLPs that travel downwards, upwards instead.
        # In this case to ensure that selections are performed accurately, the other particles in the event must also be reflected.
        # Rotation by pi in phi is equivalent to a reflection in x and y.
        # Due to the phi invariance, this will artificial enhance our sample as required. So the final raw numbers for selected events must be halved.
        # By the phi definition 0 is the x-axis and increases up to pi anticlockwise, or decreases down to -pi clockwise
        foldEventsList = df[(df.PID.isin([sampleInfo["LLPid"]])) & (df["phi"] < 0)]["eventNumber"].tolist()
        print(foldEventsList)
        temp = df[df.eventNumber.isin(foldEventsList)].apply(doPhiFold, args=(foldEventsList,),axis=1)
        
        df.loc[temp.index,["px","py","phi","prodVertex","decayVertex"]] = pd.DataFrame(temp.tolist(), 
                                                                                       columns=["px","py","phi","prodVertex","decayVertex"],
                                                                                       index=temp.index)
        print("Max phi after folding:", df[df.PID.isin([sampleInfo["LLPid"]])]["phi"].max())
        print("Min folding:", df[df.PID.isin([sampleInfo["LLPid"]])]["phi"].min())
        print("phi > pi after folding:")
        print(df[(df.PID.isin([sampleInfo["LLPid"]])) & (df["phi"] > np.pi)][["eventNumber", "PID", "phi"]])


        if verbose:
            phiFoldEndTime = datetime.datetime.now()
            print(f"Done Phi Folding! ({phiFoldEndTime})")
            print(f"Took {phiFoldEndTime - phiFoldStartTime}")

    #======================================================================#
    # Divide the total Dataframe into several smaller dataframes           #
    #   - e.g. LLPs, LLP children, chargedFinalStates, neutralFinalStates  #
    #======================================================================#
    if verbose:
        sampleDFstartTime = datetime.datetime.now()
        print(f"Creating sample dataframes... ({sampleDFstartTime})")
    
    sampleDFs = createSampleDataFrames(df, sampleInfo["LLPid"], selection["minPt"])
    myplot.plotStandardPlots(sampleDFs, outDir=outputDir, suffix=f"{suffix}_Run{runID}")

    if verbose:
        sampleDFendTime = datetime.datetime.now()
        print(f"Done creating sample dataframes... ({sampleDFendTime})")
        print(f"Took: {sampleDFendTime - sampleDFstartTime}")

    #===================================================================================#
    # Reweight the decay positions of the LLP(s) based on a given lifetime if specified #
    #===================================================================================#
    if lifetimeReweighting:
        if verbose:
            reweightingStartTime = datetime.datetime.now()
            print(f"Reweighting LLPs by lifetime... ({reweightingStartTime})")
        #LLPlifetime = 0.00001

        if lifetime == -99:
            LLPlifetime = getLifetimeFromCSV(sampleInfo, config)

        else:
            LLPlifetime = args.lifetime

    

        if LLPlifetime=="":
            print(f"Unable to obtain a lifetime, skipping file")
            failedFiles = fileName
            return {}
        
        print(LLPlifetime)

        reweightDecayPositions(sampleDFs["LLPs"], sampleDFs["LLPchildren"], LLPlifetime, sampleInfo["LLPid"], seed = config["seed"]*int(runID))
        
        if verbose:
            reweightingEndTime = datetime.datetime.now()
            print(f"Finished reweighting LLPs by lifetime! ({reweightingEndTime})")
            print(f"Took: {reweightingEndTime - reweightingStartTime}")

    #==============================================#
    # Create Jet objects using the fastJet package #
    #==============================================#
    sampleDFs["finalStatePromptJets"] = createJetDF(list(set(df['eventNumber'])), sampleDFs["chargedFinalStates"], sampleDFs["neutralFinalStates"])

    #=================================#
    #   Apply Selection to the LLPs   #
    #=================================#
    if verbose:
        selectionStartTime = datetime.datetime.now()
        print(f"Applying Selection... ({selectionStartTime})")

    if cutDictDir!="":
        if not os.path.exists(cutDictDir):
            os.makedirs(cutDictDir)
    else:
        cutDictDir = f"./CutDicts/Run{runID}/"

    cutDict =  applySelection(sampleDFs, selection=selection)
    saveCutDict(cutDict, runID, truerunID, sampleInfoDicts, outputDir=cutDictDir, suffix=suffix, verbose=verbose, saveAsPickle=cacheSelectedDF)
    sampleDFs["LLPs_selected"] = cutDict["finalDF"]

    if verbose:
        selectionEndTime = datetime.datetime.now()
        print(f"Finished applying Selection... ({selectionEndTime})")
        print(f"Took: {selectionEndTime - selectionStartTime}")

    #============================================================#
    #   Debug test to see the unweighted lifetime distributions  #
    #============================================================#
    #NOTE: This probably could be removed? It was used to compare the simulation sample lifetimes to those from Sofie and HNLCalc's predictions
    createCTauFile(sampleDFs["LLPs"], runID, sampleInfoDicts, suffix=f"{suffix}_Run{runID}", outputDir="./")
    if verbose:
        print(f"LLP Mass: {np.mean(sampleDFs['LLPs']['mass'])}")
        print(f"Before Selection:")
        print(f"\tAverage LLP Decay Distance (mm): {sampleDFs['LLPs']['decayVertexDist'].mean()}")
        print(f"\tAverage LLP Weighted Decay Distance (mm): {sampleDFs['LLPs']['decayVertexDist_weighted'].mean()}")
        print(f"\tAverage LLP ctau (mm): {sampleDFs['LLPs']['ctau'].mean()}")
        print(f"\tAverage LLP Weighted ctau (mm): {sampleDFs['LLPs']['ctau_weighted'].mean()}")
        print(f"After Selection:")
        print(f"\tAverage LLP Decay Distance (mm): {sampleDFs['LLPs_selected']['decayVertexDist'].mean()}")
        print(f"\tAverage LLP Weighted Decay Distance (mm): {sampleDFs['LLPs_selected']['decayVertexDist_weighted'].mean()}")
        print(f"\tAverage LLP ctau (mm): {sampleDFs['LLPs_selected']['ctau'].mean()}")


    makePlots = True
    if makePlots:
        if verbose:
            plotStartTime = datetime.datetime.now()
            print(f"Plotting... ({plotStartTime})")
        print(f"Plotting to {outputDir}")

        # Plot the Number of events that pass each stage of the cutFlow
        myplot.plotCutFlow(cutDict["cutFlow"], logScale = True, outputDir=outputDir, suffix=f"{suffix}_Run{runID}")
        
        # Plot each of the columns of the dataframe for each item in the sampleDFs
        myplot.plotStandardPlots(sampleDFs, outDir=outputDir, suffix=f"{suffix}_Run{runID}")

        # MET derived by different methods -- debugging tool. 
        #   - NOTE: This could be commented out
        myplot.plotMETcomparison(sampleDFs, outputDir=outputDir, suffix=f"{suffix}_Run{runID}")

        # Plot the LLP decay positions compared to the ATLAS/ANUBIS geometry at different stages of selection
        #   - Split off the Geometry positions into a subfolder
        geoOutputDir = f"{outputDir}/geoPlots/"
        if not os.path.exists(geoOutputDir):
            os.makedirs(geoOutputDir)

        myplot.plotPositionsDuringSelection(sampleDFs["LLPs"], cutDict, selection, nBins=[100,100,100], decayVertex="decayVertex_weighted", 
                                            outDir=geoOutputDir, suffix=suffix)
        if verbose:
            plotEndTime = datetime.datetime.now()
            print(f"Finished plotting... ({plotEndTime})")
            print(f"Took: {plotEndTime - plotStartTime}")

    print(f"Finished processing file")
    if verbose:
        endTime = datetime.datetime.now()
        print(f"End Time: {endTime}")
        print(f"Overall Took: {endTime - start}")

    return {"sampleDFs": sampleDFs, "cutDict": cutDict, "unknownPIDs": unknownPIDs, "failedFiles": failedFiles}



if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-c","--coupling", nargs="+", type=float, default = [1.0])
    parser.add_argument("-m","--mass", nargs="+", type=float, default = [1.0])
    parser.add_argument("-s","--sample", nargs="+", type=str, default = ["CCDY_eev"])
    parser.add_argument("--runID", nargs="+", type=float, default = [])
    #cuz looping with lifetime so will have higher number but then have to index csv so cant go out of bound
    parser.add_argument("--truerunID", nargs="+", type=float, default = [])
    parser.add_argument("--minRunID", default = "")
    parser.add_argument("--maxRunID", default = "")
    parser.add_argument("--lifetime", type=float, default = 0.00000003)
    parser.add_argument("--outDir", type=str, default = "selectedPlots") 
    parser.add_argument("--suffix", type=str, default = "") 
    parser.add_argument("--phiFold", action="store_true", default = True) 
    parser.add_argument("--remake", action="store_true", default = True) 
    parser.add_argument("--noCacheDFs", action="store_true", default = False) 
    parser.add_argument("--noLifeReweight", action="store_true", default = False) 
    parser.add_argument("--config", type=str, default = "configs/config_HNL_m1_c1.yaml", 
                        help = "Configuration file containing necessary things for the run", required=True )
    args = parser.parse_args()

    # Create outDir if it doesn't exist
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)
    outputDir = args.outDir
    

    args.phiFold = True
    if args.phiFold:
        print("phifold suffix check")
        args.suffix+="_PhiFolded"

    # Process the Configuration file
    config = h.parseConfig(args.config)

    # Process the Simulation Information table
    sampleInfoDicts = h.parseCSV(config["scanFile"])
    runIDsByMassCoupling = h.groupMassAndCoupling(sampleInfoDicts) # var[sample][mass][coupling] = list of associated runIDs

    # Create/Load a version of the ATLAS Cavern and ANUBIS Geometry 
    geometryPath = config["geometryPath"] + "geoATLAS_ANUBIS.pickle"
    cavern = h.createGeometry(geometryPath)

    # Hard-coding some inputs for testing
    #runID = '599' # Used for electrons
    print(sampleInfoDicts)
    print("Available runIDs:", sampleInfoDicts.keys())



    unknownPIDs=[]
    failedFiles={}

    runID = args.runID
    truerunID = args.truerunID
    #disabling for loops for now as trying to loop batches instead of in the files
    #for runID in sampleInfoDicts.keys():
    print(runID)
    
    #for runID in sampleInfoDicts.keys():
    print("test1")

    #for non batch run testing
    runID = [4]
    truerunID = [4]


    runID = runID[0]
    runID = str(int(runID))

    truerunID = str(int(truerunID[0]))
    sampleInfo = sampleInfoDicts[truerunID]
    """
    if int(runID) not in args.runID and len(args.runID)!=0:
        continue

    if ( (args.minRunID !="" and int(runID) < int(args.minRunID)) or
        (args.maxRunID !="" and int(runID) > int(args.maxRunID)) ):
        continue
    """
    print(f"== Run ID: {runID}")
    print(sampleInfo)

    outputDir=args.outDir+f"/Run{runID}/"
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    # selection procedure does the procedure on a per file basis for better versatility
    #   - outputs dict with: 'sampleDFs' pandas dataframes for different particle types e.g. LLPs, LLPchildren, finalStateChargedTracks etc, 
    #                        'cutDict' a dictionary of the LLP events per cut, the associated indices from sampleDFs and the final selected dataframe,
    #                        'failedFiles' a list of files that couldn't be opened,
    #                        'unknownPIDs' a list of PIDs not included in the charge module to determine the charge of the particle.
    selectionDict =  selectionProcedure(runID, truerunID, config, sampleInfo, cavern, lifetime=args.lifetime, outDir=outputDir, suffix= args.suffix, 
                                        phiFold=args.phiFold, remakeDFs=args.remake, cacheDFs=not args.noCacheDFs, 
                                        lifetimeReweighting=not args.noLifeReweight, cacheSelectedDF=not args.noCacheDFs, 
                                        makePlots=False, verbose=True)
    
    unknownPIDs.extend(selectionDict["unknownPIDs"])
    failedFiles[runID] = selectionDict["failedFiles"]

    print(f"Unknown PIDs: {list(set(unknownPIDs))}")
    print("======================================================\n\n")
    print("======================================================")


print(json.dumps(failedFiles, indent=2))


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
    # plt.show()

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
    

"""
# Example of how to call this function:
# Assuming 'LLPs', 'LLPchildren', and 'finalStates' are your dataframes
select_event(sampleDFs['LLPs'], sampleDFs["LLPchildren"], sampleDFs["finalStates"])
"""