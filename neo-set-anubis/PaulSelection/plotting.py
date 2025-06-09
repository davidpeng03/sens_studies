import sys, os
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd
import helper as h
pd.set_option('use_inf_as_na', True)
from matplotlib.backends.backend_pdf import PdfPages

# For plotting assume that we have pandas dataframes from the selectLLPs.py script, which have the following columns available in them typically:
# ["eventNumber","particleIndex","px","py","pz","pt","E","mass","prodVertex","prodVertexDist","decayVertex","decayVertexDist",
#  "boost","phi","eta","MET","theta","beta","PID","charge","nParents","nChildren","parentIndices","childrenIndices",
#  "weight","status","ctau"]

def plotColumnAsHist(df, columnName, nBins=40, axisRange=None, normalise=False, logScale=False, pp="", outputDir="", suffix=""):
    fig, ax = plt.subplots(figsize=(14, 10))
    
    xData = (df[columnName].dropna()).to_numpy()
    binSize = abs(max(xData) - min(xData))/nBins if nBins!=None and len(xData)!=0 else "Bin"

    ax.set_xlabel(columnName, fontsize = 20)
    ax.set_ylabel(f'Counts / {binSize}', fontsize=20)
    ax.hist(xData, bins=nBins, range=axisRange, density = normalise, 
            color = "coral", ec = "red", label=columnName)
    if logScale:
        ax.set_yscale('log')
    plt.grid()              

    plotName=outputDir+f"/LLPDF_{columnName}{suffix}.pdf"
    if pp == "":
        plt.savefig(plotName)
    else:
        pp.savefig(fig, bbox_inches="tight")
    plt.close(fig)

def plotPositionsAsHist(df, columnName, nBins=40, axisRange=None, normalise=False, logScale=False, pp="", outputDir="", suffix=""):
    posData = (df[columnName].dropna()).to_numpy()
    pos={}
    pos["x"] = [x[0] for x in posData]
    pos["y"] = [x[1] for x in posData]
    pos["z"] = [x[2] for x in posData]
    #pos["t"] = [x[3] for x in posData]

    plotNames=[]
    for label, position in pos.items():
        binSize = abs(max(position) - min(position))/nBins if nBins!=None and len(position)!=0 else "Bin"

        fig, ax = plt.subplots(figsize=(14, 10))
        ax.set_xlabel(f"{columnName}, {label}", fontsize = 20)
        ax.set_ylabel(f'Counts / {binSize}', fontsize=20)
        ax.hist(position, bins=nBins, range=axisRange, density = normalise, 
                color = "coral", ec = "red", label=columnName)
        plt.grid()              
        if logScale:
            ax.set_yscale('log')

        plotName=outputDir+f"/LLPDF_{columnName}_{label}Coord{suffix}.pdf"
        if pp=="":
            plt.savefig(plotName)
        else: 
            pp.savefig(fig, bbox_inches="tight") # Save to a PdfPage
        plt.close(fig)

        plotNames.append(plotName)


def plotColumnsAs2DHist(df, columnNames, nBins=40, axisRange=None, normalise=False, logScale=[False, False], pp="", outputDir="", suffix=""):
    fig, ax = plt.subplots(1, figsize=(14, 10), dpi=100)

    if len(columnNames)!=2:
        print("This function requires ONLY two columns provided.")
        return -1

    if len(nBins) not in [0,2]:
        print("This function requires either two values in nBins ([nBinX, nBinY]) or to be empty.")
        return -1
    
    df = df.dropna(subset=columnNames)
    xData = df[columnNames[0]].to_numpy()
    yData = df[columnNames[1]].to_numpy()
    binSizeX = "Bin"
    binSizeY = "Bin"
    if len(nBins)!=0:
        binSizeX = abs(max(xData) - min(xData)) / nBins[0] if nBins[0] !=None else "Bin"
        binSizeY = abs(max(yData) - min(yData)) / nBins[1] if nBins[1] !=None else "Bin"

    counts, xedges, yedges, im = ax.hist2d(xData, yData, nBins, range=axisRange, cmin=1)
    fig.colorbar(im, ax=ax) 
    plt.xlabel(f"{columnNames[0]} / {binSizeX}")     
    plt.ylabel(f"{columnNames[1]} / {binSizeY}")     
    
    if logScale[0]:
        ax.set_xscale('log')
    if logScale[1]:
        ax.set_yscale('log')

    plt.grid()              
    fig.tight_layout()        
    if pp =="":
        fig.savefig(f"{outputDir}/AllEventTimes_Summary{suffix}.pdf", bbox_inches="tight")
    else: 
        pp.savefig(fig, bbox_inches="tight") # Save to a PdfPage
    plt.close(fig)   

def plotStandardPlots(sampleDFs, outDir="", suffix=""):

    for dfType, df in sampleDFs.items():
        columns = df.columns.values.tolist()
        print(suffix)
        with PdfPages(f"{outDir}/DF{dfType}_Plots{suffix}.pdf") as pdf:
            for col in columns:
                print(col)
                if "Indices" in col or col in ["intersectionsWithANUBIS"]:
                    continue

                if "Vertex" in col and "Dist" not in col:
                    plotPositionsAsHist(df, col, pp=pdf)
                else:
                    if "Vertex" in col: 
                        doLog=True
                    else:
                        doLog=False
                    plotColumnAsHist(df, col, logScale=doLog, pp=pdf)
                    

    plt.close('all') # Ensure all plots are closed to help with memory usage


def plotDecayVertexPositions(df, geometry, ranges={"x": [], "y": [], "z": []}, nBins=[100,100,100], decayVertex="decayVertex", outDir="", suffix=""):
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    hitsDecay = {"passed": [], "failed": []}
    for vertex in (df[decayVertex]).tolist():
        # Convert coordinates to be in the centre of Cavern coord system in /m 
        X, Y, Z = geometry.coordsToOrigin(vertex[0]*1E-3, vertex[1]*1E-3, vertex[2]*1E-3)
        hitsDecay["passed"].append([X,Y,Z])

    defaultCavernBounds = {"x": [-18,18], "y": [-15,25], "z": [-30,30]}

    bins={"nX": nBins[0], "nY": nBins[1], "nZ": nBins[2]}
    for coord in ["x","y","z"]:
        if ranges[coord]== None:
            bins[f"range{coord.upper()}"] = None
        elif len(ranges[coord])==0:
            bins[f"range{coord.upper()}"] = defaultCavernBounds[coord]
        elif len(ranges["x"])==2:
            bins[f"range{coord.upper()}"] = ranges[coord]
        else:
            bins[f"range{coord.upper()}"] = None
    hitsDecay["bins"] = bins

    tempRanges={"xy": {"x": bins["rangeX"], "y": bins["rangeY"]},
            "xz": {"x": bins["rangeX"], "z": bins["rangeZ"]},
            "zy": {"z": bins["rangeZ"], "y": bins["rangeY"]},
            "3D": {"x": bins["rangeX"], "y": bins["rangeY"], "z": bins["rangeZ"]}
    }
    
    geometry.plotFullCavern(hitsDecay, simpleAnubisRPCs=geometry.ANUBIS_RPCs, plotATLAS=True, plotAcceptance=True, ranges=tempRanges, suffix=suffix, outDir=outDir)
    plt.close('all') # Ensure all plots are closed to help with memory usage

def plotPositionsDuringSelection(LLPs, cutDict, selection, nBins=[100,100,100], decayVertex="decayVertex", outDir="", suffix=""):
    # Plot the Positions before the geometric cuts 
    plotDecayVertexPositions(LLPs[LLPs.index.isin(cutDict["cutIndices"]["nLLP_LLPdecay"])], selection["geometry"], ranges={"x": [-40,40], "y": [-40,40], "z": [-40,40]}, 
                             nBins=nBins, outDir=outDir, suffix=f"BeforeGeoCut{suffix}", decayVertex=decayVertex)

    # Plot the Positions after requiring the hits to be within the Cavern
    plotDecayVertexPositions(LLPs[LLPs.index.isin(cutDict["cutIndices"]["nLLP_InCavern"])], selection["geometry"], nBins=nBins, 
                             outDir=outDir, suffix=f"AfterInCavernCut{suffix}", decayVertex=decayVertex)

    # Plot the Positions after requiring the hits to be within the Cavern and not in ATLAS
    plotDecayVertexPositions(LLPs[LLPs.index.isin(cutDict["cutIndices"]["nLLP_NotInATLAS"])], selection["geometry"], nBins=nBins, 
                             outDir=outDir, suffix=f"AfterNotInATLASCut{suffix}", decayVertex=decayVertex)

    # Plot the Positions after requiring the hits to be within the Cavern, not in ATLAS, and with LLP tracks pointing towards ANUBIS
    plotDecayVertexPositions(LLPs[LLPs.index.isin(cutDict["cutIndices"]["nLLP_Geometry"])], selection["geometry"], nBins=nBins, 
                             outDir=outDir, suffix=f"AfterGeoCut{suffix}", decayVertex=decayVertex)

    # Plot the Positions after requiring the hits to be within the Cavern, not in ATLAS, with LLP tracks pointing towards ANUBIS, and tracks from the LLP decay products
    plotDecayVertexPositions(LLPs[LLPs.index.isin(cutDict["cutIndices"]["nLLP_Tracker"])], selection["geometry"], nBins=nBins, 
                             outDir=outDir, suffix=f"AfterChildCut{suffix}", decayVertex=decayVertex)

    # Plot the Positions after requiring the hits to be within the Cavern, not in ATLAS, with LLP tracks pointing towards ANUBIS, with tracks from the LLP decay products and MET requirement
    plotDecayVertexPositions(LLPs[LLPs.index.isin(cutDict["cutIndices"]["nLLP_MET"])], selection["geometry"], nBins=nBins, 
                             outDir=outDir, suffix=f"AfterMETCut{suffix}", decayVertex=decayVertex)

    # Plot the Positions after requiring the hits to be within the Cavern, not in ATLAS, with LLP tracks pointing towards ANUBIS, with, tracks from the LLP decay products, with MET requirement and Iso requirement 
    if "nLLP_IsoJet" in cutDict["cutIndices"]:
        plotDecayVertexPositions(LLPs[LLPs.index.isin(cutDict["cutIndices"]["nLLP_IsoJet"])], selection["geometry"], nBins=nBins, 
                             outDir=outDir, suffix=f"AfterIsoJetCut{suffix}", decayVertex=decayVertex)
    
    

    if "nLLP_IsoCharged" in cutDict["cutIndices"]:
        plotDecayVertexPositions(LLPs[LLPs.index.isin(cutDict["cutIndices"]["nLLP_IsoCharged"])], selection["geometry"], nBins=nBins, 
                             outDir=outDir, suffix=f"AfterIsoChargedCut{suffix}", decayVertex=decayVertex)
    
    

    if "nLLP_IsoAll" in cutDict["cutIndices"]:
        plotDecayVertexPositions(LLPs[LLPs.index.isin(cutDict["cutIndices"]["nLLP_IsoAll"])], selection["geometry"], nBins=nBins, 
                             outDir=outDir, suffix=f"AfterIsoCut{suffix}", decayVertex=decayVertex)
    
    else:
        print("Warning: 'nLLP_IsoAll' not found in cutDict['cutIndices']")

    # Plot the other properties in the DataFrame before and after the geometric acceptance 
    plotStandardPlots({"LLP_BeforeGeoCut": LLPs[LLPs.index.isin(cutDict["cutIndices"]["nLLP_LLPdecay"])], 
                       "LLP_AfterGeoCut": LLPs[LLPs.index.isin(cutDict["cutIndices"]["nLLP_Tracker"])]}, outDir=outDir)

    plt.close('all') # Ensure all plots are closed to help with memory usage

def plotCutFlow(cutDict, logScale = False, outputDir="", suffix=""):

    x,y = [],[]
    xw,yw = [],[]
    for cut, val in cutDict.items():
        if "weight" in cut:
            xw.append(cut)
            yw.append(val)
        else:
            x.append(cut)
            y.append(val)
    
    fontsize = 25

    fig, ax = plt.subplots(1, figsize=(10, 10), dpi=100)
    p = plt.bar(x, y, width=0.6, label = "CutFlow", fill=False)
    ax.set_ylabel('N events', fontsize=fontsize)
    ax.tick_params(axis='both', which='major', labelsize=fontsize)  # tick label size

    if logScale:
        ax.set_yscale('log')
    else:
        ax.set_ylim(bottom=0)
    plt.xticks(rotation=50)
    plt.tight_layout()

    plotName = outputDir + f"/CutFlow{suffix}.pdf" 
    plt.savefig(plotName)
    plt.close(fig)

    fig, ax = plt.subplots(1, figsize=(10, 10), dpi=100)
    p = plt.bar(xw, yw, width=0.6, label = "CutFlow_weighted", fill=False)
    ax.set_ylabel('Weighted Number of Events')
    if logScale:
        ax.set_yscale('log')
    else:
        ax.set_ylim(bottom=0)
    plt.xticks(rotation=50, fontsize=fontsize)  # x-ticks font size & rotation
    plt.tight_layout()
    
    plotName = outputDir + f"/CutFlow_weighted{suffix}.pdf" 
    plt.savefig(plotName)
    plt.close(fig)

def plotMETcomparison(sampleDFs, outputDir="", suffix=""):
    METxplot = plotColumnAsHist(sampleDFs["LLPs"], "METx", nBins = 20, normalise=True, logScale=False, outputDir=outputDir, suffix=suffix)
    METyplot = plotColumnAsHist(sampleDFs["LLPs"], "METy", nBins = 20, normalise=True, logScale=False, outputDir=outputDir, suffix=suffix)
    METplot = plotColumnAsHist(sampleDFs["LLPs"], "MET", nBins = 20, normalise=True, logScale=False, outputDir=outputDir, suffix=suffix)

    # Debug MET by getting the pt of the neutrinos and LLPs directly
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

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.set_xlabel("Direct MET", fontsize = 20)
    ax.set_ylabel(f'Counts / Bin', fontsize=20)
    ax.hist(testMET, density = False, bins = 16, range=(-10,150), color = "coral", ec = "red", fill=False, label="DirectMET")
    ax.hist(sampleDFs["LLPs"]["MET"].to_numpy(), range=(-10,150), density = False, bins = 16, color = "turquoise", ec = "teal", fill=False, label="MET")
    plt.grid()              
    plt.legend(loc="best")
    plotName=outputDir+f"/LLPDF_MET_FromLLPsAndNeutrinos{suffix}.pdf"
    plt.savefig(plotName)
    plt.close(fig)
