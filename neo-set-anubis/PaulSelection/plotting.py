import sys, os
import matplotlib.pyplot as plt
import numpy as np 
import pandas as pd

# For plotting assume that we have pandas dataframes from the selectLLPs.py script, which have the following columns available in them typically:
# ["eventNumber","particleIndex","px","py","pz","pt","E","mass","prodVertex","prodVertexDist","decayVertex","decayVertexDist",
#  "boost","phi","eta","MET","theta","beta","PID","charge","nParents","nChildren","parentIndices","childrenIndices",
#  "weight","status","ctau"]

def plotColumnAsHist(df, columnName, nBins=None, axisRange=None, normalise=False, logScale=False, outputDir="", suffix="", metcutoff = False, expfit = False):
    fig, ax = plt.subplots(figsize=(14, 10))
    
    xData = df[columnName].to_numpy()
    binSize = abs(max(xData) - min(xData))/nBins if nBins!=None else "Bin"

    ax.set_xlabel(columnName, fontsize = 20)
    ax.set_ylabel(f'Counts / {binSize}', fontsize=20)
    ax.hist(xData, bins=nBins, range=axisRange, density = normalise, 
            color = "coral", ec = "red", label=columnName)
    if logScale:
        ax.set_yscale('log')
    
    if metcutoff:
        plt.axvline(x=30, color='r', linestyle='--', label="MET Cut-off")

    if expfit:

        import scipy as sp
        def exp_decay(x, a, b):
            return a * np.exp(-b * x)
        counts, bin_edges = np.histogram(xData, bins=16, range=(-10,150), density=False)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Compute bin centers
        valid = bin_centers > 0
        popt, pcov = sp.optimize.curve_fit(exp_decay, bin_centers, counts, p0=(max(counts), 0.1))

        x_fit = np.linspace(min(bin_centers), max(bin_centers), 100)
        y_fit = exp_decay(x_fit, *popt)
        plt.plot(x_fit, y_fit, 'r-', linewidth=2, label="Exponential Fit")
    
    plt.grid()              

    plotName=outputDir+f"/LLPDF_{columnName}{suffix}.pdf"
    plt.savefig(plotName)
    plt.close(fig)

    return plotName


def plotColumnsAs2DHist(df, columnNames, nBins=None, axisRange=None, normalise=False, logScale=[False, False], outputDir="", suffix=""):
    fig, ax = plt.subplots(1, figsize=(14, 10), dpi=100)

    if len(columnNames)!=2:
        print("This function requires ONLY two columns provided.")
        return -1

    if len(nBins) not in [0,2]:
        print("This function requires either two values in nBins ([nBinX, nBinY]) or to be empty.")
        return -1
    
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
    fig.savefig(f"{outputDir}/AllEventTimes_Summary{suffix}.pdf", bbox_inches="tight")
    plt.close(fig)   


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

    fig, ax = plt.subplots(1, figsize=(10, 10), dpi=100)
    p = plt.bar(x, y, width=0.6, label = "CutFlow", fill=False)
    ax.set_ylabel('N events')
    if logScale:
        ax.set_yscale('log')
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
    plt.xticks(rotation=50)
    plt.tight_layout()

    plotName = outputDir + f"/CutFlow_weighted{suffix}.pdf" 
    plt.savefig(plotName)
    plt.close(fig)
