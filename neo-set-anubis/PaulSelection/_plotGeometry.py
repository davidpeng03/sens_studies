import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import numpy as np

# Define a set of helper plotting functions for the ATLASCavern class

def plotCavernXY(self, ax, plotATLAS=False, plotAcceptance=False): 
    # Get the Cavern ceiling data grid
    cavernArch = self.createCavernVault(doPlot=False)

    # Cavern Boundaries
    ax.plot([self.CavernX[0], self.CavernX[0]], [self.CavernY[0], self.CavernY[1]], c="paleturquoise")
    ax.plot([self.CavernX[1], self.CavernX[1]], [self.CavernY[0], self.CavernY[1]], c="paleturquoise")
    ax.plot([self.CavernX[0], self.CavernX[1]], [self.CavernY[0], self.CavernY[0]], c="paleturquoise")
    # Cavern Ceiling
    ax.plot(cavernArch["x"], cavernArch["y"], c="paleturquoise")
    # Access Shafts
    ax.plot([self.PX14_Centre["x"]-self.PX14_Radius,self.PX14_Centre["x"]-self.PX14_Radius],
                              [self.PX14_LowestY, self.PX14_Centre["y"]+self.PX14_Height], c="red", label="PX14", alpha=0.5)
    ax.plot([self.PX14_Centre["x"]+self.PX14_Radius,self.PX14_Centre["x"]+self.PX14_Radius],
                              [self.PX14_LowestY, self.PX14_Centre["y"]+self.PX14_Height], c="red", label="PX14", alpha=0.5)

    ax.plot([self.PX16_Centre["x"]-self.PX16_Radius,self.PX16_Centre["x"]-self.PX16_Radius],
                              [self.PX16_LowestY, self.PX16_Centre["y"]+self.PX16_Height], c="turquoise", label="PX16", alpha=0.5)
    ax.plot([self.PX16_Centre["x"]+self.PX16_Radius,self.PX16_Centre["x"]+self.PX16_Radius],
                              [self.PX16_LowestY, self.PX16_Centre["y"]+self.PX16_Height], c="turquoise", label="PX16", alpha=0.5)

    # Mark the Cavern Centre, IP, and Centre of Curvature for the ceiling
    ax.scatter(0, 0, c="r", marker = "x", label="Cavern Centre")
    ax.annotate("Centre", (0,0))
    ax.scatter(self.IP["x"], self.IP["y"], c="g", marker = "o", label="IP")
    ax.annotate("IP", (self.IP["x"], self.IP["y"]))
    ax.scatter(self.centreOfCurvature["x"], self.centreOfCurvature["y"], c="b", marker = "D", label="Ceiling Centre of Curvature")
    ax.annotate("Centre of Curvature (Ceiling)", (self.centreOfCurvature["x"], self.centreOfCurvature["y"]))

    if plotATLAS:
       ax.add_patch( plt.Circle((self.ATLAS_Centre["x"], self.ATLAS_Centre["y"]), self.radiusATLAS, color="blue", fill=False, ls="--") )

    if plotAcceptance:
        # Plot a rough impression of the Acceptance
        ax.plot([self.IP["x"], self.CavernX[0]], [self.IP["y"], self.CavernY[1]], c="k", alpha=0.25, linestyle="--")
        ax.plot([self.IP["x"], self.CavernX[1]], [self.IP["y"], self.CavernY[1]], c="k", alpha=0.25, linestyle="--")

    ax.set_xlim(-18,18)
    ax.set_ylim(-15,25)

def plotCavernXZ(self, ax, plotATLAS=False): 
    # Cavern Boundaries
    ax.plot( [self.CavernX[0], self.CavernX[1]], [self.CavernZ[0], self.CavernZ[0]], c="paleturquoise")
    ax.plot( [self.CavernX[0], self.CavernX[1]], [self.CavernZ[1], self.CavernZ[1]], c="paleturquoise")
    ax.plot( [self.CavernX[0], self.CavernX[0]], [self.CavernZ[0], self.CavernZ[1]], c="paleturquoise")
    ax.plot( [self.CavernX[1], self.CavernX[1]], [self.CavernZ[0], self.CavernZ[1]], c="paleturquoise")
    # Access Shafts
    ax.add_patch( plt.Circle((self.PX14_Centre["x"], self.PX14_Centre["z"]), self.PX14_Radius, color="red", fill=False, ls="--") )
    ax.add_patch( plt.Circle((self.PX16_Centre["x"], self.PX16_Centre["z"]), self.PX16_Radius, color="turquoise", fill=False, ls="--") )
    
    # Mark the Cavern Centre, IP, and Centre of Curvature for the ceiling
    ax.scatter(0, 0, c="r", marker = "x", label="Cavern Centre")
    ax.annotate("Centre", (0,0))
    ax.scatter(self.IP["x"], self.IP["z"], c="g", marker = "o", label="IP")
    ax.annotate("IP", (self.IP["x"], self.IP["z"]))
    ax.plot([self.centreOfCurvature["x"]]*2, [self.CavernZ[0], self.CavernZ[1]], c="b", label="Ceiling Centre of Curvature")
    ax.annotate("Centre of Curvature (Ceiling)", (self.centreOfCurvature["x"],self.CavernZ[0]/4))

    if plotATLAS:
        ax.plot( [self.ATLAS_Centre["x"]-self.radiusATLAS,self.ATLAS_Centre["x"]+self.radiusATLAS], [self.ATLAS_Z[0],self.ATLAS_Z[0]], c="b")
        ax.plot( [self.ATLAS_Centre["x"]-self.radiusATLAS,self.ATLAS_Centre["x"]-self.radiusATLAS], [self.ATLAS_Z[0],self.ATLAS_Z[1]], c="b")
        ax.plot( [self.ATLAS_Centre["x"]-self.radiusATLAS,self.ATLAS_Centre["x"]+self.radiusATLAS], [self.ATLAS_Z[1],self.ATLAS_Z[1]], c="b")
        ax.plot( [self.ATLAS_Centre["x"]+self.radiusATLAS,self.ATLAS_Centre["x"]+self.radiusATLAS], [self.ATLAS_Z[0],self.ATLAS_Z[1]], c="b")

    ax.set_xlim(-18,18)
    ax.set_ylim(-30,30)

def plotCavernZY(self, ax, plotATLAS=False, plotAcceptance=False): 
    # Cavern Boundaries
    ax.plot([self.CavernZ[0], self.CavernZ[1]], [self.CavernY[1], self.CavernY[1]], c="paleturquoise", ls="--")
    ax.plot([self.CavernZ[0], self.CavernZ[0]], [self.CavernY[0], self.archRadius+self.centreOfCurvature["y"]], c="paleturquoise")
    ax.plot([self.CavernZ[1], self.CavernZ[1]], [self.CavernY[0], self.archRadius+self.centreOfCurvature["y"]], c="paleturquoise")
    ax.plot([self.CavernZ[0], self.CavernZ[1]], [self.CavernY[0], self.CavernY[0]], c="paleturquoise")
    # Cavern Ceiling
    ax.plot([self.CavernZ[0], self.CavernZ[1]], [self.archRadius+self.centreOfCurvature["y"], self.archRadius+self.centreOfCurvature["y"]], 
                              c="paleturquoise")
    # Access Shafts
    ax.plot([self.PX14_Centre["z"]-self.PX14_Radius, self.PX14_Centre["z"]-self.PX14_Radius],
            [self.PX14_Centre["y"], self.PX14_Centre["y"]+self.PX14_Height], c="red", label="PX14", alpha=0.5)
    ax.plot([self.PX14_Centre["z"]+self.PX14_Radius,self.PX14_Centre["z"]+self.PX14_Radius],
            [self.PX14_Centre["y"], self.PX14_Centre["y"]+self.PX14_Height], c="red", label="PX14", alpha=0.5)
    ax.plot([self.PX16_Centre["z"]-self.PX16_Radius,self.PX16_Centre["z"]-self.PX16_Radius],
            [self.PX16_Centre["y"], self.PX16_Centre["y"]+self.PX16_Height], c="turquoise", label="PX16", alpha=0.5)
    ax.plot([self.PX16_Centre["z"]+self.PX16_Radius,self.PX16_Centre["z"]+self.PX16_Radius],
            [self.PX16_Centre["y"], self.PX16_Centre["y"]+self.PX16_Height], c="turquoise", label="PX16", alpha=0.5)

    # Mark the Cavern Centre, IP, and Centre of Curvature for the ceiling
    ax.scatter(0, 0, c="r", marker = "x", label="Cavern Centre")
    ax.annotate("Centre", (0,0))
    ax.scatter(self.IP["z"], self.IP["y"], c="g", marker = "o", label="IP")
    ax.annotate("IP", (self.IP["z"], self.IP["y"]))
    ax.plot( [self.CavernZ[0], self.CavernZ[1]], [self.centreOfCurvature["y"]]*2, c="b", label="Ceiling Centre of Curvature")
    ax.annotate("Centre of Curvature (Ceiling)", (self.CavernZ[0]/2, self.centreOfCurvature["y"]))

    if plotATLAS:
        ax.plot( [self.ATLAS_Z[0],self.ATLAS_Z[0]], [self.ATLAS_Centre["y"]-self.radiusATLAS,self.ATLAS_Centre["y"]+self.radiusATLAS], c="b")
        ax.plot( [self.ATLAS_Z[0],self.ATLAS_Z[1]], [self.ATLAS_Centre["y"]-self.radiusATLAS,self.ATLAS_Centre["y"]-self.radiusATLAS], c="b")
        ax.plot( [self.ATLAS_Z[1],self.ATLAS_Z[1]], [self.ATLAS_Centre["y"]-self.radiusATLAS,self.ATLAS_Centre["y"]+self.radiusATLAS], c="b")
        ax.plot( [self.ATLAS_Z[1],self.ATLAS_Z[1]], [self.ATLAS_Centre["y"]+self.radiusATLAS,self.ATLAS_Centre["y"]+self.radiusATLAS], c="b")

    if plotAcceptance:
        # Plot a rough impression of the Acceptance
        ax.plot([self.IP["z"], self.CavernZ[0]], [self.IP["y"], self.CavernY[1]], c="k", alpha=0.25, linestyle="--")
        ax.plot([self.IP["z"], self.CavernZ[1]], [self.IP["y"], self.CavernY[1]], c="k", alpha=0.25, linestyle="--")

    ax.set_xlim(-30,30)
    ax.set_ylim(-15,25)

def plotCavern3D(self, ax, plotATLAS=False, plotAcceptance=False): 
    # Get the Cavern ceiling data grid
    cavernArch = self.createCavernVault(doPlot=False)
    # Get the Access Shafts data grid 
    accessShafts = self.createAccessShafts()
    # Cavern Boundaries
    cavernBounds = { "x": np.linspace(self.CavernX[0], self.CavernX[1],100),
                     "y": np.linspace(self.CavernY[0], self.CavernY[1],100),
                     "z": np.linspace(self.CavernZ[0], self.CavernZ[1],100),}

    #3D Plot
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # Cavern Boundaries
    cavX_Y, cavY_X = np.meshgrid(cavernBounds["x"], cavernBounds["y"])
    cavY_Z, cavZ_Y = np.meshgrid(cavernBounds["y"], cavernBounds["z"])
    cavX_Z, cavZ_X = np.meshgrid(cavernBounds["x"], cavernBounds["z"])
    # XY Faces
    ax.plot_surface(cavX_Y, self.CavernZ[0]*np.ones(cavX_Y.shape), cavY_X, rstride=4, cstride=4, alpha=0.25)
    ax.plot_surface(cavX_Y, self.CavernZ[1]*np.ones(cavX_Y.shape), cavY_X, rstride=4, cstride=4, alpha=0.25)
    # YZ Faces
    ax.plot_surface(self.CavernX[0]*np.ones(cavY_Z.shape), cavZ_Y, cavY_Z, rstride=4, cstride=4, alpha=0.25)
    ax.plot_surface(self.CavernX[1]*np.ones(cavY_Z.shape), cavZ_Y, cavY_Z, rstride=4, cstride=4, alpha=0.25)
    # XZ Face
    ax.plot_surface(cavX_Z, cavZ_X, self.CavernY[0]*np.ones(cavX_Z.shape), rstride=4, cstride=4, alpha=0.25)
    # IP
    ax.scatter([self.IP["x"]], [self.IP["z"]], [self.IP["y"]], c="g", marker = "o", label="IP") 

    # Access Shafts
    for shaft in ["PX14", "PX16"]:
        ax.plot_surface(accessShafts[shaft]["x"], accessShafts[shaft]["z"], accessShafts[shaft]["y"], rstride=4, cstride=4, alpha=0.25, label=shaft)

    # Ceiling
    xx, zz = np.meshgrid(cavernArch["x"], cavernArch["z"])

    mask = np.logical_and(np.sqrt((xx-self.PX14_Centre["x"])**2 + (zz-self.PX14_Centre["z"])**2) > self.PX14_Radius, 
                            np.sqrt((xx-self.PX16_Centre["x"])**2 + (zz-self.PX16_Centre["z"])**2) > self.PX16_Radius)

    xx[~mask] = np.nan
    zz[~mask] = np.nan
    yy = np.sqrt(np.power(self.archRadius,2) - np.power(xx-self.centreOfCurvature["x"],2)) + self.centreOfCurvature["y"]

    ax.plot_surface(xx, zz, yy, rstride=4, cstride=4, alpha=0.25)

    if plotATLAS:
        atlasX =  np.linspace(self.ATLAS_Centre["x"]-self.radiusATLAS,self.ATLAS_Centre["x"]+self.radiusATLAS, 100)
        atlasZ =  np.linspace(self.ATLAS_Z[0], self.ATLAS_Z[1], 100)
        
        atlasXX, atlasZZ = np.meshgrid(atlasX, atlasZ)
        atlasYY = np.sqrt(np.power(self.radiusATLAS,2) - np.power(atlasXX-self.ATLAS_Centre["x"],2))
        ax.plot_surface(atlasXX, atlasZZ, atlasYY + self.ATLAS_Centre["y"], rstride=4, cstride=4, alpha=0.4, color="gray")
        ax.plot_surface(atlasXX, atlasZZ, -atlasYY + self.ATLAS_Centre["y"], rstride=4, cstride=4, alpha=0.4, color="gray")

    if plotAcceptance:
        # Plot a rough impression of the Acceptance
        ax.plot([self.IP["x"], self.CavernX[0]], [self.IP["z"], self.CavernZ[0]], [self.IP["y"],self.obtainCavernYFromX(self.CavernX[0])], 
                 c="k", alpha=0.25, linestyle="--")
        ax.plot([self.IP["x"], self.CavernX[0]], [self.IP["z"], self.CavernZ[1]], [self.IP["y"],self.obtainCavernYFromX(self.CavernX[0])], 
                 c="k", alpha=0.25, linestyle="--")
        ax.plot([self.IP["x"], self.CavernX[1]], [self.IP["z"], self.CavernZ[0]], [self.IP["y"], self.obtainCavernYFromX(self.CavernX[1])], 
                c="k", alpha=0.25, linestyle="--")
        ax.plot([self.IP["x"], self.CavernX[1]], [self.IP["z"], self.CavernZ[1]], [self.IP["y"],self.obtainCavernYFromX(self.CavernX[1])], 
                c="k", alpha=0.25, linestyle="--")

    ax.set_zlabel("y /m")
    ax.set_xlim(-30,30)
    ax.set_ylim(-30,30)
    ax.set_zlim(-30,30)

#---------------------------------------#
#- Set of functions to plot RPC Layers -#
#---------------------------------------#
# ANUBISrpcs is a list of RPClayers, each contain a list of RPCs in the format:
#   - "corners": 8 (x,y,z) coordinates corresponding to their corners,
#   - "midPoint": The midPoint of the RPC in (x,y,z), 
#   - "LayerID" and "RPCid": A Layer ID and RPC ID to uniquely identify the RPC
#   - "plane": A Sympy plane in the eta-phi plane that passes through the midpoint
def plotRPCsXY(self, ax, ANUBISrpcs):
    LayerColours = ["violet", "coral", "green", "maroon", "chartreuse", "gray"]

    nLayer=0
    for rpcLayer in ANUBISrpcs:
        tempRPCList = self.convertRPCList(rpcLayer)
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

def plotRPCsXZ(self, ax, ANUBISrpcs):
    LayerColours = ["violet", "coral", "green", "maroon", "chartreuse", "gray"]

    nLayer=0
    for rpcLayer in ANUBISrpcs:
        if nLayer!=0:
            continue # For clarity only plot the first RPC layer
        tempRPCList = self.convertRPCList(rpcLayer)
        # Plot midpoints
        ax.scatter([x[0] for x in tempRPCList["midPoint"]], [z[2] for z in tempRPCList["midPoint"]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
        # Plot RPC Boundaries
        for i in range(len(tempRPCList["corners"])):
            c = tempRPCList["corners"][i]
            # For the Top of the RPC: use 4--5, 5--7, 4--6 & 6--7 for Bottom of RPC.
            ax.plot( [c[0][0], c[1][0]], [c[0][2],c[1][2]], c=LayerColours[nLayer], label=f"Layer {nLayer}")

def plotRPCsZY(self, ax, ANUBISrpcs):
    LayerColours = ["violet", "coral", "green", "maroon", "chartreuse", "gray"]

    nLayer=0
    for rpcLayer in ANUBISrpcs:
        if nLayer!=0:
            continue #Plot only the first layer for clarity
        tempRPCList = self.convertRPCList(rpcLayer)
        # Plot midpoints
        ax.scatter([z[2] for z in tempRPCList["midPoint"]], [y[1] for y in tempRPCList["midPoint"]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
        # Plot RPC Boundaries
        for i in range(len(tempRPCList["corners"])):
            c = tempRPCList["corners"][i]
            # For left x corner of the RPC: use 2--3, 2--6, 3--7 & 6--7 for Right of RPC.
            ax.plot( [c[0][2], c[1][2]], [c[0][1],c[1][1]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
            ax.plot( [c[0][2], c[4][2]], [c[0][1],c[4][1]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
            ax.plot( [c[1][2], c[5][2]], [c[1][1],c[5][1]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
            ax.plot( [c[4][2], c[5][2]], [c[4][1],c[5][1]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
        nLayer+=1

def plotRPCs3D(self, ax, ANUBISrpcs):
    LayerColours = ["violet", "coral", "green", "maroon", "chartreuse", "gray"]

    nLayer=0
    for rpcLayer in ANUBISrpcs:
        tempRPCList = self.convertRPCList(rpcLayer)
        # Plot midpoints
        ax.scatter([x[0] for x in tempRPCList["midPoint"]], [z[2] for z in tempRPCList["midPoint"]], [y[1] for y in tempRPCList["midPoint"]], 
                   c=LayerColours[nLayer], label=f"Layer {nLayer}")

        # Plot RPC Boundaries
        for i in range(len(tempRPCList["corners"])):
            c = tempRPCList["corners"][i]
            # For left x corner of the RPC: use 2--3, 2--6, 3--7 & 6--7 for Right of RPC.
            ax.plot_surface(c[:][0], c[:][2], c[:][1], color = LayerColours[nLayer], label=f"Layer {nLayer}", 
                            alpha=0.4, rstride=5, cstride=4)
        nLayer+=1

#----------------------------------------------#
#- Set of functions to plot Simple RPC Layers -#
#----------------------------------------------#
# Here ANUBISrpcs contain a radial distance and the angular coverage.
# This has the form of: RPCs={"r": [[minR, maxR],...], "theta": [[minTheta, maxTheta],...], "phi": [[minPhi, maxPhi],...]}
def plotSimpleRPCsXY(self, ax, ANUBISrpcs):
    LayerColours = ["violet", "coral", "green", "maroon", "chartreuse", "gray"]
    for idx in range(len(ANUBISrpcs["r"])):
        ax.add_patch( matplotlib.patches.Arc((self.centreOfCurvature["x"], self.centreOfCurvature["y"]),
                                            width = 2*(ANUBISrpcs["r"][idx][0]), height = 2*(ANUBISrpcs["r"][idx][0]), 
                                            theta1=min(ANUBISrpcs["theta"]["CoC"][idx])*(180/np.pi), theta2=max(ANUBISrpcs["theta"]["CoC"][idx])*(180/np.pi),
                                            color=LayerColours[idx], fill=False, ls="--") )
        ax.add_patch( matplotlib.patches.Arc((self.centreOfCurvature["x"], self.centreOfCurvature["y"]),
                                            width = 2*(ANUBISrpcs["r"][idx][1]), height = 2*(ANUBISrpcs["r"][idx][1]), 
                                            theta1=min(ANUBISrpcs["theta"]["CoC"][idx])*(180/np.pi), theta2=max(ANUBISrpcs["theta"]["CoC"][idx])*(180/np.pi),
                                            color=LayerColours[idx], fill=False, ls="--") )


#--------------------------------------------#
#-     Plot Hits as 2D or 3D histograms     -#
#--------------------------------------------#
def plotHitsHist(self, axis, hits, binDict={"rangeX": [-20,20], "rangeY": [-25,25], "nX": 80, "nY": 100}):
    counts, xedges, yedges, im = axis.hist2d(hits[0], hits[1], range=(binDict["rangeX"],binDict["rangeY"]), bins = (binDict["nX"], binDict["nY"]), cmin=1)
    cbar = plt.colorbar(im, fraction=0.046, pad=0.04, ax=axis)

def plotHitsScatter(self, ax, hits, styleDict={"colour": "k", "marker": "."}):
    ax.scatter(hits[0], hits[1], c=styleDict["colour"], marker=styleDict["marker"])
