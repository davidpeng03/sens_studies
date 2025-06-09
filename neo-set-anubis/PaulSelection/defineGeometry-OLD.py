import sys, os
import numpy as np
import matplotlib
#matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.ion()
from sympy import Point3D, Line3D, Plane
import json
import pickle
import helper as h

#=========================================================#
# NOTE: The ATLAS Coordinate system is assumed throughout #
#=========================================================#
#   A side of ATLAS = +ve z
#   C side of ATLAS = -ve z
#   proANUBIS is on side A.

class ATLASCavern():
    # Define the parameters for the ATLAS cavern using technical drawings:
    #   - LHCJUX150002: Axial View (xy) of the ATLAS Cavern
    #   - ATFIU___0004: Side View (zy) of the ATLAS Cavern with the Access Shafts

    def __init__(self):
        #=========================#
        #   Cavern Dimensions     #
        #=========================#
        # NOTE: On the ground floor of the cavern the length is exactly 30m, 
        #       however above 25.6m from the floor the cavern walls are thinner.
        #       Since ANUBIS Will be above this level, we take this larger size for our X Length
        self.CavernXLength = 31.0 #metres
        self.CavernX = [-self.CavernXLength/2, self.CavernXLength/2]

        # NOTE: Y here does not include the vaulted ceiling - this will be included separately.
        # Taken from Axial View of the Cavern.
        self.CavernYLength = 27.5 #metres
        self.CavernY = [-self.CavernYLength/2, self.CavernYLength/2]
        
        self.CavernZLength = 52.8 #metres  #Potentially 53m
        self.CavernZ = [-self.CavernZLength/2, self.CavernZLength/2]

        #=========================#
        #   Cavern Boundaries     #
        #=========================#
        self.CavernBounds = {"x": self.CavernX, "y": self.CavernY, "z": self.CavernZ}
        self.CavernCorners = {"X0Y0Z0": [self.CavernX[0],self.CavernY[0],self.CavernZ[0]],
                              "X0Y0Z1": [self.CavernX[0],self.CavernY[0],self.CavernZ[1]],
                              "X0Y1Z0": [self.CavernX[0],self.CavernY[1],self.CavernZ[0]],
                              "X0Y1Z1": [self.CavernX[0],self.CavernY[1],self.CavernZ[1]],
                              "X1Y0Z0": [self.CavernX[1],self.CavernY[0],self.CavernZ[0]],
                              "X1Y0Z1": [self.CavernX[1],self.CavernY[0],self.CavernZ[1]],
                              "X1Y1Z0": [self.CavernX[1],self.CavernY[1],self.CavernZ[0]],
                              "X1Y1Z1": [self.CavernX[1],self.CavernY[1],self.CavernZ[1]],
        }

        #=============================================================================#
        # The Interaction point does not align with the centre of the ATLAS Cavern.
        # Give the relative position to the centre of the ATLAS Cavern:
        self.IP = {"x": -1.7, #metres
                   "y": -self.CavernYLength/2 + 11.37, #metres #From https://core.ac.uk/download/pdf/44194071.pdf Page 14
                   "z": 0} #metres
        #=============================================================================#

        #===========================#
        #   Cavern Ceiling Arch     #
        #===========================#
        # The profile of the ATLAS Ceiling has been measured in: https://edms.cern.ch/document/2149688
        #   Gives the Equation of the cylinder to be: 20^2 = (x-1.7)^2 + (y-3.52)^2 relative to the IP
        #   Centre of Curvature of the ceiling doesn't align with the centre of the cavern 
        self.archRadius = 20 #metres, radius of curvature of the ceiling
        self.arcLength = self.archRadius * ( 2 * np.arcsin(self.CavernXLength / (2*self.archRadius)) )
        self.centreOfCurvature = {"x": 0, #metres, relative to cavern centre
                                  "y": self.IP["y"] + 3.52} #metres
        
        # Difference between the theoretical y level based on CavernYLength and that given by the ceiling profile relative to Cavern Centre 
        self.archOffset = self.CavernYLength/2 - np.sqrt(np.power(self.archRadius,2) - np.power(self.CavernX[0],2)) - self.centreOfCurvature["y"]
        
        #===========================#
        #  Service Shaft Parameters #
        #===========================#
        self.PX14_Centre = {"x": 0, #metres - aligns with cavern axis from edms 2149688
                            "y": self.centreOfCurvature["y"] + self.archRadius, #metres
                            "z": 13.5} #metres - From Axial View and edms 2149688

        self.PX14_Radius = 18/2 # metres
        self.PX14_Height = 57.85 # metres 
        self.PX14_LowestY = np.sqrt(np.power(self.archRadius,2) - (np.power(self.PX14_Radius - self.PX14_Centre["x"],2)))

        self.PX16_Centre = {"x": 0, #metres - aligns with cavern axis from edms 2149688
                            "y": self.centreOfCurvature["y"] + self.archRadius, #metres
                            "z": -17.7} #metres - From Axial View and edms 2149688

        self.PX16_Radius = 12.6/2 # metres
        self.PX16_Height = 57.85 # metres
        self.PX16_LowestY = np.sqrt(np.power(self.archRadius,2) - (np.power(self.PX16_Radius - self.PX16_Centre["x"],2)))

        self.angles = self.calculateAngles()

        #===========================#
        # ATLAS Experimental Bounds #
        #===========================#
        self.radiusATLAS = 12 # metres, ATLAS experiment envelope (Estimated from Fig 1.3 from https://cds.cern.ch/record/2285580 (Oleg had 10.5m)
        self.radiusATLAStracking = 7 # metres, ATLAS effective vertexing radius
        self.ATLAS_ZLength = 44 # metres, ATLAS experiment envelope (Estimated from Fig 1.3 from https://cds.cern.ch/record/2285580 
        self.ATLAS_Z = [-self.ATLAS_ZLength/2 - self.IP["z"], self.ATLAS_ZLength/2 - self.IP["z"]] #ATLAS Z min and max in metres relative to Cavern centre
        self.ATLAS_Centre = self.IP # ATLAS Centre coincides with the IP

    #------------------------------#
    #     Helper Functions         #
    #------------------------------#
    # Convert cartesian coordinates in terms of the Cavern Centre to be in terms of the IP instead
    #   - Useful as simulations are relative to the IP
    def cavernCentreToIP(self, x, y, z):
        return (x + self.IP["x"], y + self.IP["y"], z + self.IP["z"])

    # Convert cartesian coordinates in terms of the IP to be in terms of the Cavern Centre instead
    def IPTocavernCentre(self, x, y, z):
        return (x - self.IP["x"], y - self.IP["y"], z - self.IP["z"])

    # Convert cartestian coordinates to cylindrical coordinates, where the circular face is in xy
    def cartToCyl(self, x, y, z):
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        return (r, theta, z)

    # Convert cylindrical coordinates to cartestian coordinates, where the circular face is in xy
    def cylToCart(self, r, theta, z):
        return (r*np.cos(theta), r*np.sin(theta), z)

    # Convert spherical coordinates to cartestian coordinates
    def cartToSph(self, x, y, z):
        r = np.sqrt( (x**2) + (y**2) + (z**2) )
        theta = np.arccos( z / r )
        phi = np.sign(y)*np.arccos(x / (np.sqrt( (x**2) + (y**2) ) ) )
        return (r, theta, phi)

    def sphToCart(self, r, theta, phi):
        x = r*np.sin(theta)*np.cos(phi)
        y = r*np.sin(theta)*np.sin(phi)
        z = r*np.cos(theta)
        return (x,y,z)

    def onCeiling(self, x, y, z):
        ceiling = np.power(x-self.centreOfCurvature["x"],2) + np.power(y-self.centreOFCurvature["y"],2)
        if (abs(self.archRadius - ceiling) < 0.01) and  (z > self.CavernZ[0] and z < self.CavernZ[1]):
            return (True, ceiling)
        else:
            return (False, ceiling)

    # Get the y coordinate of the cavern boundaries for a given x coordinate
    def obtainCavernYFromX(self, x): # Relative to Cavern Centre
        return self.centreOfCurvature["y"] + np.sqrt(self.archRadius**2 - np.power( (x-self.centreOfCurvature["x"]),2))

    # Get the x coordinate of the cavern boundaries for a given y coordinate
    def obtainCavernXFromY(self, y): # Relative to Cavern Centre
        return self.centreOfCurvature["x"] + np.sqrt(self.archRadius**2 - np.power( (y-self.centreOfCurvature["y"]),2))

    # Each point is a 3D position list
    def createVector(self, point1, point2):
        diff = np.array(point1) - np.array(point2)
        return diff / np.linalg.norm(diff)

    #============================================================#
    # Function to Determine if a point lies in the ATLAS Cavern  #
    #============================================================#
    def inCavern(self, x, y, z, maxRadius="", origin=[]):
        # Radius is assumed to be from the given origin (centre of curvature unless otherwise stated)
        # changing Radius to be IP point
        if len(origin)==0:
            origin = (self.centreOfCurvature["x"], self.centreOfCurvature["y"], 0)
            #origin = (self.IP["x"], self.IP["y"], self.IP["z"])
        if maxRadius=="":
            radialAcceptance=True # No radius considered

        else:
            # Radial distance in cylindrical coordinates relative to the given origin
            r = np.linalg.norm( (x - origin[0], y - origin[1]) )
            radialAcceptance = r < maxRadius

        #Assume x, y, z provided relative to the Cavern Centre
        if ((x > self.CavernX[0] and x < self.CavernX[1]) and #Within X Dimensions
            (y > self.CavernY[0] and y < self.obtainCavernYFromX(x)) and #Within Y Dimensions
            (z > self.CavernZ[0] and z < self.CavernZ[1]) and #Within Z Dimensions
            (radialAcceptance ) # Within a Radius if specified
            ):
            return True
        else:
            return False

    def inATLAS(self, x, y, z, trackingOnly=False, origin=[], verbose=False):
        #Assume x, y, z provided relative to the Cavern Centre otherwise use the given origin
        if len(origin)==0:
            origin = (self.IP["x"], self.IP["y"], self.IP["z"])

        # Radial distance in cylindrical coordinates relative to the given origin
        #   - ATLAS is defined as a cylinder with a defined radius here.
        r = np.linalg.norm( (x - origin[0], y - origin[1]) )

        if trackingOnly:
            rTarget = self.radiusATLAStracking #Effective vertexing radius
        else:
            rTarget = self.radiusATLAS # Entire ATLAS Envelope

        if ( (r < rTarget) and (z > self.ATLAS_Z[0] and z < self.ATLAS_Z[1]) ):
            return True
        else:
            return False

    def intersectANUBISstations(self, x, y, z, ANUBISstations, origin=[]):
        # (x,y,z) is the position of a particle
        # ANUBISstations is a dictionary of RPCs, with a list of: 
        #   - "corners": The (x,y,z) positions of its 8 corners,
        #   - "midPoint": The (x,y,z) position of the midpoint,
        #   - "plane": A sympy Plane Object which goes through the midpoint
        #   - "LayerID" and "RPCid": which combine to form a unique ID to identify it
        # Origin allows the origin of the particle (x,y,z) to be set to determine its direction. If empty, assumed to originate from IP
        
        if len(origin)==0:
            origin = (self.IP["x"], self.IP["y"], self.IP["z"])

        coordSphere = self.cartToSph(x, y, z)
        
        direction = Line3D(Point3D(origin), Point3D((x,y,z)))

        nIntersections=0
        intersections=[]
        for i in range(len(ANUBISstations["plane"])):
            # Reduce function calls by only checking cases where the hit is within an angular separation of 0.1 of the RPC midpoint
            midSphere = self.cartToSph(ANUBISstations["midPoint"][i][0], ANUBISstations["midPoint"][i][1], ANUBISstations["midPoint"][i][2])
            if np.sqrt((midSphere[1]-coordSphere[1])**2 + (midSphere[2]-coordSphere[2])**2) > 0.1:
                continue

            tempIntersections = ANUBISstations["plane"][i].intersection(direction)

            if len(tempIntersections)!=0:
                #nIntersections+=len(tempIntersections) # Increment the number of intersections
                for intersect in tempIntersections:
                    if ((ANUBISstations["corners"][i][0][0] > intersect[0] and intersect[0] <  ANUBISstations["corners"][i][0][0]) and
                        (ANUBISstations["corners"][i][4][1] > intersect[1] and intersect[1] <  ANUBISstations["corners"][i][1][0]) and
                        (ANUBISstations["corners"][i][1][2] > intersect[2] and intersect[2] <  ANUBISstations["corners"][i][2][0])
                        ):
                        intersections.append((intersect[0], intersect[1], intersect[2])) # Store the Intersection points
                        nIntersections+=1

        return nIntersections, intersections
    #select lps uses
    # phi - x y angle
    # eta = y z angle
    #
    #
    #
    def intersectANUBISstationsSimple(self, x, y, z, ANUBISstations, origin=[], verbose=False):
        # (x,y,z) is the position of a particle
        # Origin allows the origin of the particle (x,y,z) to be set to determine its direction. If empty, assumed to originate from IP
        # ANUBISstations is a dictionary of a list of RPC layers represented by spherical shells, containing:
        #   - "r": List of [minRadius, maxRadius] (Representing the thickness of an RPC layer) relative to the centre of curvature.
        #   - "theta": List of [minTheta, maxTheta] (Representing the start and end of the chambers in angular coverage.
        #   - "phi": List of [minPhi, maxPhi]
        
        if len(origin)==0:
            origin = (self.IP["x"], self.IP["y"], self.IP["z"])
            #x,y,z = self.cavernCentreToIP(x, y, z)

        r, theta, phi = self.cartToSph(x - origin[0], y - origin[1], z - origin[2])
        #direction = Line3D(Point3D(origin), Point3D((x,y,z)))
        # Get the radial distance relative to the centre of curvature.
        radialDist = np.linalg.norm( ( x - self.centreOfCurvature["x"], y - self.centreOfCurvature["y"]))

        if verbose:
            print(f"(X,Y,Z): ({x},{y},{z}), (r,theta,phi): ({r},{theta},{phi})")
        
        nIntersections=0
        intersections=[]
        for i in range(len(ANUBISstations["r"])):
            if verbose:
                print(f"Station {i}: (r,theta,phi): ({ANUBISstations['r'][i]},{ANUBISstations['phi'][i]},{ANUBISstations['theta'][i]})")
                print(f"\t r Requirement: {radialDist > ANUBISstations['r'][i][1]}")
                print(f"\t phi Requirement: {( (ANUBISstations['phi'][i][0] < phi) and (phi < ANUBISstations['phi'][i][1]))}")
                print(f"\t theta Requirement: {((ANUBISstations['theta'][i][0] < theta) and (theta < ANUBISstations['theta'][i][1]))}")

            if radialDist > ANUBISstations["r"][i][1]: # ANUBIS r limits are in cylindrical coordinates
                if verbose:
                    print(f"Skipping Station...")
                continue
            if ( ( (ANUBISstations["phi"][i][0] < phi) and (phi < ANUBISstations["phi"][i][1]) ) and
                 ( (ANUBISstations["theta"][i][0] < theta) and (theta < ANUBISstations["theta"][i][1]) )):
                 intersections.append(self.sphToCart(ANUBISstations["r"][i][0],theta,phi))
                 nIntersections+=1

        return nIntersections, intersections

    def intersectANUBISstationsSimple2(self,x,y,z, r, theta, phi, ANUBISstations, origin=[], verbose=False):
        # (x,y,z) is the position of a particle
        # Origin allows the origin of the particle (x,y,z) to be set to determine its direction. If empty, assumed to originate from IP
        # ANUBISstations is a dictionary of a list of RPC layers represented by spherical shells, containing:
        #   - "r": List of [minRadius, maxRadius] (Representing the thickness of an RPC layer) relative to the centre of curvature.
        #   - "theta": List of [minTheta, maxTheta] (Representing the start and end of the chambers in angular coverage.
        #   - "phi": List of [minPhi, maxPhi]
        
        if len(origin)==0:
            origin = (self.IP["x"], self.IP["y"], self.IP["z"])
            #x,y,z = self.cavernCentreToIP(x, y, z)

        #r, theta, phi = self.cartToSph(x - origin[0], y - origin[1], z - origin[2])
        #direction = Line3D(Point3D(origin), Point3D((x,y,z)))
        # Get the radial distance relative to the centre of curvature.
        #radialDist = np.linalg.norm( ( x - self.centreOfCurvature["x"], y - self.centreOfCurvature["y"]))
        radialDist = np.linalg.norm( (x - origin[0], y - origin[1]))

        if verbose:
            print(f"(X,Y,Z): ({x},{y},{z}), (r,theta,phi): ({r},{theta},{phi})")
        
        nIntersections=0
        intersections=[]
        for i in range(len(ANUBISstations["r"])):
            if verbose:
                print(f"Station {i}: (r,theta,phi): ({ANUBISstations['r'][i]},{ANUBISstations['phi'][i]},{ANUBISstations['theta'][i]})")
                print(f"\t r Requirement: {radialDist > ANUBISstations['r'][i][1]}")
                print(f"\t phi Requirement: {( (ANUBISstations['phi'][i][0] < phi) and (phi < ANUBISstations['phi'][i][1]))}")
                print(f"\t theta Requirement: {((ANUBISstations['theta'][i][0] < theta) and (theta < ANUBISstations['theta'][i][1]))}")

            if radialDist > ANUBISstations["r"][i][1]: # ANUBIS r limits are in cylindrical coordinates
                if verbose:
                    print(f"Skipping Station...")
                continue
            if ( ( (ANUBISstations["phi"][i][0] < phi) and (phi < ANUBISstations["phi"][i][1]) ) and
                 ( (ANUBISstations["theta"][i][0] < theta) and (theta < ANUBISstations["theta"][i][1]) )):
                 intersections.append(self.sphToCart(ANUBISstations["r"][i][0],theta,phi))
                 nIntersections+=1

        return nIntersections, intersections
    
    def SolidAngle(self, a, b, d):
        # Solid Angle of a rectangular Pyramid (See https://vixra.org/pdf/2001.0603v2.pdf, equation 27)
        alpha = a / (2*d)
        beta = b / (2*d)
        return 4*np.arctan( (alpha*beta) / np.sqrt(1 + alpha**2 + beta**2)) #sr

    def calculateAngles(self, relToCentre=False):
        # All Angles calculated relative to the IP unless relToCentre specified, then done relative to Cavern Centre
        if relToCentre:
            refX, refY, refZ = self.centreOfCurvature["x"], self.centreOfCurvature["y"], 0
            print("Angles are being calculated relative to the Cavern Centre...")
        else:
            refX, refY, refZ = self.IP["x"], self.IP["y"], self.IP["z"]
            print("Angles are being calculated relative to the IP...")

        # Phi in ATLAS coords (XY) relative to the x axis, phi=0 at y=0, +ve x and increases anti-clockwise
        phi=[]
        for i in [0,1]:
            # Corner of the cavern in XY space where Z is set to 0
            vec1 = self.createVector([self.CavernX[i],self.CavernY[1]], [refX, refY])
            vec2 = self.createVector([self.CavernX[1], refY], [refX, refY])
            tempPhi = np.dot(vec1, vec2) / ( np.linalg.norm(vec1) * np.linalg.norm(vec2))
            phi.append(np.arccos(np.clip(tempPhi, -1, 1))) # In radians

        # Polar Angle in ATLAS coords (YZ), 0 when aligned with +ve z-axis
        theta=[]
        for i in [0,1]:
            # Corner of the cavern in YZ space where X is set to 0
            vec1 = self.createVector([self.CavernZ[i],self.obtainCavernYFromX(0)] , [refZ, refY])
            vec2 = self.createVector([self.CavernZ[1], refY], [refZ, refY])
            tempTheta = np.dot(vec1, vec2) / ( np.linalg.norm(vec1) * np.linalg.norm(vec2))
            theta.append(np.arccos(np.clip(tempTheta, -1, 1))) # In radians

        # Eta in ATLAS coords (YZ) 
        eta = [ (x/abs(x)) * -np.log(np.tan(abs(x)/2)) for x in theta]

        # Solid Angle
        a = self.CavernZLength
        b = self.CavernXLength
        d = self.archRadius + self.centreOfCurvature["y"]
        # On Axis case: solidAngle = self.SolidAngle(a,b,d)
        print(f"On Axis SolidAngle: {self.SolidAngle(a,b,d)}")
        # IP is Off-axis of the centre of the projected rectangle of the ceiling, 
        # use A and B to scale as in https://vixra.org/pdf/2001.0603v2.pdf, equation 34
        A = self.CavernZLength/2 
        B = abs(refX)
        solidAngle = (self.SolidAngle(2*(a-A),2*(b-B),d) + self.SolidAngle(2*A,2*(b-B),d) + self.SolidAngle(2*(a-A),2*B,d) + self.SolidAngle(2*A,2*B,d) ) /4

        return {"phi": phi, "phiRange": abs(phi[0])+abs(phi[1]),\
                "theta": theta, "thetaRange": abs(theta[1] - theta[0]),\
                "eta": eta, "etaRange": abs(eta[0])+abs(eta[1]),\
                "solidAngle": solidAngle} 

    #--------------------------#
    #   Plotting Functions     #
    #--------------------------#
    # Create a grid of points within the Access Shafts 
    def createAccessShafts(self):
        theta = np.linspace(0, 2*np.pi, 100)
        y_PX14 = np.linspace(self.PX14_Centre["y"],self.PX14_Centre["y"]+self.PX14_Height,100)
        gridTheta_PX14, gridY_PX14 = np.meshgrid(theta, y_PX14)
        gridX_PX14 = self.PX14_Radius*np.cos(gridTheta_PX14) + self.PX14_Centre["x"]
        gridZ_PX14 = self.PX14_Radius*np.sin(gridTheta_PX14) + self.PX14_Centre["z"]

        y_PX16 = np.linspace(self.PX16_Centre["y"],self.PX16_Centre["y"]+self.PX16_Height,100)
        gridTheta_PX16, gridY_PX16 = np.meshgrid(theta, y_PX16)
        gridX_PX16 = self.PX16_Radius*np.cos(gridTheta_PX16) + self.PX16_Centre["x"]
        gridZ_PX16 = self.PX16_Radius*np.sin(gridTheta_PX16) + self.PX16_Centre["z"]

        return { "PX14": {"x": gridX_PX14, "y": gridY_PX14, "z": gridZ_PX14},\
                 "PX16": {"x": gridX_PX16, "y": gridY_PX16, "z": gridZ_PX16} }

    # Plot the Cavern Ceiling -- deprecated with plotFullATLASCavern
    def createCavernVault(self, doPlot=True):
        archX =  np.linspace(self.CavernX[0], self.CavernX[1], 100)
        archY = []
        for x in archX:
            archY.append(np.sqrt(np.power(self.archRadius,2) - np.power(x,2)) + self.centreOfCurvature["y"])
        
        archZ =  np.linspace(self.CavernZ[0], self.CavernZ[1], 100)

        if doPlot:
            #XY at Z=0
            fig, ax = plt.subplots(1, figsize=(16, 10), dpi=100)
            ax.scatter(archX, archY, c="paleturquoise")
            plt.xlabel("x /m")
            plt.ylabel("y /m")
            plt.title("ATLAS Cavern Ceiling at z=0m")
            plt.tight_layout()
            plt.savefig("cavernCeiling_XY.pdf")
            plt.close()

            #XZ at Y = CavernYLength/2 
            fig, ax = plt.subplots(1, figsize=(16, 10), dpi=100)
            ax.scatter(archX, archZ, c="paleturquoise")
            plt.Circle((self.PX14_Centre["x"], self.PX14_Centre["z"]), self.PX14_Radius, color="blue", fill=False)
            ax.annotate("PX14", xy = (self.PX14_Centre["x"], self.PX14_Centre["z"]), fontsize=20, ha="center")
            plt.Circle((self.PX16_Centre["x"], self.PX16_Centre["z"]), self.PX16_Radius, color="blue", fill=False)
            ax.annotate("PX16", xy = (self.PX16_Centre["x"], self.PX16_Centre["z"]), fontsize=20, ha="center")
            plt.xlabel("x /m")
            plt.ylabel("z /m")
            plt.title(f"ATLAS Cavern Ceiling at y={self.CavernYLength/2}m")
            plt.tight_layout()
            plt.savefig("cavernCeiling_XZ.pdf")
            plt.close()

            #3D Plot
            fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
            
            xx, zz = np.meshgrid(archX, archZ)

            mask = np.logical_and(np.sqrt((xx-self.PX14_Centre["x"])**2 + (zz-self.PX14_Centre["z"])**2) > self.PX14_Radius, 
                                    np.sqrt((xx-self.PX16_Centre["x"])**2 + (zz-self.PX16_Centre["z"])**2) > self.PX16_Radius)

            xx[~mask] = np.nan
            zz[~mask] = np.nan
            yy = np.sqrt(np.power(self.archRadius,2) - np.power(xx-self.centreOfCurvature["x"],2))
            #yy[mask] = np.nan
            
            ax.plot_surface(xx, zz, yy, rstride=4, cstride=4, alpha=0.25)
            plt.xlabel("x /m")
            plt.ylabel("z /m")
            ax.set_zlabel("y /m")
            plt.tight_layout()
            plt.savefig("cavernCeiling_XYZ.pdf")
            plt.close()

        return {"x": archX, "y": archY, "z": archZ}

    def convertRPCList(self, rpcList):
        # Convert a list of RPCs with the following info:
        #   - corners: 8 (x,y,z) coordinates corresponding to their corners, 
        #   - midPoint: The midPoint of the RPC in (x,y,z), 
        #   - "LayerID" and "RPCid": A Layer ID and RPC ID to uniquely identify the RPC
        #   - "plane": A Sympy plane in the eta-phi plane that passes through the midpoint
        # To a total set of lists for each entry
        
        corners, midPoints, layerIDs, RPCIDs, planes = [], [], [], [], []
        for rpc in rpcList:
            corners.append(rpc["corners"])
            midPoints.append(rpc["midPoint"])
            layerIDs.append(rpc["LayerID"])
            RPCIDs.append(rpc["RPCid"])
            planes.append(rpc["plane"])

        return {"corners": corners, "midPoint": midPoints, "LayerID": layerIDs, "RPCid": RPCIDs, "plane": planes} 

    # Plot all features of the ATLAS Cavern, plus additional features if provided: e.g. ANUBIS.
    def plotFullCavern(self, hits={}, ANUBISrpcs=[], plotRPCs={"xy": True, "xz": False, "zy": False, "3D": False},
                             plotFailed=True, plotATLAS=False, suffix="", outDir="/usera/dp728/run_dir/plots"):
        # ANUBISrpcs is a list of RPClayers, each contain a list of RPCs in the format:
        #   - "corners": 8 (x,y,z) coordinates corresponding to their corners,
        #   - "midPoint": The midPoint of the RPC in (x,y,z), 
        #   - "LayerID" and "RPCid": A Layer ID and RPC ID to uniquely identify the RPC
        #   - "plane": A Sympy plane in the eta-phi plane that passes through the midpoint
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        # Get the Cavern ceiling data grid
        cavernArch = self.createCavernVault(doPlot=False)
        # Get the Access Shafts data grid 
        accessShafts = self.createAccessShafts()
        # Cavern Boundaries
        cavernBounds = { "x": np.linspace(self.CavernX[0], self.CavernX[1],100),
                         "y": np.linspace(self.CavernY[0], self.CavernY[1],100),
                         "z": np.linspace(self.CavernZ[0], self.CavernZ[1],100),}

        tempSuffix= { "xy": suffix, "xz": suffix, "zy": suffix, "3D": suffix}

        if len(ANUBISrpcs)!=0:
            LayerColours = ["violet", "coral", "green", "maroon", "chartreuse", "gray"]
            
            for label, decision in plotRPCs.items():
                if decision:
                    tempSuffix[label]+=f"_withANUBIS_{len(ANUBISrpcs)}Layers"

        #XY
        from matplotlib.colors import LogNorm

        fig, ax = plt.subplots(1, figsize=(10, 16), dpi=100)
        # Plot x,y positions of hits that pass and fail geometric acceptance.
        if len(hits)!=0:
            if plotFailed:
                ax.scatter([x[0] for x in hits["failed"]], [y[1] for y in hits["failed"]], c="red", marker="x")
            #ax.scatter([x[0] for x in hits["passed"]], [y[1] for y in hits["passed"]], c="lime", marker="^")
            counts, xedges, yedges, im = ax.hist2d([x[0] for x in hits["passed"]], [y[1] for y in hits["passed"]], 
                                                   range=(hits["bins"]["rangeX"],hits["bins"]["rangeY"]), bins = (hits["bins"]["x"], hits["bins"]["y"]), cmin=1,
        norm=LogNorm(vmin=1, vmax=1000) )
            cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
            # Plot a rough impression of the Acceptance
            ax.plot([self.IP["x"], self.CavernX[0]], [self.IP["y"], self.CavernY[1]], c="k", alpha=0.25, linestyle="--")
            ax.plot([self.IP["x"], self.CavernX[1]], [self.IP["y"], self.CavernY[1]], c="k", alpha=0.25, linestyle="--")

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
        
        if len(ANUBISrpcs)!=0 and plotRPCs["xy"]: #Include ANUBIS positions
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
            
            # Naive envelopes
            ax.add_patch( matplotlib.patches.Arc((self.centreOfCurvature["x"], self.centreOfCurvature["y"]),
                                                width = 2*(self.archRadius), height = 2*(self.archRadius), 
                                                theta1=min(self.angles["theta"])*(180/np.pi), theta2=max(self.angles["theta"])*(180/np.pi),
                                                color=LayerColours[0], fill=False, ls="--") )
            ax.add_patch( matplotlib.patches.Arc((self.centreOfCurvature["x"], self.centreOfCurvature["y"]),
                                                width = 2*(self.archRadius-0.06), height = 2*(self.archRadius-0.06),
                                                theta1=min(self.angles["theta"])*(180/np.pi), theta2=max(self.angles["theta"])*(180/np.pi),
                                                color=LayerColours[0], fill=False, ls="--") )

            ax.add_patch( matplotlib.patches.Arc((self.centreOfCurvature["x"], self.centreOfCurvature["y"]),
                                                width = 2*(self.archRadius-0.4), height = 2*(self.archRadius-0.4),
                                                theta1=min(self.angles["theta"])*(180/np.pi), theta2=max(self.angles["theta"])*(180/np.pi),
                                                color=LayerColours[1], fill=False, ls="--") )
            
            ax.add_patch( matplotlib.patches.Arc((self.centreOfCurvature["x"], self.centreOfCurvature["y"]),
                                                width =2*(self.archRadius-0.4-0.06), height = 2*(self.archRadius-0.4-0.06),
                                                theta1=min(self.angles["theta"])*(180/np.pi), theta2=max(self.angles["theta"])*(180/np.pi),
                                                color=LayerColours[1], fill=False, ls="--") )

            ax.add_patch( matplotlib.patches.Arc((self.centreOfCurvature["x"], self.centreOfCurvature["y"]),width = 2*(self.archRadius-1), height = 2*(self.archRadius-1),
                                                theta1=min(self.angles["theta"])*(180/np.pi), theta2=max(self.angles["theta"])*(180/np.pi),
                                                color=LayerColours[2], fill=False, ls="--") )
            ax.add_patch( matplotlib.patches.Arc((self.centreOfCurvature["x"], self.centreOfCurvature["y"]),
                                                width = 2*(self.archRadius-1-0.06), height = 2*(self.archRadius-1-0.06),
                                                theta1=min(self.angles["theta"])*(180/np.pi), theta2=max(self.angles["theta"])*(180/np.pi),
                                                color=LayerColours[2], fill=False, ls="--") )

        # Plot the ATLAS experiment boundaries
        if plotATLAS:
            ax.add_patch( plt.Circle((self.ATLAS_Centre["x"], self.ATLAS_Centre["y"]), self.radiusATLAS, color="blue", fill=False, ls="--") )

        # Mark the Cavern Centre, IP, and Centre of Curvature for the ceiling
        ax.scatter(0, 0, c="r", marker = "x", label="Cavern Centre")
        ax.annotate("Centre", (0,0))
        ax.scatter(self.IP["x"], self.IP["y"], c="g", marker = "o", label="IP")
        ax.annotate("IP", (self.IP["x"], self.IP["y"]))
        ax.scatter(self.centreOfCurvature["x"], self.centreOfCurvature["y"], c="b", marker = "D", label="Ceiling Centre of Curvature")
        ax.annotate("Centre of Curvature (Ceiling)", (self.centreOfCurvature["x"], self.centreOfCurvature["y"]))

        plt.xlabel(f"x /m")
        plt.ylabel(f"y /m")
        plt.title("ATLAS Cavern")
        plt.savefig(f"{outDir}/ATLASCavern_XY_withShafts{tempSuffix['xy']}.pdf")
        plt.gca().set_aspect('equal')
        ax.set_ylim(top=24)
        plt.savefig(f"{outDir}/ATLASCavern_XY{tempSuffix['xy']}.pdf")
        ax.set_xlim(-16,-5)
        ax.set_ylim(10,20)
        plt.savefig(f"{outDir}/ATLASCavern_XY_Zoomed{tempSuffix['xy']}.pdf")
        plt.close()

        #XZ
        fig, ax = plt.subplots(1, figsize=(14, 14), dpi=100)
        # Plot x,z positions of hits that pass and fail geometric acceptance.
        if len(hits)!=0:
            if plotFailed:
                ax.scatter([x[0] for x in hits["failed"]], [z[2] for z in hits["failed"]], c="red", marker="x")
            counts, xedges, yedges, im = ax.hist2d([x[0] for x in hits["passed"]], [z[2] for z in hits["passed"]], 
                                                    range=(hits["bins"]["rangeX"],hits["bins"]["rangeZ"]), bins = (hits["bins"]["x"], hits["bins"]["z"]), cmin=1)
            plt.colorbar(im)
            # Plot a rough impression of the Acceptance
            ax.plot([self.IP["x"], self.CavernX[0]], [self.IP["z"], self.CavernZ[0]], c="k", alpha=0.25, linestyle="--")
            ax.plot([self.IP["x"], self.CavernX[1]], [self.IP["z"], self.CavernZ[1]], c="k", alpha=0.25, linestyle="--")
        # Cavern Boundaries
        ax.plot( [self.CavernX[0], self.CavernX[1]], [self.CavernZ[0], self.CavernZ[0]], c="paleturquoise")
        ax.plot( [self.CavernX[0], self.CavernX[1]], [self.CavernZ[1], self.CavernZ[1]], c="paleturquoise")
        ax.plot( [self.CavernX[0], self.CavernX[0]], [self.CavernZ[0], self.CavernZ[1]], c="paleturquoise")
        ax.plot( [self.CavernX[1], self.CavernX[1]], [self.CavernZ[0], self.CavernZ[1]], c="paleturquoise")
        # Access Shafts
        ax.add_patch( plt.Circle((self.PX14_Centre["x"], self.PX14_Centre["z"]), self.PX14_Radius, color="blue", fill=False, ls="--") )
        ax.add_patch( plt.Circle((self.PX16_Centre["x"], self.PX16_Centre["z"]), self.PX16_Radius, color="blue", fill=False, ls="--") )
        # Mark the Cavern Centre, IP, and Centre of Curvature for the ceiling
        ax.scatter(0, 0, c="r", marker = "x", label="Cavern Centre")
        ax.annotate("Centre", (0,0))
        ax.scatter(self.IP["x"], self.IP["z"], c="g", marker = "o", label="IP")
        ax.annotate("IP", (self.IP["x"], self.IP["z"]))
        ax.plot([self.centreOfCurvature["x"]]*2, [self.CavernZ[0], self.CavernZ[1]], c="b", label="Ceiling Centre of Curvature")
        ax.annotate("Centre of Curvature (Ceiling)", (self.centreOfCurvature["x"],self.CavernZ[0]/4))

        if len(ANUBISrpcs)!=0 and plotRPCs["xz"]: #Include ANUBIS positions
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
                    ax.plot( [c[0][0], c[2][0]], [c[0][2],c[2][2]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
                    ax.plot( [c[2][0], c[3][0]], [c[2][2],c[3][2]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
                    ax.plot( [c[1][0], c[3][0]], [c[1][2],c[3][2]], c=LayerColours[nLayer], label=f"Layer {nLayer}")
                nLayer+=1

        if plotATLAS: #xz coordinates for atlas
            ax.plot( [self.ATLAS_Centre["x"]-self.radiusATLAS,self.ATLAS_Centre["x"]+self.radiusATLAS], [self.ATLAS_Z[0],self.ATLAS_Z[0]], c="gray")
            ax.plot( [self.ATLAS_Centre["x"]-self.radiusATLAS,self.ATLAS_Centre["x"]-self.radiusATLAS], [self.ATLAS_Z[0],self.ATLAS_Z[1]], c="blue")
            ax.plot( [self.ATLAS_Centre["x"]-self.radiusATLAS,self.ATLAS_Centre["x"]+self.radiusATLAS], [self.ATLAS_Z[1],self.ATLAS_Z[1]], c="green")
            #ax.plot( [self.ATLAS_Centre["x"]-self.radiusATLAS,self.ATLAS_Centre["x"]-self.radiusATLAS], [self.ATLAS_Z[0],self.ATLAS_Z[1]], c="gray")
            ax.plot( [self.ATLAS_Centre["x"]+self.radiusATLAS,self.ATLAS_Centre["x"]+self.radiusATLAS], [self.ATLAS_Z[0],self.ATLAS_Z[1]], c="red")

        #ax.set_xlim(-25,25)
        #ax.set_ylim(-25,25)
        plt.xlabel(f"x /m")
        plt.ylabel(f"z /m")
        plt.title("ATLAS Cavern")
        plt.savefig(f"{outDir}/ATLASCavern_XZ{tempSuffix['xz']}.pdf")
        plt.close()

        #ZY
        fig, ax = plt.subplots(1, figsize=(16, 10), dpi=100)
        # Plot z,y positions of hits that pass and fail geometric acceptance.
        if len(hits)!=0:
            if plotFailed:
                ax.scatter([z[2] for z in hits["failed"]], [y[1] for y in hits["failed"]], c="red", marker="x")
            #ax.scatter([z[2] for z in hits["passed"]], [y[1] for y in hits["passed"]], c="lime", marker="^")
            counts, xedges, yedges, im = ax.hist2d([z[2] for z in hits["passed"]], [y[1] for y in hits["passed"]], 
                                                    range=(hits["bins"]["rangeZ"],hits["bins"]["rangeY"]), bins = (hits["bins"]["z"], hits["bins"]["y"]), cmin=1)
            plt.colorbar(im,ax=ax)
            # Plot a rough impression of the Acceptance
            ax.plot([self.IP["z"], self.CavernZ[0]], [self.IP["y"], self.CavernY[1]], c="k", alpha=0.25, linestyle="--")
            ax.plot([self.IP["z"], self.CavernZ[1]], [self.IP["y"], self.CavernY[1]], c="k", alpha=0.25, linestyle="--")
        # Cavern Boundaries
        ax.plot([self.CavernZ[0], self.CavernZ[1]], [self.CavernY[1], self.CavernY[1]], c="paleturquoise", ls="--")
        ax.plot([self.CavernZ[0], self.CavernZ[0]], [self.CavernY[0], self.archRadius+self.centreOfCurvature["y"]], c="paleturquoise")
        ax.plot([self.CavernZ[1], self.CavernZ[1]], [self.CavernY[0], self.archRadius+self.centreOfCurvature["y"]], c="paleturquoise")
        ax.plot([self.CavernZ[0], self.CavernZ[1]], [self.CavernY[0], self.CavernY[0]], c="paleturquoise")
        # Cavern Ceiling
        ax.plot([self.CavernZ[0], self.CavernZ[1]], [self.archRadius+self.centreOfCurvature["y"], self.archRadius+self.centreOfCurvature["y"]], c="paleturquoise")
        # Access Shafts
        ax.plot([self.PX14_Centre["z"]-self.PX14_Radius, self.PX14_Centre["z"]-self.PX14_Radius],
                [self.PX14_Centre["y"], self.PX14_Centre["y"]+self.PX14_Height], c="paleturquoise", label="PX14", alpha=0.5)
        ax.plot([self.PX14_Centre["z"]+self.PX14_Radius,self.PX14_Centre["z"]+self.PX14_Radius],
                [self.PX14_Centre["y"], self.PX14_Centre["y"]+self.PX14_Height], c="paleturquoise", label="PX14", alpha=0.5)
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

        if len(ANUBISrpcs)!=0 and plotRPCs["zy"]: #Include ANUBIS positions
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

        if plotATLAS:
            ax.plot( [self.ATLAS_Z[0],self.ATLAS_Z[0]], [self.ATLAS_Centre["y"]-self.radiusATLAS,self.ATLAS_Centre["y"]+self.radiusATLAS], c="grey")
            ax.plot( [self.ATLAS_Z[0],self.ATLAS_Z[1]], [self.ATLAS_Centre["y"]-self.radiusATLAS,self.ATLAS_Centre["y"]-self.radiusATLAS], c="grey")
            ax.plot( [self.ATLAS_Z[1],self.ATLAS_Z[1]], [self.ATLAS_Centre["y"]-self.radiusATLAS,self.ATLAS_Centre["y"]+self.radiusATLAS], c="grey")
            ax.plot( [self.ATLAS_Z[0],self.ATLAS_Z[1]], [self.ATLAS_Centre["y"]+self.radiusATLAS,self.ATLAS_Centre["y"]+self.radiusATLAS], c="gray")

        plt.xlabel(f"z /m")
        plt.ylabel(f"y /m")
        plt.title("ATLAS Cavern")
        plt.savefig(f"{outDir}/ATLASCavern_ZY{tempSuffix['zy']}.pdf")
        plt.close()

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

        if len(ANUBISrpcs)!=0 and plotRPCs["3D"]: #Include ANUBIS positions
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

        # Plot x,y,z positions of hits that pass and fail geometric acceptance.
        if len(hits)!=0:
            if plotFailed:
                ax.scatter([x[0] for x in hits["failed"]], [z[2] for z in hits["failed"]], [y[1] for y in hits["failed"]], c="red", marker="x")
            ax.scatter([x[0] for x in hits["passed"]], [z[2] for z in hits["passed"]], [y[1] for y in hits["passed"]], c="lime", marker="^")

        if plotATLAS:
            atlasX =  np.linspace(self.ATLAS_Centre["x"]-self.radiusATLAS,self.ATLAS_Centre["x"]+self.radiusATLAS, 100)
            atlasZ =  np.linspace(self.ATLAS_Z[0], self.ATLAS_Z[1], 100)
            
            atlasXX, atlasZZ = np.meshgrid(atlasX, atlasZ)
            atlasYY = np.sqrt(np.power(self.radiusATLAS,2) - np.power(atlasXX-self.ATLAS_Centre["x"],2))
            ax.plot_surface(atlasXX, atlasZZ, atlasYY + self.ATLAS_Centre["y"], rstride=4, cstride=4, alpha=0.4, color="gray")
            ax.plot_surface(atlasXX, atlasZZ, -atlasYY + self.ATLAS_Centre["y"], rstride=4, cstride=4, alpha=0.4, color="gray")

        if len(hits)!=0:
            # Plot a rough impression of the Acceptance
            ax.plot([self.IP["x"], self.CavernX[0]], [self.IP["z"], self.CavernZ[0]], [self.IP["y"],self.obtainCavernYFromX(self.CavernX[0])], 
                     c="k", alpha=0.25, linestyle="--")
            ax.plot([self.IP["x"], self.CavernX[0]], [self.IP["z"], self.CavernZ[1]], [self.IP["y"],self.obtainCavernYFromX(self.CavernX[0])], 
                     c="k", alpha=0.25, linestyle="--")
            ax.plot([self.IP["x"], self.CavernX[1]], [self.IP["z"], self.CavernZ[0]], [self.IP["y"], self.obtainCavernYFromX(self.CavernX[1])], 
                    c="k", alpha=0.25, linestyle="--")
            ax.plot([self.IP["x"], self.CavernX[1]], [self.IP["z"], self.CavernZ[1]], [self.IP["y"],self.obtainCavernYFromX(self.CavernX[1])], 
                    c="k", alpha=0.25, linestyle="--")

        plt.xlabel("x /m")
        plt.ylabel("z /m")
        ax.set_zlabel("y /m")
        ax.set_xlim(-30,30)
        ax.set_ylim(-30,30)
        plt.tight_layout()
        plt.savefig(f"{outDir}/ATLASCavern_XYZ{tempSuffix['3D']}.pdf")
        plt.close()



    def ANUBIS_RPC_positions(self, RPCx=1, RPCy=0.06, RPCz=1.8, overlapAngleXY=0, overlapZ=0, startTheta=-999, layerRadius=0, ID=0):
        #RPC{x,y,z} are the lengths in those directions in the ATLAS Cavern coordinate system used here
        RPCs=[]

        if layerRadius==0:
            layerRadius = self.archRadius - (RPCy/2) # RPC attaches to the ceiling.

        if startTheta==-999:
            # Assume starting on the corner of the ceiling 
            totalTheta = 2*np.arcsin( self.CavernXLength/(2*layerRadius)) # Relative to centre of curvature.
            theta= -totalTheta/2
        else:
            # Assume full ceiling coverage starting from the initial theta value: which is 0 at the vertical x=0
            halfTheta = np.arcsin( self.CavernXLength/(2*layerRadius)) # Relative to centre of curvature.

            if halfTheta < startTheta:
                raise Exception(f"Attempting to start beyond the limit of the ATLAS Cavern at theta of {startTheta}")

            totalTheta = halfTheta - startTheta
            theta = startTheta
            

        # XY positions
        #   - Find the angular size of the RPC to segment the ceiling with
        #   - Each RPC position includes the two corners of the RPC and the midpoint. 
        #   - Effectively unravelling the ceiling into a straight line.
        dTheta = 2*np.arcsin( RPCx/(2*layerRadius))
        Nsegments = np.ceil( totalTheta / (dTheta - overlapAngleXY) ) #Not sure this is correct 

        breakCon = True
        xPos, yPos = [], []
        nIter=0
        while breakCon:
            # Step through theta in units of dTheta - overlap angle
            # Get (x,y) coord of centre of angular segment and save to list
            if nIter == 0:
                overlapTheta = 0
            else:
                overlapTheta = overlapAngleXY 

            tempTheta = theta + (dTheta - overlapTheta)/2 # Midpoint of RPC segment
            if tempTheta > totalTheta/2 or (theta+dTheta) > totalTheta:
                breakCon=False
                break

            bottomOfRPC = layerRadius * np.cos(dTheta/2) - RPCy
            # Convert to cartesian positions with x = r*sin(theta), y = r*cos(theta) relative to Cavern centre
            tempXPos = {"c1": layerRadius*np.sin(theta) + self.centreOfCurvature["x"], 
                        "c2": layerRadius*np.sin(theta + (dTheta - overlapTheta)) + self.centreOfCurvature["x"],
                        "c3": bottomOfRPC*np.sin(theta) + self.centreOfCurvature["x"], 
                        "c4": bottomOfRPC*np.sin(theta + (dTheta - overlapTheta)) + self.centreOfCurvature["x"],
                        "mid": (layerRadius*np.cos(dTheta/2) - (0.5*RPCy))*np.sin(tempTheta) + self.centreOfCurvature["x"],
                        "midTop": layerRadius*np.sin(theta + ((dTheta/2) - overlapTheta)) + self.centreOfCurvature["x"],
            }

            tempYPos = {"c1": layerRadius*np.cos(theta) + self.centreOfCurvature["y"], 
                        "c2": layerRadius*np.cos(theta + (dTheta - overlapTheta)) + self.centreOfCurvature["y"],
                        "c3": bottomOfRPC*np.cos(theta) + self.centreOfCurvature["y"], 
                        "c4": bottomOfRPC*np.cos(theta + (dTheta - overlapTheta)) + self.centreOfCurvature["y"],
                        "mid": (layerRadius*np.cos(dTheta/2) - (0.5*RPCy))*np.cos(tempTheta) + self.centreOfCurvature["y"],
                        "midTop": layerRadius*np.cos(theta + ((dTheta/2) - overlapTheta)) + self.centreOfCurvature["y"],
                        }
            

            xPos.append(tempXPos)
            yPos.append(tempYPos)

            # increment to the next corner of the RPC segments
            theta+=(dTheta - overlapAngleXY)
            
            if theta > totalTheta/2:
                breakCon=False
                break

            nIter+=1

        print(f"There were {nIter} Iterations...")

        #Z Positions:
        breakCon2 = True
        zPos=[]
        Z = self.CavernZ[0]
        nIter2=0
        while breakCon2:
            if nIter2 == 0:
                overlap = 0
            else:
                overlap = overlapZ

            tempZ = Z + ((RPCz-overlap)/2) #Midpoint of RPC
            if tempZ > self.CavernZ[1] or (Z + (RPCz-overlap)) > self.CavernZ[1]:
                breakCon2 = False
                break
            
            zPos.append({"c1": Z, "c2": Z + RPCz - overlap, "mid": tempZ})

            # increment to next RPC segment
            Z+= RPCz - overlap
            if Z > self.CavernZ[1]:
                breakCon2 = False
                break

            nIter2+=1

        print(f"There were {nIter2} Iterations...")


        # Create a list of RPC stations, each with: 
        #   - corners: 8 (x,y,z) coordinates corresponding to their corners, 
        #   - midPoint: The midPoint of the RPC in (x,y,z), 
        #   - "LayerID" and "RPCid": A Layer ID and RPC ID to uniquely identify the RPC
        #   - "plane": A Sympy plane in the eta-phi plane that passes through the midpoint
        posRPC={"x": {"c1": [], "c2": [], "c3": [], "c4": [], "mid": []}, 
                "y": {"c1": [], "c2": [], "c3": [], "c4": [], "mid": []},
                "z": {"c1": [], "c2": [], "mid": []},
        }
        rpcID=0
        for coordZ in zPos:
            for coordX, coordY in zip(xPos, yPos):
                #   Below is a representation of which corner in (x,y,z) each element of the corners list corresponds to.
                #
                #             1--------3
                #            /|       /|
                #           / |      / |
                #          0--------2  |
                #          |  5-----|--7
                # y        | /      | /
                # | z      |/       |/
                # |/       4--------6
                # o--x
                corners=[ (coordX["c1"], coordY["c1"], coordZ["c1"]),
                          (coordX["c1"], coordY["c1"], coordZ["c2"]),
                          (coordX["c2"], coordY["c2"], coordZ["c1"]),
                          (coordX["c2"], coordY["c2"], coordZ["c2"]),
                          (coordX["c3"], coordY["c3"], coordZ["c1"]),
                          (coordX["c3"], coordY["c3"], coordZ["c2"]),
                          (coordX["c4"], coordY["c4"], coordZ["c1"]),
                          (coordX["c4"], coordY["c4"], coordZ["c2"]), ]

                midPoint = (coordX["mid"], coordY["mid"], coordZ["mid"])

                # Define a plane using three points on the top of the RPC segment
                plane = Plane( (coordX["c1"], coordY["c1"],coordZ["mid"]),
                               (coordX["c2"], coordY["c2"],coordZ["mid"]),
                               (coordX["midTop"], coordY["midTop"],coordZ["mid"]))

                RPCs.append( {"corners": corners, "midPoint": midPoint, "plane": plane, "RPCid": rpcID, "LayerID": ID} )
                rpcID+=1

        return RPCs

    def createSimpleRPCs(self, radii, RPCthickness=0.06, origin=[]):
        cav = ATLASCavern()
        if len(origin)==0: # Assume origin of IP unless otherwise stated
            origin = (self.IP["x"], self.IP["y"], self.IP["z"])

        RPCs={"r": [], "theta": [], "phi": []}

        # Radii are assumed to be relative to the centre of curvature
        for r in radii:
            RPCs["r"].append( [r - RPCthickness, r] )
            #RPCs["theta"].append([min(cav.angles["theta"]), max(cav.angles["theta"])])
            RPCs["theta"].append([min(self.angles["theta"]), max(self.angles["theta"])])
            
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


if __name__=="__main__":
    import argparse
    import datetime
    parser = argparse.ArgumentParser(description='Provide an input file')
    parser.add_argument('--remake', action='store_true')
    parser.add_argument('--suffix', type=str, default="")
    parser.add_argument('--hitfile', type=str, help='Path to CSV file with hits')
    # FourVector is x, y, z, t - in mm and s
    args = parser.parse_args()

    print(datetime.datetime.now())

    cav = ATLASCavern()
    cav.createCavernVault()

    print(f"Cavern Bounds: {cav.CavernBounds}")
    print(f"Cavern Corners: {cav.CavernCorners}")
    print(f"IP Location: {cav.IP}")
    print(f"Centre of Curvature: {cav.centreOfCurvature}")
    print(cav.angles)

    # Create ANUBIS RPC Layers
    #   - Afterwards save to a pickle file. 
    #   - Reload RPCs from pickle file if it exists.

    pickleDir = "/usera/dp728/run_dir/pickles"
    if not os.path.exists(pickleDir):
        os.makedirs(pickleDir)

    print("Making RPC Layers...")
    print(datetime.datetime.now())
    if not os.path.exists(f"{pickleDir}/ANUBIS_RPCs_Layer0.pickle") or args.remake:
        #Layer of RPCs (Triplet) attached to the ceiling.
        RPC_Pos1 = cav.ANUBIS_RPC_positions(RPCx=1, RPCy=0.06, RPCz=1.8, overlapAngleXY=0, overlapZ=0, layerRadius = cav.archRadius, ID=0)
        with open(f"{pickleDir}/ANUBIS_RPCs_Layer0.pickle", "wb") as pkl:
            pickle.dump(RPC_Pos1, pkl)
    else:
        with open(f"{pickleDir}/ANUBIS_RPCs_Layer0.pickle", "rb") as pkl:
            RPC_Pos1 = pickle.load(pkl)

    if not os.path.exists(f"{pickleDir}/ANUBIS_RPCs_Layer1.pickle") or args.remake:
        #Singlet RPCs between the Triplets (40cm below Top Triplet).
        RPC_Pos2 = cav.ANUBIS_RPC_positions(RPCx=1, RPCy=0.06, RPCz=1.8, overlapAngleXY=0, overlapZ=0, layerRadius = cav.archRadius-0.40, ID=0)
        with open(f"{pickleDir}/ANUBIS_RPCs_Layer1.pickle", "wb") as pkl:
            pickle.dump(RPC_Pos2, pkl)
    else:
        with open(f"{pickleDir}/ANUBIS_RPCs_Layer1.pickle", "rb") as pkl:
            RPC_Pos2 = pickle.load(pkl)

    if not os.path.exists(f"{pickleDir}/ANUBIS_RPCs_Layer2.pickle") or args.remake:
        #Bottom Triplet RPCs (1m below Top Triplet).
        RPC_Pos3 = cav.ANUBIS_RPC_positions(RPCx=1, RPCy=0.06, RPCz=1.8, overlapAngleXY=0, overlapZ=0, layerRadius = cav.archRadius-1, ID=0)
        with open(f"{pickleDir}/ANUBIS_RPCs_Layer2.pickle", "wb") as pkl:
            pickle.dump(RPC_Pos3, pkl)
    else:
        with open(f"{pickleDir}/ANUBIS_RPCs_Layer2.pickle", "rb") as pkl:
            RPC_Pos3 = pickle.load(pkl)

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
    ANUBISstations = RPC_Pos1
    ANUBISstations.extend(RPC_Pos2)
    ANUBISstations.extend(RPC_Pos3)
    ANUBISstations = cav.convertRPCList(ANUBISstations)
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

    ANUBISstations = cav.createSimpleRPCs([cav.archRadius, cav.archRadius-0.4, cav.archRadius-1], RPCthickness=0.06)
    print(ANUBISstations)
    minStationRadius = min(min(ANUBISstations["r"]))
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
    hits_df = pd.read_csv(args.hitfile)
    import ast

    hits_df['decayVertexParsed'] = hits_df['decayVertex'].apply(ast.literal_eval)
    passedHits, failedHits = [], []
    passedevents, failedevents = [], []
    eventeta = []
    eventphi = []
    allhits = []
    for i, row in hits_df.iterrows():
        X, Y, Z = row['decayVertexParsed'][0], row['decayVertexParsed'][1], row['decayVertexParsed'][2]
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
        r, theta, phi = row['decayVertexDist'], h.to_theta(row['eta']),row['phi']


        inCavern = cav.inCavern(X, Y, Z, maxRadius=minStationRadius - 0.20)
        inATLAS = cav.inATLAS(X, Y, Z, trackingOnly=True)
        intANUBIS = cav.intersectANUBISstationsSimple2(X, Y, Z,r, theta, phi, ANUBISstations)

        
        if ( inCavern and not inATLAS and (len(intANUBIS[1]) >= 2) ):
            passedHits.append((X,Y,Z))
            passedevents.append(row["eventNumber"])


        else:
            failedHits.append((X,Y,Z))
            failedevents.append(row["eventNumber"])
        allhits.append((X,Y,Z))

        eventeta.append(row["eta"])
        eventphi.append(row["phi"])

        if i % 1000 == 0:
            print(f"Checked {i}/{len(hits_df)} hits", end="\r", flush=True)

    hitEnd2=datetime.datetime.now()
    print(hitEnd2)
    print(f"Took {hitEnd2 - hitEnd1}")
    print(f"Passed: {len(passedHits)} ({len(passedHits)/(len(passedHits)+len(failedHits))}%)  | Failed: {len(failedHits)} ({len(failedHits)/(len(passedHits)+len(failedHits))}%)")
    
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
    min_x, max_x = min(x_vals), max(x_vals)
    min_y, max_y = min(y_vals), max(y_vals)
    min_z, max_z = min(z_vals), max(z_vals)
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
    cav.plotFullCavern(ANUBISrpcs=[RPC_Pos1, RPC_Pos2, RPC_Pos3], plotATLAS=True, plotFailed=False, 
                        hits={"passed": allhits, "failed": failedHits, "bins": hitBins}, suffix=f"precut") #{args.suffix}


    cav.plotFullCavern(ANUBISrpcs=[RPC_Pos1, RPC_Pos2, RPC_Pos3], plotATLAS=True, plotFailed=False, 
                                   hits={"passed": passedHits, "failed": failedHits, "bins": hitBins}, suffix=f"_WithHits{args.suffix}")
    print(datetime.datetime.now())

    plt.figure(figsize=(10,6))
    cuteventeta = [e for e in eventeta if -15 <= e <= 15]
    plt.hist(cuteventeta, bins=100, color='purple', alpha=0.75)
    plt.xlabel("η (pseudorapidity)")
    plt.ylabel("Number of hits")
    plt.title("Pseudorapidity Distribution")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("pseudorapidity_distribution.pdf")
    #plt.show()

    plt.figure(figsize=(10,6))
    plt.hist(eventphi, bins=90, color='orange', alpha=0.75)
    plt.xlabel("ϕ (phi) [radians]")
    plt.ylabel("Number of hits")
    plt.title("Azimuthal Angle (ϕ) Distribution")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("phi_distribution.pdf")
    #plt.show()
