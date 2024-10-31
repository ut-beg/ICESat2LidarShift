#=========================================================================================================
# LidarShift.py
# Aaron Averett and Kutalmis Saylam - Bureau of Economic Geology - The University of Texas at Austin
# Copyright 2023, The University of Texas at Austin
#
# This is a Python script for adjusting the data produced by BEG's Chiroptera 4X lidar sensor during the time period between the satellite and corresponding aircraft overflight during the July of 2022 field campaign.
#=========================================================================================================

from cmath import isnan
from ctypes import ArgumentError
from re import I
from statistics import mean, stdev
import pandas as pd
import datetime
import arcpy
import os
import argparse
import io
from csv import writer as csvWriter
import math
import numpy as np
#Need this extra import to provide a local reference to numpy.power when we call it from a dataframe.query
#Apparently, you can't do np.power() in a df.query
from numpy import median, power as nppower
import datetime
from numba import jit
from scipy import stats
import multiprocessing
import traceback

#ATL03 Spot size (meters)
ATL03_SPOT_SIZE = 6

ATL03_SPOT_SIZE_SELECTION_STRING = "{0} Meters".format(ATL03_SPOT_SIZE)

#The character to use to split each line
SPLIT_CHAR=','

#Enable automatic output overwriting
arcpy.env.overwriteOutput = True

@jit("(float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64[:], float64, float64, float64, float64, float64, float64, float64, float64, float64, float64, bool_)", nopython=True)
def shiftPoints(npaTime, npaX, npaY, npaZ, npaDeltaT, npaEasting, npaNorthing, driftXVel, driftYVel, extraDriftX, extraDriftY, chiropteraTotalTime, chiropteraT0, chiropteraT1, atl07T0, atl07T1, atl07TotalTime, sameDirection):
    """
    Shifts the points in the given numpy arrays based on the time, speed and direction paraemters.
    """
    
    numPoints = len(npaTime)

    for i in range(numPoints):
        t = npaTime[i]
        x = npaX[i]
        y = npaY[i]
        z = npaZ[i]
        
        
        #The time, as a fraction of the total duration of the Chiroptera overflight, at which this point was measured.
        chiropteraTimeFraction = (t - chiropteraT0) / chiropteraTotalTime

        atl07TimeFraction = 0

        #If we're going the same direction, we use the same time fraction for IceSAT2.  If we're going opposite directions, we do 1 - time fraction.
        if sameDirection:
            atl07TimeFraction = chiropteraTimeFraction
        else:
            atl07TimeFraction = 1 - chiropteraTimeFraction

        atl07TimeOnTarget = atl07T0 + (atl07TotalTime * atl07TimeFraction)

        #Delta T is the time difference between Chiroptera actually hitting this spot and the ATL07 "hypothetically" hitting it.
        deltaT = t - atl07TimeOnTarget

        #Compute X and Y components of our ice drift vector
        chiropteraDeltaX = deltaT * driftXVel + extraDriftX
        chiropteraDeltaY = deltaT * driftYVel + extraDriftY
         
        #Compute the final X and Y position of our point, adjusted for reversal of the sea ice drift.           
        xFinal = x - chiropteraDeltaX
        yFinal = y - chiropteraDeltaY
        zFinal = z

        #Compute a time that represents when IceSAT2 would have hit this point, if it had a Chiroptera installed on it.
        tFinal = t + deltaT
        
        npaDeltaT[i] = deltaT
        npaEasting[i] = xFinal
        npaNorthing[i] = yFinal
        

def shiftFeatureClass(featureClassPath, outFeatureClassPath, driftMetersPerSecond, driftBearing, extraDrift, chiropteraT0, chiropteraT1, atl07T0, atl07T1, atl07TotalTime, sameDirection, outSR):
    """
    Performs the shift operation on the given feature class and writes the output to the outFeatureClass variable
    """
    
    print("Working on shift.  Drift speed is {0} m/s.  Drift bearing is {1}.".format(driftMetersPerSecond, driftBearing))
    #Convert our drift bearing to radians, relative to due east.
    driftPolarDir = 360 + (90.0 - driftBearing)

    driftPolarDirRad = math.radians(driftPolarDir)

    #Calculate the X and Y components of our drift velocity vector
    driftXVel = math.cos(driftPolarDirRad) * driftMetersPerSecond
    driftYVel = math.sin(driftPolarDirRad) * driftMetersPerSecond

    #Calculate the X and Y components of our extra 
    extraDriftX = math.cos(driftPolarDirRad) * extraDrift
    extraDriftY = math.sin(driftPolarDirRad) * extraDrift
    
    columns = ["Shape@X", "Shape@Y", "c", "GPS_time", "chiroptera_z"] #List the fields you want to include. I want all columns except the geometry
    df = pd.DataFrame(data=arcpy.da.SearchCursor(featureClassPath, columns), columns=columns)
    
    npaTime = df["GPS_time"].to_numpy()
    npaX = df["Shape@X"].to_numpy()
    npaY = df["Shape@Y"].to_numpy()
    npaZ = df["chiroptera_z"].to_numpy()
    npaC = df["c"].to_numpy()
    
    arrayLen = len(npaTime)
    
    npaDeltaT = np.zeros(arrayLen)
    npaEasting = np.zeros(arrayLen)
    npaNorthing = np.zeros(arrayLen)
    
    chiropteraTotalTime = chiropteraT1 - chiropteraT0    

    timeBefore = datetime.datetime.now()

    shiftPoints(npaTime, npaX, npaY, npaZ, npaDeltaT, npaEasting, npaNorthing, driftXVel, driftYVel, extraDriftX, extraDriftY, chiropteraTotalTime, chiropteraT0, chiropteraT1, atl07T0, atl07T1, atl07TotalTime, sameDirection)
    
    timeAfter = datetime.datetime.now()
    
    timeDiff = timeAfter - timeBefore
    
    print(timeDiff)

    ret = pd.DataFrame()
    
    ret["c"] = npaC
    ret["t"] = npaTime
    ret["x"] = npaEasting
    ret["y"] = npaNorthing
    ret["z"] = npaZ
    """
    ret["dt"] = npaDeltaT
    ret["x_old"] = npaX
    ret["y_old"] = npaY
    """ 
    
    return ret
    
def calcChiropteraElevationForATL03Point(chiropteraDataFrame, atl03X, atl03Y, spotSize):
    """
    Determines the elevation to be used as the corresponding "chiroptera elevation" for the given ATLAS sensor point.
    This works by first querying the dataframe for any points falling within the spot size radius of the given X and Y, which theoretically represent the center of the area illuminated by the ATLAS sensor's laser.
    """
    
    shape1 = chiropteraDataFrame.shape

    relevantPoints = chiropteraDataFrame.query("sqrt(@nppower((`x` - @atl03X), 2) + @nppower((`y` - @atl03Y), 2)) <= @spotSize")

    shape2 = relevantPoints.shape

    maxZ = relevantPoints["z"].mean()

    return maxZ, shape2[0]

def compareToATL03(chiropteraDataFrame, atl03DataFrame, cMaxX, cMinX, cMaxY, cMinY):
    """
    Performs the comparison of the given set of chiroptera data to the given set of ATLAS sensor data.
    """

    print("Comparing with ATL03 data...")

    #Array to hold our difference values
    chiropteraElevs = []
    atl03Elevs = []
    pointCounts=[]
   
    atl03InBoundsDataFrame = atl03DataFrame.query("`SHAPE@X` >= @cMinX and `SHAPE@X` <= @cMaxX and `SHAPE@Y` >= @cMinY and `SHAPE@Y` <= @cMaxY")

    for ind, row in atl03InBoundsDataFrame.iterrows():

        x = row["SHAPE@X"]
        y = row["SHAPE@Y"]
        z = row["z"]
       
        chiropteraElev, pointCount = calcChiropteraElevationForATL03Point(chiropteraDataFrame, x, y, ATL03_SPOT_SIZE)

        if not isnan(chiropteraElev) and not isnan(z) and not isnan(pointCount):

            chiropteraElevs.append(chiropteraElev)
            atl03Elevs.append(z)
            pointCounts.append(pointCount)
            

    return chiropteraElevs, atl03Elevs, pointCounts

def processWithSpeedAndBearing(argSet):
    """
    Performs the iterative process of first shifting the given chiroptera data by the given speed and distance parameters, and then performs the comparison and linear regression operations.
    This function takes its parameters as a dictionary because it is meant to be launched via multiprocessing library.
    """
    
    ret = None

    try:
        
        curSpeed = argSet["curSpeed"]
        curBearing = argSet["curBearing"]
        chiropteraDfArray = argSet["chiropteraDfArray"]
        atl03DfArray = argSet["atl03DfArray"]
        outputDataPath = argSet["outputDataPath"]
        chiropteraT0 = argSet["chiropteraT0"]
        chiropteraT1 = argSet["chiropteraT1"]
        atl07T0 = argSet["atl07T0"]
        atl07T1 = argSet["atl07T1"]
        sameDirection = argSet["sameDirection"]
        exportShiftedData = argSet["exportShiftedData"]
        
        #If we don't have a workspace path set, set it to our output path.
        if arcpy.env.workspace == None:
            arcpy.env.workspace = outputDataPath

        sr_utm = arcpy.SpatialReference(32616)

        atl07TotalTime = atl07T1 - atl07T0
    
        chiropteraUtmFCPath = os.path.join("memory", "chiropterautmdata")
        arcpy.da.NumPyArrayToFeatureClass(chiropteraDfArray, chiropteraUtmFCPath, ("Easting", "Northing"), sr_utm)
    
        atl03FCPath = os.path.join("memory", "atl03utmdata")
        arcpy.da.NumPyArrayToFeatureClass(atl03DfArray, atl03FCPath, ("x", "y"), sr_utm)

        print("Beginning loop for bearing: {0} and speed: {1}".format(curBearing, curSpeed))

        loopStartTime = datetime.datetime.now()
        
        driftMetersPerSecond = curSpeed / 60.0

        fcName = "shifted_{0}_{1}".format(str(curBearing), str(curSpeed))
        fcName = fcName.replace(".","_")
    
        shiftedFcPath = os.path.join("memory", fcName)
    
        #Perform the time/speed/distance derived geometric shift of the chiroptera data, using the given parameters.
        dfShifted = shiftFeatureClass(chiropteraUtmFCPath, shiftedFcPath, driftMetersPerSecond, curBearing, 0, chiropteraT0, chiropteraT1, atl07T0, atl07T1, atl07TotalTime, sameDirection, sr_utm)
        
        chiropteraDataFrame = dfShifted
        
        cMaxX = chiropteraDataFrame["x"].max()
        cMinX = chiropteraDataFrame["x"].min()
        cMaxY = chiropteraDataFrame["y"].max()
        cMinY = chiropteraDataFrame["y"].min()

        #Load up the ATL03 data
        npaATL03 = arcpy.da.FeatureClassToNumPyArray(atl03FCPath, ["SHAPE@X", "SHAPE@Y", "z"])
        atl03DataFrame = pd.DataFrame(npaATL03)
    
        chiropteraZs, atl03Zs, pointCounts = compareToATL03(chiropteraDataFrame, atl03DataFrame, cMaxX, cMinX, cMaxY, cMinY)

        dfZs = pd.DataFrame(columns=["c","a"])

        dfZs["c"] = chiropteraZs
        dfZs["a"] = atl03Zs

        zsFullPath = os.path.join(outputDataPath, "{0}_zs.csv".format(fcName))
        dfZs.to_csv(zsFullPath)
        
        if exportShiftedData:
            exportShiftedDataPath = os.path.join(outputDataPath, "{0}.csv".format(fcName))
            chiropteraDataFrame.to_csv(exportShiftedDataPath, index=False, header=False)
    
        resultRow = [curBearing, driftMetersPerSecond * 60.0, None, None, None, None, None, None, None, 0, None, None, None, None, None]
    
        if len(chiropteraZs) > 0 and len(atl03Zs) > 0 and len(chiropteraZs) == len(atl03Zs):
            
            chiropteracount_min = min(pointCounts)
            chiropteracount_max = max(pointCounts)
            chiropteracount_mean = mean(pointCounts)
            chiropteracount_median = median(pointCounts)
            chiropteracount_stdev = stdev(pointCounts)
    
            slope, intercept, r_value, p_value, std_err = stats.linregress(atl03Zs, chiropteraZs)

            rsquared = r_value ** 2
  
            npacz = np.array(chiropteraZs)
            npaatlz = np.array(atl03Zs)

            rmse = math.sqrt(np.square((np.subtract(npaatlz, npacz))).mean())
            print("R^2: {0}, RMSE: {1}".format(rsquared, rmse))
        
            ret = [curBearing, driftMetersPerSecond * 60.0, slope, intercept, r_value, p_value, std_err, rsquared, rmse, len(chiropteraZs), chiropteracount_min, chiropteracount_max, chiropteracount_mean, chiropteracount_median, chiropteracount_stdev]

    
    
        loopEndTime = datetime.datetime.now()
        print("Loop run time was: {0}".format((loopEndTime - loopStartTime)))
    
        if arcpy.Exists(chiropteraUtmFCPath):
            arcpy.management.Delete(chiropteraUtmFCPath)
        
        if arcpy.Exists(atl03FCPath):
            arcpy.management.Delete(atl03FCPath)
            
    except Exception as e:
        print(e)
        
        st = traceback.format_exc()
        
        print(st)
        exit()
    
    return ret
    


def main():
    #Parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="Path to chiroptera data file")
    parser.add_argument("-sif", type=str, required=True, help="Path to ATL03 data file")
    parser.add_argument("-w", type=str, required=True, help="Workspace")
    parser.add_argument("-si", type=float, required=True, help="Lowest drift speed in meters per *minute*")
    parser.add_argument("-sf", type=float, required=True, help="Highest drift speed in meters per minute")
    parser.add_argument("-ss", type=float, required=True, help="Drift speed step in meters per minute")
    parser.add_argument("-bi", type=float, required=True, help="Initial drift bearing in degrees")
    parser.add_argument("-bf", type=float, required=True, help="Final drift bearing in degrees")
    parser.add_argument("-bs", type=int, required=True, help="Drift bearing step")
    parser.add_argument("-e", type=float, default=0, required=False, help="Optional extra fixed drift along bearing.  Default 0.")
    parser.add_argument("-t0", type=float, required=True, help="Section start time in GPS week seconds for ATL07")
    parser.add_argument("-t1", type=float, required=True, help="Section end time in GPS week seconds for ATL07")
    parser.add_argument("-ct0", type=float, required=True, help="Section start time in GPS week seconds Chiroptera")
    parser.add_argument("-ct1", type=float, required=True, help="Section end time in GPS week seconds Chiroptera")
    parser.add_argument("-o", type=str, required=True, help="Output file path")
    parser.add_argument("-m", type=str, required=True, help="Number of parallel instances to run")
    parser.add_argument("-sd", type=str, required=True, help="Whether satellite and aircraft pass are in the same direction.  y for same direction, n for opposing directions")
    parser.add_argument("-ed", type=str, required=True, help="Whether to export the actual shifted points.  Y for yes, n for no.")

    args = parser.parse_args()

    iSensorType = 0

    #Input file path
    inputDataPath = args.i
    workspaceDataPath = args.w
    outputDataPath = args.o
    inputSatelliteDataPath = args.sif

    #Drift speed in meters per minute - this is converted to m/s later
    initialDriftMetersPerMinute = args.si
    finalDriftMetersPerMinute = args.sf
    driftSpeedStep = args.ss

    initialDriftBearing = args.bi
    finalDriftBearing = args.bf
    driftBearingStepSize = args.bs

    if initialDriftBearing < 0 or initialDriftBearing > 360:
        raise ArgumentError("Initial drift bearing is out of bounds.  Must be between 0 and 360.")
    
    if finalDriftBearing < 0 or finalDriftBearing > 360:
        raise ArgumentError("Final drift bearing is out of bounds.  Must be between 0 and 360.")

    if driftBearingStepSize < 0.1 or driftBearingStepSize > 360:
        raise ArgumentError("Drift bearing steps must be at least 0.1 and less than 360")
    
    if initialDriftBearing > finalDriftBearing:
        raise ArgumentError("Initial drift bearing must not be less than final drift bearing.")
    
    if driftSpeedStep <= 0:
        raise ArgumentError("Drift speed step must be greater than zero.") 
    
    #An extra arbitrary adjustment in the drift bearing direction (in m), to allow adjustment for the limitations in accuracy of this approach.
    #We can get in the ballpark by mathematically reversing the drift, but 
    extraDrift = args.e

    #Whether or not we need to export the full shifted data.
    exportShiftedData = False
    if args.ed == "y":
        exportShiftedData = True

    atl07T0 = -1
    atl07T1 = -1

    #This indicates whether the overflights were conducted in the same or opposing direction
    sameDirection = False
    
    if args.sd == "y":
        sameDirection = True

    #ATL07 pass parameters for line 1
    sr_utm = arcpy.SpatialReference(32616)

    #Beginning and end time of IceSAT2 overflight of this area (in GPS seconds)
    atl07T0 = args.t0
    atl07T1 = args.t1

    chiropteraT0 = args.ct0
    chiropteraT1 = args.ct1
    chiropteraTotalTime = chiropteraT1 - chiropteraT0

    chiropteraFile = inputDataPath
    atl03File = inputSatelliteDataPath

    #Initialize the path to our geodatabase, which is required for some of the arcpy operations.
    gdbPath = workspaceDataPath

    #Set our default workspace path to be the same as our geodatabase.
    arcpy.env.workspace = gdbPath

    #Create a spatial reference object for the UTM zone where our data were collected.
    sr_utm = arcpy.SpatialReference(32616)

    chiropteraFCPath = os.path.join("memory", "chiropteraData")
    chiropteraUtmFCPath = os.path.join(outputDataPath, "chiropteraUtmData")

    arcpy.env.overwriteOutput = True
    

    #Compute total duration of IceSAT2 pass
    atl07TotalTime = atl07T1 - atl07T0

    startTime = datetime.datetime.now()

    print("Starting up at {0}".format(startTime))

    print("Ingesting ATL03 data")

    atl03Df = pd.read_fwf(atl03File,header=None)
    atl03FCPath = os.path.join("memory", "atl03Data")
    atl03DfWithNames= atl03Df.rename(columns={0: "x", 1:"y", 2:"z"})
    atl03DfArray = atl03DfWithNames.to_records(index=False)
    arcpy.da.NumPyArrayToFeatureClass(atl03DfArray, atl03FCPath, ("x", "y"), sr_utm)

    #Load up the original, un-shifted chiroptera data
    rawChiropteraDataDf = pd.read_fwf(inputDataPath, header=None)
    rawChiropteraDataWithNamesDf = rawChiropteraDataDf.rename(columns={0: "c", 1: "GPS_time", 2: "Easting", 3: "Northing", 4: "chiroptera_z"})
    chiropteraDfArray = rawChiropteraDataWithNamesDf.to_records(index=False)

    #Load the ATL03 data into a dataframe, which we'll use later.
    npaATL03 = arcpy.da.FeatureClassToNumPyArray(atl03FCPath, ["SHAPE@X", "SHAPE@Y", "z"])
    atl03DataFrame = pd.DataFrame(npaATL03)

    #Get ready to do the iterative process...
    curBearing = initialDriftBearing
    curSpeed = initialDriftMetersPerMinute

    bearingStep = driftBearingStepSize

    dfParamSets = pd.DataFrame(columns=["driftbearing", "driftspeed", "slope","intercept", "r","p", "std_err","rsquared", "rmse", "intersectCount", "chiropteracount_min", "chiropteracount_max", "chiropteracount_mean", "chiropteracount_median", "chiropteracount_stdev"])

    runArgSets = []

    #Iterate over bearing steps
    while curBearing <= finalDriftBearing:
        #Iterate over speed steps.
        while curSpeed <= finalDriftMetersPerMinute:
        
            #Compose a set of arguments to be used for this iteration.  Note that we'll actually launch this later.
            runArgSet = {
                "curSpeed": curSpeed,
                "curBearing": curBearing,
                "chiropteraDfArray": chiropteraDfArray,
                "outputDataPath": outputDataPath,
                "atl03DfArray": atl03DfArray,
                "chiropteraT0": chiropteraT0,
                "chiropteraT1": chiropteraT1,
                "atl07T0": atl07T0,
                "atl07T1": atl07T1,
                "sameDirection": sameDirection,
                "exportShiftedData": exportShiftedData
            }
        
            runArgSets.append(runArgSet)

            #Increment our speed.
            curSpeed = curSpeed + driftSpeedStep
    
        #Increment our bearings
        curBearing = curBearing + bearingStep
        curSpeed = initialDriftMetersPerMinute

    #Do the heavy lifting.  Two methods are provided.  One uses multi-processing, and the other does everything in sequence.  The multi processing version is obviously much faster, but may be problematic on some systems.
    #It's also tough to debug.

    if int(args.m) > 1:
        #Multi-processing code path
        finalMulti = int(args.m)
        p = multiprocessing.Pool(finalMulti)
        resultRows = p.map(processWithSpeedAndBearing, runArgSets)
    
        for resultRow in resultRows:
            dfParamSets.loc[-1] = resultRow
            dfParamSets.index = dfParamSets.index + 1
            dfParamSets.sort_index()

    else:
        #Sequential version.  Slower, but reliable.
        for argSet in runArgSets:
            resultRow = processWithSpeedAndBearing(argSet)
        
            dfParamSets.loc[-1] = resultRow
            dfParamSets.index = dfParamSets.index + 1
            dfParamSets.sort_index()
        
    #Write out the results report file.
    paramSetsCsvPath = os.path.join(outputDataPath, "paramresults.csv")

    dfParamSets.to_csv(paramSetsCsvPath)

    finishTime = datetime.datetime.now()
    print("Finished at {0}".format(finishTime))
    
    timeDiff = finishTime - startTime

    print("running time is: {0}".format(timeDiff))



#Required to keep the multi-processing library from trying to launch the whole shebang again when the processor threads are called.
if __name__ == '__main__':
	#Run the program
	main()