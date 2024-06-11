import arcpy
import argparse
import pandas as pd
import math
from pathlib import Path
import os
import numpy as np

DEFAULT_POINT_CLASS=9

def main():
    #Parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="Path to chiroptera data file")

    args = parser.parse_args()
    
    inputDataPath = args.i

    rawChiropteraDataDf = pd.read_csv(inputDataPath)
    
    chiropteraWgsFCPath = "memory\\chiropteradataWgs"

    sr_wgs = arcpy.SpatialReference(4326)
    sr_utm = arcpy.SpatialReference(32616)
    
    npaRawChiropteraData = rawChiropteraDataDf.to_records(index=False)

    arcpy.da.NumPyArrayToFeatureClass(npaRawChiropteraData, chiropteraWgsFCPath, ("longitude", "latitude"), sr_wgs)
        
    arcpy.management.AddField(chiropteraWgsFCPath, "z_tide_free", "FLOAT")
    arcpy.management.AddField(chiropteraWgsFCPath, "c", "LONG")
    
    with arcpy.da.UpdateCursor(chiropteraWgsFCPath, ["SHAPE@X", "SHAPE@Y", "corr_height_rel006", "z_tide_free", "c"]) as ucursor:
        for row in ucursor:
            lat = row[1]
            lng = row[0]
            zMeanTide = row[2]
                        
            latRadians = math.radians(lat)
            geoid_free2Mean = (0.1287 - (0.3848 * math.pow(math.sin(latRadians), 2)))
            zMssTideFree = zMeanTide + geoid_free2Mean
            zFinal = zMssTideFree
            
            row[3] = zFinal
            row[4] = DEFAULT_POINT_CLASS
            
            ucursor.updateRow(row)
            
    chiropteraUtmFCPath = "memory\\chiropteradataUtm" 

    arcpy.management.Project(chiropteraWgsFCPath, chiropteraUtmFCPath, sr_utm)
    
    npaChiropteraUtmData = arcpy.da.FeatureClassToNumPyArray(chiropteraUtmFCPath, ["c", "GPS_time",  "SHAPE@X", "SHAPE@Y", "z_tide_Free"])
    
    filenameStem = Path(inputDataPath).stem
    
    finalFileName = filenameStem + "_ellipsoid_utm.txt"
    
    dirname = os.path.dirname(inputDataPath)
    
    finalOutputPath = os.path.join(dirname, finalFileName)
    
    np.savetxt(finalOutputPath, npaChiropteraUtmData, fmt='%f')

#Required to keep the multi-processing library from trying to launch the whole shebang again when the processor threads are called.
if __name__ == '__main__':
	#Run the program
	main()