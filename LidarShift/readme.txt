#=========================================================================================================
# Readme for LidarShift.py
# Aaron Averett and Kutalmis Saylam - Bureau of Economic Geology - The University of Texas at Austin
# Copyright 2023, The University of Texas at Austin
#=========================================================================================================

This directory contains the source code for LidarShift.py, a Python script for adjusting the data produced by BEG's Chiroptera 4X lidar sensor during the time period between the satellite and corresponding aircraft overflight during the July of 2022 field campaign.

The script works by calculating the time difference between each point in the supplied Chiroptera data file and the corresponding time at which the ICESat-2 satellite passed over that location.  It then calculates a vector based on the input speed and bearing for the X and Y components of that point that corrects its position to what it would have been at the time of the satellite overflight.

The X and Y components of the correction vector are calculated using the following equations, with intputs being the time difference (t), drift velocity in m/s (V) and drift bearing (θ).

X_C = X_0 + (t * cos⁡(θ) * V)

Y_C = Y_0 + (t * sin⁡(θ) * V)

This process is done iteratively, once for each combination of drift speed and bearing that falls within the given ranges.  The results of each shift operation are then compared to the values found in the supplied ATLAS sensor data file, and a linear regression is performed.  The results of that comparison are finally written to a report file in .csv format.

A complete set of example data files and a batch file containing the appropriate commands to launch the script are included.

Requirements:

This script requires the Python environment installed with ESRI ArcGIS Pro 2.9 or later, along with the "numba" Python package.  For information on installation of additional packages into ArcGIS Pro's Python environment, see here:
https://pro.arcgis.com/en/pro-app/3.0/arcpy/get-started/add-a-package.htm

Installation:

Place LidarShift.py in a directory you have access to, and create at least one directory on a writeable drive to use as the output directory, given as -o at the command line.  You must also ensure that you have an ArcGIS Python environment configured with the numba package.

Finally, use ArcGIS Pro to create an empty file geodatabase, and copy the path, to be used as the -w argument.

Usage (all parameters are required):

propy.bat LidarShift.py 
	-i [Input file - this is a fixed-width text file containing the Chiroptera point data]
	-sif [Satellite input file - this is a fixed-width text file containing the satellite point data]
	-w [ArcGIS Workspace Path - the path to an ArcGIS File Geodatabase] 
	-si [Initial drift speed]
	-sf [Final drift speed]
	-ss [Drift speed step in meters per minute]
	-bi [Initial compass bearing]
	-bf [Final compass bearing.  Note - you may not get an iteration that uses this exactly if the difference between -bf and -bi is not evenly divisble by -bs]
	-bs [Compass bearing step size, in degrees, for each iteration]
	-t0 [GPS time for beginning of satellite pass]
	-t1 [GPS time for end of satellite pass]
	-ct0 [GPS time for beginning of Chiroptera pass]
	-ct1 [GPS time for end of Chiroptera pass]
	-o [Path to output directory]
	-m [Number of simultaneous iterations to be run, using Python multiprocessing library.  Should be less than number of hardware CPU cores.  Enter 1 to run without multiprocessing.] 
	-sd [Whether the satellite and aircraft passes are flown in the same direction.  y for same direction, n for opposing directions]
	-ed [y/n, whether to save the shifted data, or discard it after the comparison.  y for yes, n for no]

Chiroptera Fixed Width Format:

This is fixed-width text file with no column headings.  Each row represents a single lidar return, with five paramaters, one for each column.  The column widths do not need to be any specific size, but should be consistent from row to row throughout the file.

The five parameters are:

Class: LAS file point class
Time: GPS time
X: UTM easting
Y: UTM northing
Z: Height measurement

Example of Chiroptera fixed width format:

9.000000 233764.125109 577869.070000 9518705.140000 0.022123
9.000000 233764.125115 577869.320000 9518704.660000 -0.007882
9.000000 233764.125122 577869.570000 9518704.169900 -0.037887
9.000000 233764.125129 577869.800000 9518703.680000 -0.027891
9.000000 233764.125135 577870.030000 9518703.190000 0.002104

Satellite Fixed Width Format:

This is a fixed-width text file with no column headings.  Each row represents a single point, with four parameters, one for each column.  

The three parameters are:

X: UTM easting
Y: UTM northing
Z: Height measurement

Example of Satellite Fixed Width Format:

578527.39 9520949.78 0.30
578526.28 9520946.33 0.23
578524.92 9520942.15 0.35
578523.35 9520937.31 0.09
578522.24 9520933.89 0.00
