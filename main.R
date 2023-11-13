
# Steepest paths tool
# v0.10 (2023-11-11)
# main R script executing preparatory operations


################################ Description ###################################
###
# This program implements a simplified model of forest fire propagation in
# mountainous regions by calculating the steepest paths, i.e. the 3D polylines
# that ascent the slope following the direction of the maximum gradient and
# moving on land covered by woody vegetation (trees or shrubs), solely on the
# basis of a digital elevation model depicting the terrain topography and a 
# vector map of the forested areas. On mountainsides and in the absence of
# strong winds, the steepest paths can be considered as a first approximation of
# the likeliest fire spread paths, that is, the trajectories on the ground
# surface with the highest probability of forest fire propagation.
# 
#
# 1. General procedure
# ````````````````````
# Fires are modelled in such a way that they proceed uphill in direction of the
# maximal slope (i.e. the gradient). In order to more accurately model the
# behaviour of forest fires in the wild a separate procedure is adopted upon
# reaching a local maximum. In this case, the ascent stops and the surroundings
# (see rmax below) are searched for higher points to which the fire might jump.
# This is repeated until there is no higher point available. Three kinds of
# lines are returned by the algorithm:
# tot: The total trajectory of the fire ignoring the fact that it may exit
#      from the forested region.
# blocked: The trajectory of the fire up until the first non-forested region
#          of length at least maxnoforest.
# last_forest: The trajectory of the fire up until the last point at which
#              the fire leaves the forest.
#
# 2. Technical details
# ````````````````````
# The above algorithm is implemented by using the explicit Euler method (speed
# is much more relevant than precision) to solve the differential equation
# dx/dt=\nabla*h where h is the height at position x. The ode solver is
# interrupted when the speed of the fire falls below vmin. At this point, the
# search procedure outlined in the previous section is performed. A circle of
# radius rmax is drawn around the current position and only points lying within
# this circle are considered. Each jump must increase the height of the fire by
# more than mindeltah meters, otherwise the calculation of the trajectory is
# terminated. This is to prevent erratic behaviour in flat topography which is
# unrealistic. Note that the functionality provided in `calc_lines.r` could also
# be provided in the form of a function. However, we strongly advise against such
# an approach. Namely, should the procedure fail at any point within the function
# (e.g. if the Spatialite module was not installed correctly), then the local
# context of the function will be destroyed meaning that all progress up until
# that point is lost. We therefore prefer the "source()" solution, which also gives
# more descriptive error messages containing specific line numbers.
#
# 2.1 Input Data
# ^^^^^^^^^^^^^^
# Two kinds of input data are required to run the simulation. Firstly, a
# shapefile containing the polygons representing zones covered by forest must be
# supplied. This may be done via the "forest_file" parameter below. Secondly, a
# digital elevation model (DEM) must be supplied in .tif format containing
# information about the elevation of the terrain within the region of interest.
# This DEM should be a plane grid with square pixels and use a projected
# coordinate system (CRS) which form a Cartesian reference system. Please make
# sure that the shapefile and the DEM have the same well-defined CRS.
# The script is set up to work with coordinate systems that use the meter as
# reference unit of measure. If you select a CRS with another unit you probably
# have to readjust the script in `calc_lines.r`.
#
# 2.2 Output Data
# ^^^^^^^^^^^^^^^
# A spatialite database is created containing the output of the algorithm. Documentation
# for this database can be found in `README.md`. Note that other tables and columns
# *may* be present. However, these are relevant only for calculation and should
# be disregarded if you do not know what they are.



############################# Required packages ################################
###
# Install required packages
# install.packages("RSQLite")
# install.packages("raster")
# install.packages("sf")
# install.packages("terra")
# install.packages("DBI")
# install.packages("tidyr")
# install.packages("abind")
# install.packages("prodlim")
# install.packages("deSolve")
# install.packages("dplyr")
# install.packages("ctmcmove")
# install.packages("lwgeom")
# install.packages("devtools")
# library(devtools)
# devtools::install_version("rgeos", version="0.6-4", repos="http://cran.us.r-project.org")

# Load required packages
library(RSQLite)
library(raster)
library(sf)
library(terra)
library(DBI)
library(tidyr)
library(abind)     # Probably not used
library(prodlim)   # Probably not used
library(deSolve)   # Probably not used
library(dplyr)
library(ctmcmove)  # Probably not used
library(lwgeom)    # Probably not used
library(rgeos)



########################### Special requirements ###############################
###
# Please also ensure that the SQLite and SpatiaLite system libraries are installed.

# Define here the path to the folder containing the dynamic link library files
# (e.g., "mod_spatialite.dll", "libcrypto-3-x64.dll", "libsqlite3-0.dll")
# that allow to dynamically loading SpatiaLite as an extension module.
# Please note: folder path should be written with backslash in Windows style
# but since in R the \ is a special character, reserved to escape the character
# that follows it, you have to write it twice.
spatialite_path <- "C:\\spatialite"

# The following line of code will then add the path to this important folder to
# the 'PATH' section of the environment variables of the current R session
Sys.setenv(PATH=paste(Sys.getenv('PATH'), spatialite_path, sep=";"))



############################ Input/output files ################################
###
# Define here the coordinate system CRS as epsg code (numeric)
epsg_code <- 2056

# Filename (and path) of the shapefile of polygons defining the boundaries and
# location of the areas which are covered by forests.
forest_file <- "Forest.shp"

# The filename (and path) of the digital elevation model (DEM) to use.
# The script is conceived to import the DEM in .tif format (i.e., a GeoTIFF with
# location information embedded within the TIFF file). But probably the "rast"
# function of the "terra" package (see below) can accept also other formats.
# Note: For the algorithm to run correctly, the DEM must have equal vertical and
# horizontal resolution (the pixels should be squares).
dem_file <- "DEM.tif"

# The filename (and path) of the output SpatiaLite database.
db_file <- "output.sqlite"



################################# Parameters ###################################
###
# The maximum radius in meters within which to search for a higher point
# upon reaching a local maximum.
rmax <- 100
# When looking for a higher point upon reaching a local maximum, only points
# are considered which are more than mindeltah meters higher than the current
# position.
mindeltah <- 2
# The maximum distance in meters which a fire can travel outside of a forest.
maxnoforest <- 40
# When calculating the trajectory of the fire, the time is incremented by
# tstep (code internal units) every time a new point is calculated. Increasing
# this decreases the precision of the resulting trajectories but increases
# computation speed.
tstep <- 10

# Whether or not to include the heights of the points in the output (slows down
# the calculation)
calc_heights <- T
# Amount of time which passes between two jump checks (i.e., a value of 100 means
# that 100 steps of length tstep are executed between jumps).
lapseduration <- 400
# Only calculate each nth point. If this is set to 1, one ignition point per pixel
# in the DEM is simulated.
calc_nth <- 10
# Maximum number of lapses without significant altitude increase (prevents infinite
# loops due to machine imprecision and slowly crawling fires)
max_noinc <- 2
# Minimum distance in meters that the fire must move between two jump checks in
# order to be able to continue without a jump (i.e., if the fire moves less
# than this, a jump is performed and if no jump is possible the fire is
# extinguished).
mindist <- 2



################################ Load input data ###############################
###
# Load shapefile forest
# forest <- sf::read_sf(dsn= forest_file_path, layer=forest_file_name)
forest <- sf::read_sf(forest_file)

# Load the digital elevation model
# dem <- terra::rast(paste0(dem_file_path, "/", dem_file_name))
dem <- terra::rast(dem_file)
names(dem) <- "h"
set.values(dem)

# Run the external core R script which calculates the steepest path lines.
# Or alternatively open the second script manually and execute it one step at a
# time in order to have more control over the different calculation processes. 
source("calc_lines.r")
