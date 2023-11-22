# steepest_paths_tool

Description
-----------
The **Steepest paths tool** is an R script which computes three-dimensional polyline vectors representing a first approximation of the trajectories of propagation of forest fires along mountain slopes in the absence of strong winds, operating solely on the basis of a digital elevation model and a vector map of the forested areas.

One of the main purposes of this tool is to provide a simplified representation of the potential forest fire spread on mountain slopes, which can be implemented from readily available GIS data (i.e., a DEM and a shapefile representing forest cover) without taking into consideration more complex data such as the influence of particular weather conditions or the presence of different fuel types in the territory.
In this sense the tool is rather economical in terms of input data and computational resources and may therefore be applied even over large areas. 

The outputs of the tool may for instance be used to assess the potential fire propagation across a landscape or to estimate the risk for large fires to occur on a regional level, by at the same time taking into account the spatial distribution of the forest cover and the shape of the terrain surface.

The implementation is based on the assumption that in mountainous environments the spread of a fire is mainly determined by the combined influence of topography and fuel connectivity, that is, by the spatial interplay between the forest and the slope. In other words, we consider the spatial arrangement of the forests and open areas on the slope in relation to the orientation of the lines of maximum gradient as a key factor for the assessment of the potential of spreading of forest fires along mountain slopes.

Definitions
----------------------
Regarding the spread of forest fires on mountain slopes, we define the _steepest path_ to be the 3D polyline or linestring which starts at the point of ignition and then proceeds uphill, following the direction of the maximum gradient of the terrain surface (i.e. gaining the most elevation with the shortest possible path), while always remaining within the land covered by woody vegetation (trees or shrubs), and finally terminating where the path reaches a major summit (e.g., a mountain top) or encounters a large open area which can prevent the further spread of fire (e.g., grassland, farmland, rock faces).

Outputs
----------------------
The output of the steepest paths tool is a SpatiaLite database which can be read using either the open source geographic information system QGIS or the R statistical software (with “rgdal” or “sf” packages).
This database contains 3 vector layers of three-dimensional polylines, which can be described as follow:
* **lines**: The entire computed line from the starting point (representing a possible ignition point of a forest fire) to a main mountain summit is reached, calculated solely on the basis of the morphology of the terrain and without taking into consideration the intersection with the forest cover.
* **blocked_lines**: The computed line up until the point at which the trajectory leaves the forest for a horizontal distance greater than a threshold value for the first time (see the _maxnoforest_ variable which is set to 40 meters by default in the main R script).
* **last_forest_lines**: The computed line up until the point at which the trajectory leaves the forested area for the last time.

Common attributes assigned to each resulting vector feature:
* **id**: A unique number identifying each line. This number is determined by the starting point of the line which in our case is intended to represent a possible ignition point of a forest fire. These starting points are generated at regular distances on the digital elevation model by dividing its extent longitudinally and latitudinally by a multiple of the pixel size (see the _calc_nth_ variable which is set to 10 by default in the main R script) so that each of these points is located precisely midway between four pixel centroids. For each ignition point the steepest paths algorithm produce three types of polylines (lines, blocked_lines, last_forest_lines) which share the same origin and therefore have the same identification number.
* **valid**: Can the computed steepest path be considered valid or not? If the line reached the edge of the DEM raster at some point during its evolution, this value is 0, else 1.
* **end_contained_in_forest**: If the endpoint of the line is contained within the forest, this value is 1, else 0.
* **intersects_forest**: If the line intersects the forest, this value is 1, else 0.
* **length**: The two-dimensional length of the line (in meters).
* **length3d**: The three-dimensional length of the line (in meters).

Challenges
----------------------
For an overview of the main challenges that were faced in the attempt to optimise the steepest path algorithm, read the [document](https://github.com/Insubric/steepest_paths_tool/blob/master/Challenges.md) dedicated to this topic.

Comparison with other tools
----------------------
The _r3.flow_ function of the [GRASS GIS](https://grass.osgeo.org/) software uses a digital elevation model as input data to compute 3D flow lines both downstream (with parameter _direction_ set to "down") and upstream (with parameter _direction_ set to "up").
Upstream 3D flow lines resulting from this GRASS GIS function are quite similar to the steepest paths.

Compared to an established and highly reliable tool such as this GRASS GIS function, our ad hoc project is in many ways still unfinished and there is certainly much room for improvement.
However, the **Steepest paths tool** has some features that may be advantageous for some users:
1) It is written in R and can therefore be more easily customised to particular needs.
2) It offers an integrated module for the continuation of the upward progression on the terrain surface when a local maximum of little importance is reached (see the paragraph devoted to "jumps" in [this document](https://github.com/Insubric/steepest_paths_tool/blob/master/Challenges.md)).
3) It computes some useful attributes for each resulting vector feature.
4) It provides the possibility to analyze the intersection of the steepest paths with the forested areas.

Requirements
----------------------
R, RTools and eventually also RStudio Desktop, as well as a selection of R packages (see the section "Required packages" in the main R script) must be installed and up to date on the target computer.

The SQLite and SpatiaLite system libraries must also be installed on your system. In Linux, these libraries may be installed using the package manager.

For Windows on the other hand, all dynamic link libraries (.dll) consituting the SpatiaLite extension module (e.g., "mod_spatialite.dll", "libcrypto-3-x64.dll", "libsqlite3-0.dll", "libgeos.dll", "libproj_9_2.dll") must be available in a folder included among the environment variables of the R session (see "Special requirements" section in the main R script).

These .dll files can be downloaded as MS Windows binaries from the [Gaia-SINS](https://www.gaia-gis.it/gaia-sins) federated projects home-page.
See for instance the content of the "mod_spatialite-5.1.0-win-amd64.7z" compressed archive file (vers. 2023-08-05). This SpatiaLite extension module is a pure loadable module lacking any explicit SQLite3 dependency (see the specific [web page](https://www.gaia-gis.it/fossil/libspatialite/wiki?name=mod_spatialite) of the Gaia-SINS for further explanations). 

Installation instructions
----------------------
Download from the "Releases" section and copy/paste to a desired folder location on your computer the two R scripts which constitute the current version (v0.10) of the Steepest paths tool:
* The main script "**main.R**" which performs all preparatory operations prior to calculations (i.e., the setting of all parameters that control the calculations and loading of input data).
* The core script "**calc_lines.R**" which calculates the steepest path lines and the intersection with forest polygons.

Add in the same folder the following two input data:
* The **digital elevation model** (DEM) in .tif format which represents the elevation of the terrain surface of your region of interest. Obviously, there must be some fairly pronounced relief (hills or mountains) since on a flat surface such as a wide plain steepest paths fail to develop. The DEM should be a plane grid with square pixels and use a projected coordinate system which form a Cartesian reference system. 
* The **shapefile of polygons** which represent the forest areas in your region of interest. This shapefile and the DEM must have the same and well-defined coordinate system (CRS). The tool is set up to work with coordinate systems that use the meter as reference unit of measure.

Getting started tips
----------------------
Launch RStudio, create a new empty project (File -> New Project) and save it in the same folder where the above-mentioned 4 items are already found.

Add the main and core scripts ("main.R", "calc_lines.R") to this newly created RStudio project (File -> Open File). If you don't want to use RStudio or even create a new project remember to set the working directory properly (function _setwd_).

Run the main R script. Attention: at least the first time, it is best to proceed one step at a time so that you understand the organization of the code and read the included comments and instructions.

In the sections "**Special requirements**" and "**Input/output files**" the user has to enter and define some file and path names as well as the coordinate system.
In the following section entitled "**Parameters**" the user can set and modify some important variables that control the computation of the steepest paths.

When executing the last line of the main R script, the function _source_ trigger the run of the joined core R script which calculates the steepest path lines. Alternatively the user can skip this last line and execute the core script manually one step at a time in order to have more control over the different calculation processes.

Be aware that the calculation can take a long time (many hours or even days) depending on the size of the DEM, the number of starting points generated and the processing resources available in the computer. We therefore recommend starting with a relatively small DEM so as to essay the calculation time. For example, using a DEM with an extent of 20 km x 20 km and a pixel size of 10 m and creating a starting point every 10 pixels (i.e., with _calc_nth_ set to 10 in the main R script) the calculation could last about twenty minutes and you already get as a result 40,000 output polyline features (number of output lines = DEM size / pixel size / how many pixels for a single starting point = 20,000 * 20,000 m / 10 * 10 m / 10 * 10 pixels = 400,000,000 / 100 / 100 = 40,000).

Contact
----------------------
Jeremy Feusi (e-mail: jeremy.feusi@wsl.ch)

Citation
----------------------
If you use the **Steepest paths tool** and publish some result based on this R script, please cite this software as follows:

Swiss Federal Institute for Forest Snow and Landscape Research WSL, 2023, Steepest paths tool (v0.10), https://github.com/Insubric/steepest_paths_tool.
