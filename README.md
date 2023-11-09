# steepest_paths_tool

Description
-----------
**Steepest paths tool** is an R script that calculates three-dimensional polyline vectors representing a first approximation of the trajectories of propagation of forest fires on mountain slopes in the absence of strong winds solely on the basis of a digital elevation model and a vector map of forested areas.

One of the main purposes of this tool is to provide a simplified representation of the possibilities of fire spread on mountain slopes which can be implemented from readily available GIS data (i.e., a DEM and a shapefile representing forest cover) without considering more complex data such as the influence of particular weather circumstances or the presence of different fuel types and conditions in the territory.
In this sense the tool is relatively economical in terms of imput data and computing resources and is therefore suitable for application over vast areas. 

Tool outputs can be used for instance to assess the potential of fire propagation across landscapes or to estimate the risk of large fires at regional level by taking into account at the same time the spatial distribution of forest cover and the shape of the terrain surface.

The tool is based on the assumption that in mountainous environment, the spread of fire is mainly determined by the combined influence of topography and fuel connectivity, that is, by the spatial arrangement of the forests and open areas on the slope in relation to the orientation of the lines of maximum gradient.

Definitions
----------------------
With reference to the spread of forest fires on mountain slopes, we define the steepest path as the 3D polyline or linestring that starts from the point of ignition, runs up the slope following the direction of the maximum gradient of the terrain surface (i.e. gaining the most elevation with the shortest possible path), always moves on land covered by woody vegetation (trees or shrubs), and stops where the path reaches a main summit (e.g., a mountain top) or encounters large open areas that can arrest the spread of fire (e.g., grassland, farmland, rock faces).

Outputs
----------------------
The output of the steepest paths tool is a SpatiaLite database which can be accessed using either the open source geographic information system QGIS or the R statistical software (with “rgdal” or “sf” packages).
This database contains 3 vector layers of three-dimensional polylines, which can be described as follow:
* **lines**: The entire computed line from the point of ignition until a main mountain summit is reached, calculated solely on the basis of the morphology of the terrain and without considering the intersections with the forest cover.
* **blocked_lines**: The computed line up until the point at which for the first time the trajectory leaves the forest for a horizontal distance greater than a threshold value (see the _maxnoforest_ variable which is set by default to 40 meters in the R script).
* **last_forest_lines**: The computed line up until the point at which the trajectory leaves the forested area for the last time.

Common attributes assigned to each resulting vector feature:
* **id**: A unique number identifying the line. This number is determined by the starting point of the line which in our case is intended to represent a possible ignition point of a forest fire. These starting points are located on a regular grid by subsampling the centroids of the pixels of the digital elevation model (see the _calc_nth_ variable which is set by default to 10 in the R script). For each ignition point the steepest paths algorithm produce three types of polylines (lines, blocked_lines, last_forest_lines) which share the same origin and have therefore the same identification number.
* **valid**: Can the steepest path obtained be considered valid or not? If the line reached the edge of the DEM raster at some point during its evolution, this value is 0, else 1.
* **end_contained_in_forest**: If the endpoint of the line is contained within the forest, this value is 1, else 0.
* **intersects_forest**: If the line intersects the forest, this value is 1, else 0.
* **length**: The 2 dimensional length of the line (in meters).
* **length3d**: The three-dimensional length of the line (in meters).

Challenges
----------------------
For an overview of the main challenges that were faced in the attempt to optimise the steepest path algorithm, read the [document](https://github.com/Insubric/steepest_paths_tool/blob/master/Challenges.md) dedicated to this topic.

Comparison with other tools
----------------------
The _r3.flow_ function of [GRASS GIS](https://grass.osgeo.org/) software uses a digital elevation model as imput data to compute 3D flow lines both downstream (with parameter _direction_ set to "down") and upstream (with parameter _direction_ set to "up").
Upstream 3D flow lines resulting from this GRASS GIS function are quite similar to the steepest paths.

However, the **Steepest paths tool** has some features that may be advantageous for some users:
1) Is written in R and can therefore be more easily customised to particular needs.
2) It offers an integrated module for not interrupting the upward progression on the terrain surface of the resulting three-dimensional polylines when they reach a local maximum of little importance in the digital elevation model.
3) It calculated a whole series of attributes and descriptors for each resulting vector feature.

Contact
----------------------
Jeremy Feusi (e-mail: jeremy.feusi@wsl.ch or jeremy@feusi.co)

Citation
----------------------
If you use the **Steepest paths tool** and publish some result based on this R script, please cite this software as follows:

Swiss Federal Institute for Forest Snow and Landscape Research WSL, 2023, Steepest paths tool (vers. 1.0), https://github.com/Insubric/steepest_paths_tool.
