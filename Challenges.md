# Challenges in optimizing the steepest paths algorithm

Introduction
-----------
With the aim to obtain a simplified model of forest fire propagation in mountainous regions, we define the steepest path as the 3D polyline or linestring that starts from the point of ignition, runs up the slope following the direction of the maximum gradient of the terrain surface (i.e. gaining the most elevation with the shortest possible path), always moves on land covered by woody vegetation (trees or shrubs), and stops where the path reaches a main summit (e.g., a mountain top) or encounters large open areas that can arrest the spread of fire (e.g., grassland, farmland, rock faces).

This document presents an overview of the technical challenges of the steepest paths algorithm. This should be considered by anyone wishing to use the code or implement a similar method.

Trajectory computation
----------------------
Denote the position of the fire at time **t** by **x(t)**. The procedure outlined in the submitted paper (Conedera et al. 2023) makes use of the fact that the path which follows the gradient of the slope is exactly the path for which

**ẋ(t)=∇h**

where **h** is the altitude. In other words, the velocity of **x** is collinear with the gradient at the current location. Starting from the user-defined point of ignition, we thus may formulate the problem of determining **x** as an initial value problem (**IVP**).
We numerically resolve this equation by means of an explicit Euler method:

**x(t+Δt)≈x(t)+Δt·∇h**

Note that this is a first order method and as such, a low degree of accuracy of the obtained results is to be expected. However, in practice this proved to be sufficient and was thus preferred over more complex methods due to its low computational complexity.

Vectorization
----------------------
Much of the code design was governed by the constraints stemming from the fact that on the one hand the **IVP** needed to be solved for a large number of initial values (several millions) while at the same time, the entire codebase was required to be written in R, in order to ensure a high degree of maintainability and customizability for a wide range of possible end-users.
This meant that much of the code was necessarily written in a vectorized form. In particular, instead of computing the trajectories separately for each point of ignition, we simultaneously update the positions for all simulated fire trajectories, leveraging the extract function from the terra package, which supports simultaneously extracting the values of a raster at multiple points.

Memory usage
----------------------
Another constraint is posed by the finite amount of memory available for computation. The large amount of three-dimensional location data produced by the algorithm may exhaust this memory when operated across large areas.

To resolve this, we incrementally transfer the computed trajectories to an SQLite database, which can either be stored in-memory (for speed) or written to the hard disk.
Subsequently, we avoid reading the data back into the memory by performing all computations directly within the database using the SpatiaLite module.

Computation of the jumps
----------------------
Usually, digital elevation models contain lots of local elevation maxima that could trap the steepest paths preventing them from developing further up the slope.

Considering that forest fires usually only stop when they reach a summit that is truly prominent in relation to the surrounding relief, we have implemented a procedure to allow the path to jump from a local maximum to the nearest point located at a higher elevation on the slope within a certain distance (see the rmax variable which is set by default to 100 meters in the R script).

Regarding these “jumps” (see Fig. 2 and related description in the submitted paper), observe that many ignition points may converge towards a common local maximum (e.g., the summit of a small hill).
Hence separately computing the point to which each fire jumps will inevitably lead to duplicate computations.

Therefore, prior to commencing the computation of the trajectories, we instead compute for each pixel the point to which a fire would jump assuming it reaches a local maximum at this location.
During the computation of the trajectories, we thus avoid the expensive computation of jump points.

Structure of the algorithm
----------------------
The structure of the algorithm can hence be summarized as follows:

1. Compute the gradient from the DEM.
2. Compute the “jump points” as described in the previous paragraph.
3. Compute the steepest paths following the procedure described in the submitted paper.
4. Truncate the steepest paths at the location at which for the first time they leave the forest for a horizontal distance greater than a threshold value (see the maxnoforest variable which is set by default to 40 meters in the R script).

Computation of the truncation of the trajectories
----------------------
Computationally, point 4 is the most expensive operation, since it requires intersecting the computed trajectories with the polygons of the forested area. This can be expensive if the study area is large and the geometry of the forest is complex.

To ameliorate this situation, instead of directly computing the intersection of each path with the entirety of the forest, we instead split the forest into squares of a given initial side length.
Subsequently, for each square, we compute the intersections with the trajectories which are completely contained within it.

Clearly, some trajectories may be too long to be contained in a single square. Therefore, we repeat the procedure with a larger square size, until all trajectories are contained in some square.
In this manner, at each step the computation becomes more expensive as the square size grows but at the same time, the number of trajectories with which the forest must be intersected diminishes, since the shorter paths will have already been intersected with squares of smaller side length.

Reference
----------------------
Conedera, Marco; Feusi, Jeremy; Pezzatti, Gianni Boris; Krebs, Patrik; 2023. Will future risk of large fires in mountain areas depend on bottom-up fire drivers? Paper submitted to _Natural Hazards_ (Springer) in November 2023.
