
# Steepest paths tool
# v0.10 (2023-11-11)
# core R script that performs all calculations


# WARNING! This script is very finicky: GIS functions in R and SQLite present
# many possible pitfalls such as differing return types, slight errors in geospatial
# operations etc. If you change something, make sure that this *really* doesn't have
# any effect (f.e. while intersecting with the complement *should* be the same as
# taking the difference, in many cases this is actually *not* the case due to
# rounding errors).


# This file should be sourced from an R script which defines the following variables:
# rmax
# mindeltah
# maxnoforest
# tstep
# forest
# calc_heights
# lapseduration
# calc_nth
# max_noinc
# dem
# mindist
# db_file
# See `main.r` for an explanation of these variables and their possible values.

# Extract the height from the DEM at location x.
h <- function(x) {
  terra::extract(x=dem,y=x,method="bilinear",ID=F)
}

# Compute the gradient of the dem.
# We need to do this manually since on large DEMs the "grad" function
# fails.

iext <- intersect(ext(dem),shift(ext(dem),dx=res(dem)[1],dy=res(dem)[2]))
gr <- (crop(dem,iext)-crop(shift(dem,dx=res(dem)[1],dy=0),iext))/res(dem)[1]
add(gr) <- (crop(dem,iext)-crop(shift(dem,dx=0,dy=res(dem)[2]),iext))/res(dem)[2]
names(gr) <- c("x","y")

gradient <- function(x) {
  terra::extract(x=gr,y=x,method="bilinear",ID=F)
}

# Intersections can be computed more easily if the entire forest is aggregated into a single
# multipolygon.
forest <- st_union(forest)
forest <- st_set_crs(forest,st_crs(dem))

# Square whose sidelength is 2*rmax (approximately a circle of radius rmax)
mat <- matrix(1,ncol=2*round((round(rmax/res(dem)[1])-1)/2)+1,nrow=2*round((round(rmax/res(dem)[2])-1)/2)+1)
# For each pixel in mat compute its (Euclidean) distance to the central pixel
rowc <- t(replicate(ncol(mat),seq(from=1,to=nrow(mat))))
colc <- replicate(nrow(mat),seq(from=1,to=ncol(mat)))
dists <- sqrt((rowc-mean(rowc))**2+(colc-mean(colc))**2)

# For each pixel on the DEM, assuming that a fire is stopped at this point, we
# calculate the pixel to which the fire should jump. This calculation only works
# when the x and y resolutions of the DEM are equal. This is done by jumping to the nearest
# point which is at least mindeltah higher than the current point.
find_jump <- function(X,dists,...) {dists[X<=X[(length(X)+1)/2]+mindeltah] <- NA;t<-arrayInd(which.min(dists),dim(dists));if(length(t)==0){c(NA,NA)} else {t[1,]-arrayInd((length(X)+1)/2,dim(dists))}}
localmaxs <- terra::focal(dem, fun="max", expand=T,w=2*round(min(res(dem))/2)+1+2,silent=F)
localmaxs <- localmaxs==dem
tmpFiles(current=F,orphan=T,old=T,remove=T)
jumppts <- terra::focal(dem, fun=find_jump, dists, w=mat, pad=TRUE, padValue=NA,silent=F)
jumppts$x <- init(jumppts[[1]],"x")
jumppts$y <- init(jumppts[[1]],"y")
jumppts$jx <- jumppts$lyr1*res(dem)[1]+jumppts$x
jumppts$jy <- -jumppts$lyr2*res(dem)[2]+jumppts$y

gc()

# Allow exporting from the database
Sys.setenv(SPATIALITE_SECURITY="relaxed")
db <- dbConnect(RSQLite::SQLite(), db_file)

# The times at which to evaluate the explicit euler method
times <-seq(from=tstep,to=lapseduration-tstep,by=tstep)
# Extract the points at which the steepest paths should start
points <- rasterToPoints(aggregate(raster(dem),fact=calc_nth,fun="mean"),spatial=T) %>% as.data.frame
# maxs contains the maximal height which the steepest path has reached so far. In
# theory, this should always be the final height but due to integration errors the current height may actually
# decrease at certain points in time. maxs is used to detect these points and halt steepest path propagation
# if necessary (or cause a jump, see below).
maxs <- points$h
if(!calc_heights) {
  points <- points %>% dplyr::select(all_of(c("x","y")))
}
points$id <- 1:nrow(points)
# Ignore all steepest paths which started at a local maximum for this timestep.
# This prevents spurious behaviour around the local maximum, where the gradient
# vanishes.
points$ignore <- terra::extract(localmaxs,points[,c("x","y")],ID=F)
points$ignore[is.na(points$ignore)] <- 1
# was_lower contains the number of times that the steepest path was lower after a timestep than
# before it (i.e., integration errors were larger than the increase in altitude.
points$was_lower <- 0

# This is the main loop which simulates the evolution of the steepest paths. I will use this opportunity to
# give a high-level overview of how this loop works:
# The evolution time of the steepest paths is split into time steps of length lapseduration*tstep. These in
# turn are split into substeps of length tstep. During each time step the steepest path evolves by following
# the gradient (using an explit Euler method) without any jumps etc. occurring. Then, after lapseduration many
# steps have passed, a number of checks are performed, jumps are made where necessary and we continue with the
# next time step as before. The procedure between the time steps is as follows:
# * For all steepest paths which did not move at least mindist, perform a jump (or terminate the path if no jump is possible).
# * For all steepest paths whose altitude did not increase by at lease
#   mindeltah, perform a jump (or terminate the path if no jump is possible).
# * For all points which did not experience a jump, perform a single iteration
#   of the explicit Euler method (so that all paths are at the same point in time
#   again).
# * Terminate all lines for which was_lower is larger or equal to max_noinc
# * If the altitude of the steepest path is lower than ints previous maximum
#   plus mindeltah, increment was_lower, else reset was_lower to 0.
# * Mark all steepest paths which are now located at a local maximum to be
#   ignored (see the description of "ignore" above).
# * Continue with the next time step or terminate the loop if all steepest paths were terminated.

it <- 1
while(T) {
  # Clear the memory
  gc()
  tmpFiles(current=F,orphan=T,old=T,remove=T)

  prevpts <- points

  # Write the points to the database
  dbWriteTable(db,"points",points %>% dplyr::mutate(t=(it-1)*lapseduration) %>% dplyr::select(-"ignore"),append=T)
  # Model the steepest path as following the gradient up until the next time step (next jump check)
  for(j in seq_along(times)) {
    print(times[j])
    # Explicit Euler
    points[points$ignore==F,c("x","y")] <- points[points$ignore==F,c("x","y")]+tstep*gradient(points[points$ignore==F,c("x","y")])
    if(calc_heights) {
      points[,"h"] <- h(points[,c("x","y")])
    }
    dbWriteTable(db,"points",points %>% dplyr::mutate(t=times[j]+(it-1)*lapseduration) %>% dplyr::select(-"ignore"),append=T)
  }

  # This array contains 0 if the steepest path does not need to jump and 3 if
  # it does need to jump (yes, I know, this convention is bizzarre for
  # historical reasons and I'm afraid to touch it).
  conn <- array(data=0,dim=nrow(points))
  #conn[which(is.na(endpts[,"x"])|is.na(endpts[,"y"]))] <- 4

  # Distance travelled since last time step
  dist <- sqrt(rowSums((points[,c("x","y")]-prevpts[,c("x","y")])**2))
  # Didn't move mindist meters
  conn[which(dist<mindist)] <- 3
  rm(dist)
  # Didn't increase by at least mindeltah (also prevents loops)
  conn[h(points[,c("x","y")])<=maxs+mindeltah] <- 3

  # Do not jump from NA points
  conn[is.na(points$h)] <- 0

  # Steepest paths which need to jump
  validpts <- points[which(conn==3),c("x","y")]

  # Points to jump to
  jnd<-terra::extract(jumppts,validpts %>% st_as_sf(coords=c("x","y"),crs=crs(dem)))
  jumps <- cbind(validpts,dplyr::select(jnd,all_of(c("jx","jy"))))
  rm(jnd)

  # Perform a jump or an Euler step.
  points[which(conn==3),c("x","y")] <- jumps[,c("jx","jy")]
  points[which(conn!=3),c("x","y")] <- points[which(conn!=3),c("x","y")]+tstep*gradient(points[which(conn!=3),c("x","y"),drop=F])
  points[points[,"was_lower"]>=max_noinc,"x"] <- NA

  # Points which are still valid
  nonnull <- which(!is.na(points[,"x"]))
  if(length(nonnull)==0) {
    break
  }
  hts <- h(points[,c("x","y")])
  points[which(hts<=maxs+mindeltah),"was_lower"] <- points[which(hts<=maxs+mindeltah),"was_lower"]+1
  points[which(hts>maxs+mindeltah),"was_lower"] <- 0
  # Update maxs
  maxs <- pmax(maxs, hts %>% unlist)[nonnull]
  points <- points[nonnull,]
  prevpts <- points
  # Update ignore
  points$ignore <- terra::extract(localmaxs,points[,c("x","y")],ID=F)
  points$ignore[is.na(points$ignore)] <- T

  # Number of steepest paths remaining
  print(length(nonnull))
  it <- it+1
}

# At this point, all steepest paths have been calculated. However, they are stored as a collection of points in the database.
# We now need to aggregate these points into linestrings and intersect them with the forest in order to obtain the blocked and
# last_forest lines.

# Load the dynamic link libraries (.dll) consituting the SpatiaLite extension module
# from the folder path "spatialite_path" defined in the main R script.
# Please note: this path should be already included among the environment variables of the R session
dbExecute(db,"select load_extension('mod_spatialite')")

# Following alternative formulations do not work:
# dbExecute(db,"select load_extension('C:/spatialite/mod_spatialite.dll')")
# dbExecute(db,"select load_extension('C:\\spatialite\\mod_spatialite.dll')")

# Initialize the SpatiaLite extension module
dbExecute(db,"select initspatialmetadata()")

# valid is 1 if the line does not touch the border.
dbExecute(db,"create table lines(id int,valid int default 0,primary key(id))")
dbExecute(db,paste0("select addgeometrycolumn('lines','line',",epsg_code, ",'LINESTRING','XYZ')"))
# Note that we should order by t. But since we are inserting the values ordered by t it is
# highly probable that they are already sorted by t. It would be good to check the impact
# of the ordering and add it if it not too expensive.
dbExecute(db,paste0("insert into lines(id,line) select id,makeline(makepointz(x,y,h,",epsg_code,")) from points group by id"))
wkb <- ext(dem) %>% as.polygons(crs=crs(dem)) %>% st_as_sf %>% st_buffer(dist=-2*max(res(dem))) %>% st_as_sfc %>% st_as_binary(hex=T)
dbExecute(db,paste0("update lines set valid=1 where contains(polyfromwkb(casttoblob('",wkb,"',1),",epsg_code,"),line)"))

# Get the maximal extent of any steepest path in x and y direction.
mxd <- ceiling(as.numeric((dbGetQuery(db,"select max(MbrMaxX(line)-MbrMinX(line)) from lines"))))
myd <- ceiling(as.numeric((dbGetQuery(db,"select max(MbrMaxY(line)-MbrMinY(line)) from lines"))))
# Add the maximal extent of each line as a column (maximum of x and y extent).
dbExecute(db,paste0("alter table lines add maxdim float"))
dbExecute(db,"update lines set maxdim=max(MbrMaxY(line)-MbrMinY(line),MbrMaxX(line)-MbrMinX(line))")

# The block of code which follows has the sole purpose of intersecting the lines with the forest in order
# to determine the intersection segments, whether the end is contained in the forest, etc. It is a mess (although
# unfortunately a necessary one), pray you do not need to change it. The general idea is to cut the forest into
# squares (initially of a side length of 2000 m) and then process all lines which are completely contained within
# this square. Since some lines may be longer than 2000 m, we then repeat this procedure for squares with a
# side length of 3000 m, etc. until all lines are processed.
dbExecute(db,paste0("alter table lines add row int"))
dbExecute(db,paste0("alter table lines add col int"))
dbExecute(db,"create index row_col_idx_lines on lines (row,col)")
forestbdry <- st_boundary(forest) %>% st_cast("LINESTRING")
forestbdrybboxs <- sapply(forestbdry,st_bbox)
forest <- forest %>% st_cast("POLYGON")
forestbboxs <- sapply(forest,st_bbox)
dbExecute(db,"create table intersections(id int,foreign key (id) references lines(id))")
dbExecute(db,paste0("select addgeometrycolumn('intersections','intersections',",epsg_code,",'MULTIPOINT','XYZ')"))
dbExecute(db,paste0("alter table lines add end_contained_in_forest bool"))
dbExecute(db,paste0("alter table lines add intersects_forest bool"))
dbExecute(db,paste0("alter table lines add processed int default 0"))
dbExecute(db,"create index processed_lines on lines (processed)")
dbExecute(db,"create index processed_row_col_lines on lines (processed,row,col)")
# Boundary lines
dbExecute(db,"delete from lines where line is null")
width=2000
height=2000
while(dbGetQuery(db,"select count(*) from lines where processed==0")[[1]]>0) {
  print(width)
  dbExecute(db,paste0("update lines set processed=1 where processed=0 and maxdim<",width/2))
  # +1 b.c. R indices start at 1
  dbExecute(db,paste0("update lines set row=floor((MbrMinY(line)-",ymin(dem),")/",height/2,")+1 where processed=1"))
  dbExecute(db,paste0("update lines set col=floor((MbrMinX(line)-",xmin(dem),")/",width/2,")+1 where processed=1"))
  idxs=dbGetQuery(db,"select col,row from lines where processed=1 and col>0 and row>0") %>% unique %>% arrange(col,row)
  xs=seq(from=xmin(dem),to=xmax(dem),by=width/2)
  ys=seq(from=ymin(dem),to=ymax(dem),by=height/2)
  for(r in rownames(idxs)) {
    i=idxs[r,"col"]
    j=idxs[r,"row"]
    print(paste0("(",i,",",j,")"))
    ijforest <- st_union(st_crop(forestbdry[(forestbdrybboxs["xmin",]<=xs[i]+width & forestbdrybboxs["xmax",]>=xs[i] & forestbdrybboxs["ymin",]<=ys[j]+height & forestbdrybboxs["ymax",]>=ys[j])],c(xmin=xs[i],ymin=ys[j],xmax=xs[i]+width,ymax=ys[j]+height)))
    if(length(st_geometry_type(ijforest))>0&as.character(st_geometry_type(ijforest))[1]!="GEOMETRYCOLLECTION") {
      ijforest <- ijforest %>% st_geometrycollection
    }
    ijforest <- ijforest %>% st_collection_extract(type="LINESTRING") %>% st_cast("MULTILINESTRING") %>% st_union

    ijforestint <- st_union(st_crop(forest[(forestbboxs["xmin",]<=xs[i]+width & forestbboxs["xmax",]>=xs[i] & forestbboxs["ymin",]<=ys[j]+height & forestbboxs["ymax",]>=ys[j])],c(xmin=xs[i],ymin=ys[j],xmax=xs[i]+width,ymax=ys[j]+height)))
    # Type problems
    if(length(st_geometry_type(ijforestint))>0&as.character(st_geometry_type(ijforestint))[1]!="GEOMETRYCOLLECTION") {
      ijforestint <- ijforestint %>% st_geometrycollection
    }
    ijforestint <- ijforestint %>% st_collection_extract() %>% st_cast("MULTIPOLYGON") %>% st_union
    # The points at which the line intersects the forest boundary.
    dbExecute(db,paste0("insert into intersections(id,intersections) select id,casttomultipoint(st_intersection(line,mlinefromwkb(casttoblob('",st_as_binary(ijforest,hex=T),"',1),",epsg_code,"))) from lines where processed=1 and col=",i," and row=",j))
    # Whether the end of the line is contained in the forest
    dbExecute(db,paste0("update lines set end_contained_in_forest=st_contains(mpolyfromwkb(casttoblob('",st_as_binary(ijforestint,hex=T),"',1),",epsg_code,"),endpoint(line)) where processed=1 and col=",i," and row=",j))
    # Does the line intersect the forest?
    dbExecute(db,paste0("update lines set intersects_forest=st_intersects(mpolyfromwkb(casttoblob('",st_as_binary(ijforestint,hex=T),"',1),",epsg_code,"),line) where processed=1 and col=",i," and row=",j))
  }
  dbExecute(db,paste0("update lines set processed=2 where processed=1"))
  width=width+1000
  height=width+1000
}

# Split the lines at the intersection points with the forest.
dbExecute(
          db,paste0(
                    "create table lines_split as
                    SELECT id, casttomultilinestring(snapandsplit(line,intersections,0.1)) as geom
                    FROM intersections inner join lines using (id)"
          )
)
dbExecute(db,paste0("SELECT RecoverGeometryColumn('lines_split','geom',",epsg_code,",'MULTILINESTRING','XYZ');"))
dbExecute(db,paste0("update lines_split set geom=casttomultilinestring(t.line) from (select * from lines) as t where t.id=lines_split.id and lines_split.geom is NULL"))
dbExecute(db,paste0("select elementarygeometries('lines_split','geom','split_segs','out_pk','out_multi_id')"))
dbExecute(db,paste0("alter table split_segs add noforest_len decimal"))

# The next block of code is the second mess. Now that we have split the lines into segments, we need to find
# out which ones lie within and which ones lie outside the forest. To do so, for each segment we compute the
# length of its intersection with the forest.
dbExecute(db,paste0("alter table split_segs add row int"))
dbExecute(db,paste0("alter table split_segs add col int"))
dbExecute(db,paste0("alter table split_segs add outside_forest bool"))
dbExecute(db,paste0("alter table split_segs add maxdim float"))

dbExecute(db,"update split_segs set maxdim=max(MbrMaxY(geom)-MbrMinY(geom),MbrMaxX(geom)-MbrMinX(geom))")
dbExecute(db,paste0("alter table split_segs add processed int default 0"))
dbExecute(db,"create index row_col_idx_segs on split_segs (row,col)")
dbExecute(db,"create index processed_split_segs on split_segs (processed)")
dbExecute(db,"create index processed_row_col_split_segs on split_segs (processed,row,col)")

width=2000
height=2000

while(dbGetQuery(db,"select count(*) from split_segs where processed=0")[[1]]>0) {
  print(width)
  dbExecute(db,paste0("update split_segs set processed=1 where processed=0 and maxdim<",width/2))
  # +1 b.c. R indices start at 1
  dbExecute(db,paste0("update split_segs set row=floor((MbrMinY(geom)-",ymin(dem),")/",height/2,")+1 where processed=1"))
  dbExecute(db,paste0("update split_segs set col=floor((MbrMinX(geom)-",xmin(dem),")/",width/2,")+1 where processed=1"))
  idxs=dbGetQuery(db,"select col,row from split_segs where processed=1 and col>0 and row>0") %>% unique %>% arrange(col,row)
  xs=seq(from=xmin(dem),to=xmax(dem),by=width/2)
  ys=seq(from=ymin(dem),to=ymax(dem),by=height/2)
  for(r in rownames(idxs)) {
    i=idxs[r,"col"]
    j=idxs[r,"row"]
    print(paste0("(",i,",",j,")"))
    ijforest <- st_union(st_crop(forest[(forestbboxs["xmin",]<=xs[i]+width & forestbboxs["xmax",]>=xs[i] & forestbboxs["ymin",]<=ys[j]+height & forestbboxs["ymax",]>=ys[j])],c(xmin=xs[i],ymin=ys[j],xmax=xs[i]+width,ymax=ys[j]+height)))
    if(length(st_geometry_type(ijforest))>0&as.character(st_geometry_type(ijforest))[1]!="GEOMETRYCOLLECTION") {
      ijforest <- ijforest %>% st_geometrycollection
    }
    ijforest <- ijforest %>% st_collection_extract(type="POLYGON") %>% st_cast("MULTIPOLYGON") %>% st_union
    if(length(ijforest)>0) {
      dbExecute(db,paste0("update split_segs set noforest_len=st_length(st_difference(geom,mpolyfromwkb(casttoblob('",st_as_binary(ijforest,hex=T),"',1),",epsg_code,"))) where processed=1 and col=",i," and row=",j))
    } else {
      dbExecute(db,paste0("update split_segs set noforest_len=st_length(geom) where processed=1 and col=",i," and row=",j))
    }
  }
  dbExecute(db,paste0("update split_segs set processed=2 where processed=1"))
  width=width+1000
  height=width+1000
}

# Now we are almost done! It remains to determine the blocked_lines and last_forest_lines and compute the lengths.

# Split the line at the first point at which it leaves the forest for at least maxnoforest meters.
dbExecute(db,paste0("create table blocked_lines as select id,line_substring(line,0,min(line_locate_point(line,startpoint(geom)))) as line from lines inner join split_segs using (id) where noforest_len>=",maxnoforest," group by id"))
dbExecute(db,paste0("SELECT RecoverGeometryColumn('blocked_lines','line',",epsg_code,",'LINESTRING','XYZ');"))
# Add all the lines which are completely contained in the forest.
dbExecute(db,paste0("insert into blocked_lines(id,line) select id,line from (select id,line,intersects_forest from lines union all select id,lines.line,intersects_forest from blocked_lines inner join lines using (id)) as t group by id having count(*)=1 and sum(intersects_forest)=1"))
# Add all lines to the last_forest_lines whose end is contained in the forest
dbExecute(db,paste0("create table last_forest_lines as select id,line from lines where end_contained_in_forest=TRUE"))
# For all other lines, add the portion of the line up until the last time it leaves the forest.
dbExecute(db,paste0("insert into last_forest_lines(id,line) select id,line_substring(line,0,max(line_locate_point(line,startpoint(geom)))) from split_segs inner join lines using (id) where end_contained_in_forest=FALSE and intersects_forest=TRUE and st_length(split_segs.geom)>=0.2 group by id"))
dbExecute(db,paste0("SELECT RecoverGeometryColumn('last_forest_lines','line',",epsg_code,",'LINESTRING','XYZ');"))
# Shorter than 0.2 and contained in noforest
dbExecute(db,paste0("insert into last_forest_lines(id,line) select id,NULL from (select id,line from lines union all select id,lines.line from last_forest_lines inner join lines using (id)) as t group by id having count(*)=1"))

# Lengths
dbExecute(db,paste0("alter table lines add length decimal"))
dbExecute(db,paste0("alter table blocked_lines add length decimal"))
dbExecute(db,paste0("alter table last_forest_lines add length decimal"))
dbExecute(db,paste0("update lines set length=st_length(casttoxy(line))"))
dbExecute(db,paste0("update blocked_lines set length=st_length(casttoxy(line))"))
dbExecute(db,paste0("update last_forest_lines set length=st_length(casttoxy(line))"))

dbExecute(db,paste0("alter table lines add length3d decimal"))
dbExecute(db,paste0("alter table blocked_lines add length3d decimal"))
dbExecute(db,paste0("alter table last_forest_lines add length3d decimal"))
dbExecute(db,paste0("update lines set length3d=st_3dlength(line)"))
dbExecute(db,paste0("update blocked_lines set length3d=st_3dlength(line)"))
dbExecute(db,paste0("update last_forest_lines set length3d=st_3dlength(line)"))

