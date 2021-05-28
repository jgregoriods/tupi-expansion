library(rasterVis)
library(rgdal)

source("utils.R")

coast <- readOGR("shp/south_america.shp")

biomes <- raster("layers/biomes.asc")
proj4string(biomes) <- wgs
speed <- frontSpeed(0.025, 50, 30)
origin <- c(-61.9574, -10.9557)

cost.null <- biomes
cost.null[values(cost.null) > 1] <- 1
iso.null <- simulateDispersal(cost.null, origin, 5000, speed)

cost.forest <- biomes
cost.forest[values(cost.forest) > 1] <- 4
iso.forest <- simulateDispersal(cost.forest, origin, 5000, speed)

nm <- raster("res/null_model.asc")
proj4string(nm) <- albers

fm <- raster("res/forest_model.asc")
proj4string(fm) <- albers

nm <- projectRaster(nm, iso.null)
fm <- projectRaster(fm, iso.null)

levelplot(stack(iso.null, iso.forest), at=500*1:10, col.regions=viridis(9),
          names.attr=c("a", "b"), scales=list(x=list(draw=FALSE),
          y=list(draw=F)), xlab="", ylab="",
          colorkey=list(height=0.75, width=1)) + layer(sp.polygons(coast))
dev.print(jpeg, "disperse.jpg", width=1200, height=900, res=300)
dev.off()

levelplot(stack(nm, fm), at=500*1:10, col.regions=viridis(9),
          names.attr=c("a", "b"), scales=list(x=list(draw=FALSE),
          y=list(draw=F)), xlab="", ylab="",
          colorkey=list(height=0.75, width=1)) + layer(sp.polygons(coast))
dev.print(jpeg, "sim.jpg", width=1200, height=900, res=300)
dev.off()
