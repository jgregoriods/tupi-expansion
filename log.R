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

plt <- levelplot(stack(iso.null, iso.forest), at=500*1:10, col.regions=viridis(9),
          names.attr=c("a", "b"), scales=list(x=list(draw=FALSE),
          y=list(draw=F)), xlab="", ylab="",
          colorkey=list(height=0.75, width=1)) + layer(sp.polygons(coast))
plot(plt)
dev.print(jpeg, "disperse.jpg", width=1200, height=900, res=300)
dev.off()

plt <- levelplot(stack(nm, fm), at=500*1:10, col.regions=viridis(9),
          names.attr=c("a", "b"), scales=list(x=list(draw=FALSE),
          y=list(draw=F)), xlab="", ylab="",
          colorkey=list(height=0.75, width=1)) + layer(sp.polygons(coast))
plot(plt)
dev.print(jpeg, "sim.jpg", width=1200, height=900, res=300)
dev.off()

lst <- list()
for (i in 1:5) {
    r1 <- raster(paste('res/null_model_',sq[i],'.asc',sep=''))
    proj4string(r1) <- albers
    lst[[i]] <- r1
    r2 <- raster(paste('res/forest_model_',sq[i],'.asc',sep=''))
    proj4string(r2) <- albers
    lst[[i+5]] <- r2
}
levelplot(stack(lst), layout=c(5,2))
