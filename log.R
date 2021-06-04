library(grid)
library(rasterVis)
library(rgdal)

source("utils.R")

scores <- testModels()
write.csv(scores, "results/dispersal_scores.csv")

coast <- readOGR("shp/south_america.shp")

biomes <- raster("layers/biomes.asc")
proj4string(biomes) <- wgs
speed <- frontSpeed(0.025, 50, 30)
origin <- c(-61.96, -10.96)

cost.null <- biomes
cost.null[values(cost.null) > 1] <- 1
iso.null <- simulateDispersal(cost.null, origin, 5000, speed)

cost.forest <- biomes
cost.forest[values(cost.forest) > 1] <- 4
iso.forest <- simulateDispersal(cost.forest, origin, 5000, speed)

nm <- raster("results/rasters/null_model.asc")
proj4string(nm) <- albers

fm <- raster("results/rasters/forest_model.asc")
proj4string(fm) <- albers

nm <- projectRaster(nm, iso.null)
fm <- projectRaster(fm, iso.null)

# ISOCHRONES FIGURE (DISPERSAL)
plt <- levelplot(stack(iso.null, iso.forest), at=500*1:10, col.regions=viridis(9),
          names.attr=c("a", "b"),
          #scales=list(x=list(draw=FALSE), y=list(draw=F)),
          scales=list(alternating=3),
          xlab="", ylab="",
          colorkey=list(height=0.7, width=1)) + layer(sp.polygons(coast))
plot(plt)

dev.print(jpeg, "img/disperse.jpg", width=1800, height=1200, res=300)
dev.off()

# ISOCHRONES FIGURE (SIM)
plt <- levelplot(stack(nm, fm), at=500*1:10, col.regions=viridis(9),
          names.attr=c("a", "b"),
          #scales=list(x=list(draw=FALSE), y=list(draw=F)),
          scales=list(alternating=3),
          xlab="", ylab="",
          colorkey=list(height=0.7, width=1)) + layer(sp.polygons(coast))
plot(plt)

dev.print(jpeg, "img/sim.jpg", width=1800, height=1200, res=300)
dev.off()

# TIME SLICE FIGURE
sq <- c(seq(4970, 1370, -900), 500)
lst <- list()
for (i in 1:6) {
    r1 <- raster(paste('results/rasters/null_model_',sq[i],'.asc',sep=''))
    proj4string(r1) <- albers
    r1 <- projectRaster(r1, iso.null, method="ngb")
    lst[[i]] <- r1
    r2 <- raster(paste('results/rasters/forest_model_',sq[i],'.asc',sep=''))
    proj4string(r2) <- albers
    r2 <- projectRaster(r2, iso.null, method="ngb")
    lst[[i+6]] <- r2
}

plt <- levelplot(stack(lst), layout=c(6,2), col.regions = gray(seq(1,0,-1)),
                 names.attr=c(paste(as.character(rep(sq, 2)), "BP")),
                 #scales=list(x=list(draw=FALSE), y=list(draw=F)),
                 scales=list(alternating=3),
                 xlab="", ylab="",
                 colorkey=FALSE) +
       layer(sp.polygons(coast))
plot(plt)

dev.print(jpeg, "img/slices.jpg", width=2800, height=1600, res=300)
dev.off()

# SCATTERPLOT
real_dates <- read.csv("sites/tupi_filtered_100b.csv")
sim_dates <- read.csv("results/sim_dates.csv")

sim_dates_null <- sim_dates[sim_dates$model == "null" & sim_dates$sim_dates > 0,]
sim_dates_forest <- sim_dates[sim_dates$model == "forest" & sim_dates$sim_dates > 0,]

par(mfrow=c(1, 2))
plot(real_dates$dist, real_dates$bp, pch=21, bg="white", xlab="distance from origin (km)", ylab="age (cal BP)", main="a")
points(sim_dates_null$dist, sim_dates_null$sim_dates, pch=19)
plot(real_dates$dist, real_dates$bp, pch=21, bg="white", xlab="distance from origin (km)", ylab="age (cal BP)", main="b")
points(sim_dates_forest$dist, sim_dates_forest$sim_dates, pch=19)

dev.print(jpeg, "img/scatterplot.jpg", width=3000, height=1500, res=300)
dev.off()
