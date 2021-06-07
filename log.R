library(grid)
library(rasterVis)
library(rgdal)
library(viridisLite)

source("utils.R")

coastline <- readOGR("shp/south_america.shp")
biomes <- raster("layers/biomes.asc")
proj4string(biomes) <- wgs

scores <- testModels()
write.csv(scores, "results/dispersal_scores.csv")

speed <- frontSpeed(0.025, 50, 30)
origin <- c(-61.96, -10.96)

cost.null <- biomes
cost.null[values(cost.null) > 1] <- 1
iso.null <- disperse(cost.null, origin, 5000, speed)

cost.forest <- biomes
cost.forest[values(cost.forest) > 1] <- 4
iso.forest <- disperse(cost.forest, origin, 5000, speed)

simulation.null <- raster("results/rasters/null_model.asc")
proj4string(simulation.null) <- albers

simulation.forest <- raster("results/rasters/forest_model.asc")
proj4string(simulation.forest) <- albers

simulation.null <- projectRaster(simulation.null, iso.null)
simulation.forest <- projectRaster(simulation.forest, iso.null)

# Figure 2
plt <- levelplot(stack(iso.null, iso.forest), at=500*1:10,
                 col.regions=viridis(9), names.attr=c("a", "b"),
                 scales=list(alternating=3), xlab="", ylab="",
                 colorkey=list(height=0.7, width=1)) +
       layer(sp.polygons(coastline))
plot(plt)

dev.print(jpeg, "img/disperse.jpg", width=1800, height=1200, res=300)
dev.off()

# Figure 3
plt <- levelplot(stack(simulation.null, simulation.forest), at=500*1:10,
                 col.regions=viridis(9), names.attr=c("a", "b"),
                 scales=list(alternating=3), xlab="", ylab="",
                 colorkey=list(height=0.7, width=1)) +
       layer(sp.polygons(coastline))
plot(plt)

dev.print(jpeg, "img/sim.jpg", width=1800, height=1200, res=300)
dev.off()

# Figure 4
real_dates <- read.csv("sites/tupi_dates.csv")
sim_dates <- read.csv("results/sim_dates.csv")

sim_dates_null <- sim_dates[sim_dates$model == "null" & sim_dates$sim_dates > 0,]
sim_dates_forest <- sim_dates[sim_dates$model == "forest" & sim_dates$sim_dates > 0,]

coordinates(real_dates) <- ~Longitude+Latitude
proj4string(real_dates) <- wgs

dates.null <- extract(iso.null, real_dates)
dates.forest <- extract(iso.forest, real_dates)

par(mfrow=c(1, 2))
plot(real_dates$dist, real_dates$bp, pch=1,
     xlab="distance from origin (km)", ylab="age (cal BP)", main="a")
points(real_dates$dist, dates.null, pch=3, col="red")
plot(real_dates$dist, real_dates$bp, pch=1,
     xlab="distance from origin (km)", ylab="age (cal BP)", main="b")
points(real_dates$dist, dates.forest, pch=3, col="red")

dev.print(jpeg, "img/scatterplot1.jpg", width=3000, height=1500, res=300)
dev.off()

par(mfrow=c(1, 2))
plot(real_dates$dist, real_dates$bp, pch=1,
     xlab="distance from origin (km)", ylab="age (cal BP)", main="a")
points(sim_dates_null$dist, sim_dates_null$sim_dates, pch=3, col="red")
plot(real_dates$dist, real_dates$bp, pch=1,
     xlab="distance from origin (km)", ylab="age (cal BP)", main="b")
points(sim_dates_forest$dist, sim_dates_forest$sim_dates, pch=3, col="red")

dev.print(jpeg, "img/scatterplot2.jpg", width=3000, height=1500, res=300)
dev.off()

# Figure 5
years <- c(seq(4970, 1370, -900), 500)
timeSlices <- list()
for (i in 1:6) {
    r.null <- raster(paste('results/rasters/null_model_', years[i], '.asc', sep=''))
    proj4string(r.null) <- albers
    r.null <- projectRaster(r.null, iso.null, method="ngb")
    timeSlices[[i]] <- r.null
    r.forest <- raster(paste('results/rasters/forest_model_', years[i], '.asc', sep=''))
    proj4string(r.forest) <- albers
    r.forest <- projectRaster(r.forest, iso.null, method="ngb")
    timeSlices[[i+6]] <- r.forest
}

plt <- levelplot(stack(timeSlices), layout=c(6,2), col.regions = gray(seq(1,0,-1)),
                 names.attr=c(paste(as.character(rep(years, 2)), "BP")),
                 scales=list(alternating=3), xlab="", ylab="", colorkey=FALSE) +
       layer(sp.polygons(coastline))
plot(plt)

dev.print(jpeg, "img/slices.jpg", width=2800, height=1600, res=300)
dev.off()
