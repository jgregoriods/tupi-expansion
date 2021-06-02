library(raster)
library(viridisLite)

res <- raster("rasters/S1biomes_0k.asc")
proj4string(res) <- CRS("+init=epsg:4326")
r2 <- raster(extent(res))
proj4string(r2) <- CRS("+init=epsg:4326")
res(r2) <- 0.25
res <- projectRaster(res, r2)
res[!is.na(values(res))] <- -1
for (i in seq(6000, 0, -1000)) {
    a <- raster(paste("rasters/S", (i / 1000) + 1, "biomes_", i / 1000, "k.asc", sep=""))
    proj4string(a) <- proj4string(res)
    a <- projectRaster(a, res, method="ngb")
    if (i == 6000) {
        res[values(a) == 1] <- i
    } else {
        a[values(a) > 1] <- 0
        dif <- a - (res > -1)
        res[values(dif) == 1] <- i
    }
}
res[values(res) == -1] <- NA
plot(res, col=plasma(6, direction=-1), breaks=seq(6000,0,-1000))
