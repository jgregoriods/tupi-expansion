library(raster)
library(viridisLite)

res <- raster("rasters/S1biomes_0k.asc")
res[!is.na(values(res))] <- -1
for (i in seq(6000, 0, -1000)) {
    a <- raster(paste("rasters/S", (i / 1000) + 1, "biomes_", i / 1000, "k.asc", sep=""))
    if (i == 6000) {
        res[values(a) == 1] <- i
    } else {
        a[values(a) > 1] <- 0
        dif <- a - (res > -1)
        res[values(dif) == 1] <- i
    }
}
res[values(res) == -1] <- NA
plot(res, col=viridis(7))
