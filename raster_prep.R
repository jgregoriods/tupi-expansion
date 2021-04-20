library(raster)
library(rgdal)

albers <- CRS("+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42
               +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs")
wgs <- CRS("+init=epsg:4326")

gmted <- raster("rasters/gmted.tif")
sam <- readOGR("shp/south_america.shp")
ele <- crop(gmted, sam)
ele.m <- projectRaster(ele, crs=albers, res=50000)
writeRaster(ele.m, "layers/ele.asc")

gc()

makeBiomeRasters <- function(type="null") {
    biomes <- vector(mode="list", length=7)
    for (i in 2:7) {
        bio <- raster(paste("rasters/S", i, "biomes_", i - 1, "k.asc", sep=""))
        proj4string(bio) <- wgs
        bio.m <- projectRaster(bio, ele.m, method="ngb")
        if (type == "null") {
            bio.m[values(bio.m) != 1 & values(bio.m) != 3 & values(bio.m) != 7] <- 0
            bio.m[values(bio.m) > 0] <- 1
        } else if (type == "dry") {
            bio.m[values(bio.m) != 1 & values(bio.m) != 7] <- 0
            bio.m[values(bio.m) > 0] <- 1
        } else if (type == "moist") {
            bio.m[values(bio.m) != 1] <- 0
        }
        if (type != "null") {
        eco <- focal(bio.m, w=matrix(1,3,3), fun=mean, na.rm=T)
        bio.m[values(bio.m) == 1 & values(eco) < 1] <- 2
        }
        biomes[[i]] <- bio.m
        gc()
    }
    for (i in 2:7) {
        writeRaster(biomes[[i]], paste("layers/veg_", type, "/veg_", i - 1, "000.asc", sep=""))
    }
}

makeBiomeRasters("null")
makeBiomeRasters("dry")
makeBiomeRasters("moist")
