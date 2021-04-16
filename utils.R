library(sp)
library(rcarbon)
library(dplyr)

filterDates <- function(sites, radius) {
    clusters <- zerodist(sites, zero=radius, unique.ID=T)
    sites$clusterID <- clusters
    sites.df <- as.data.frame(sites)
    sites.max <- as.data.frame(sites.df %>% group_by(clusterID) %>% top_n(1, bp))
    coordinates(sites.max) <- ~Longitude+Latitude
    proj4string(sites.max) <- proj4string(sites)
    return(sites.max)
}

getDist <- function(sites, origin) {
    dists <- spDistsN1(sites, origin, longlat=T)
    return(dists)
}

tupi <- read.csv("sites/tupi_all.csv")
coordinates(tupi) <- ~Longitude+Latitude
proj4string(tupi) <- CRS("+init=epsg:4326")

cal.dates <- calibrate(tupi$C14Age, tupi$C14SD, calCurves=tupi$calCurves,
                       resOffsets=tupi$resOffsets, resErrors=tupi$resErrors)
tupi$bp <- medCal(cal.dates)

tupi.f <- filterDates(tupi, 0)
tupi.f$dist <- getDist(tupi.f, tupi.f[1,])
tupi.df <- data.frame(x=tupi.f$Longitude, y=tupi.f$Latitude, bp=tupi.f$bp, dist=tupi.f$dist)
write.csv(tupi.df, "tupi_filtered.csv")
