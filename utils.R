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

bootstrapDates <- function(calDates) {
    sampledCalBP <- numeric(length(calDates))
	for (i in 1:length(calDates)) {
        calBP <- sample(calDates[i]$grids[[1]]$calBP, size = 1,
                        prob = calDates[i]$grids[[1]]$PrDens)
        sampledCalBP[i] <- calBP
    }
    return(sampledCalBP)
}

getScore <- function(filename, num_iter=100) {
    sites <- read.csv(filename)
    simDates <- sites %>% select(contains('sim'))
    cal <- calibrate(sites$C14Age, sites$C14SD, calCurves=sites$calCurves,
                     resOffsets=sites$resOffsets, resErrors=sites$resErrors)
    res <- matrix(ncol=2, nrow=ncol(simDates))
    colnames(res) <- c("mean", "sd")
    for (i in 1:ncol(simDates)) {
        cat("bootstrapping dates...\n")
        cat(paste("simulation", i, "of", ncol(simDates), "\n"))
        pb <- txtProgressBar(min = 0, max = num_iter, style = 3)
        errors <- vector(mode="numeric", length=num_iter)
        for (j in 1:num_iter) {
            sampled <- bootstrapDates(cal)
            rmse <- sqrt(sum((sampled - simDates[,i])^2) / length(sampled))
            errors[j] <- rmse
            setTxtProgressBar(pb, j)
        }
        close(pb)
        res[i,1] <- mean(errors)
        res[i,2] <- sd(errors)
    }
    res.df <- as.data.frame(res)
    write.csv(res.df, "img/res.csv")
    return(res.df)
}

if(FALSE){
tupi <- read.csv("sites/tupi_all.csv")
coordinates(tupi) <- ~Longitude+Latitude
proj4string(tupi) <- CRS("+init=epsg:4326")

cal.dates <- calibrate(tupi$C14Age, tupi$C14SD, calCurves=tupi$calCurves,
                       resOffsets=tupi$resOffsets, resErrors=tupi$resErrors)
tupi$bp <- medCal(cal.dates)

tupi.f <- filterDates(tupi, 0)
tupi.f$dist <- getDist(tupi.f, tupi.f[1,])
#tupi.df <- data.frame(x=tupi.f$Longitude, y=tupi.f$Latitude, bp=tupi.f$bp, dist=tupi.f$dist)
tupi.df <- as.data.frame(tupi.f)[c("Longitude", "Latitude", "dist", "bp", "C14Age", "C14SD", "calCurves", "resOffsets", "resErrors")]
write.csv(tupi.df, "tupi_filtered_2.csv")
}
