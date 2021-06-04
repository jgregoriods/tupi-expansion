library(sp)
library(rcarbon)
library(dplyr)
library(data.table)
library(gdistance)
library(parallel)
library(viridisLite)
library(rnaturalearth)

wgs <- CRS("+init=epsg:4326")
albers <- CRS("+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42
               +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs")

bootstrapDates <- function(calDates) {
    sampledCalBP <- numeric(length(calDates))
	for (i in 1:length(calDates)) {
        calBP <- sample(calDates[i]$grids[[1]]$calBP, size = 1,
                        prob = calDates[i]$grids[[1]]$PrDens)
        sampledCalBP[i] <- calBP
    }
    return(sampledCalBP)
}

simulateDispersal <- function(costRaster, origin, date, speed) {
    tr <- transition(costRaster, function(x) 1 / mean(x), 16)
    tr <- geoCorrection(tr)
    ac <- accCost(tr, origin)
    isochrones <- date - (ac / 1000 / speed)
    isochrones[values(isochrones) < 500] <- NA
    return(isochrones)
}

sampleDates <- function(isochrones, sites, num_iter=100, verbose=TRUE) {
    sites$simBP <- extract(isochrones, sites)
    if ((sum(is.na(sites$simBP)) / nrow(sites)) > 0.1) {return(Inf)}
    sites <- sites[!is.na(sites$simBP),]
    errors <- vector(mode="numeric", length=num_iter)
    if (verbose) {
        pb <- txtProgressBar(min = 0, max = num_iter, style = 3)
    }
    for (k in 1:num_iter) {
        rmse <- sqrt(sum((sites$bp - sites$simBP)^2) / nrow(sites))
        errors[k] <- rmse
        if (verbose) {
            setTxtProgressBar(pb, k)
        }
    }
    if (verbose) {
        close(pb)
    }
    return(mean(errors))
}

frontSpeed <- function(a, delta, T) {
    D <- (delta^2)/(4*T)
    return((2*sqrt(a*D)) / (1 + (a*T/2)))
}

testModels <- function() {
    bestParams <- c()
    bestScore <- Inf
    sites <- read.csv("sites/tupi_filtered_100b.csv")
    biomes <- raster("layers/biomes.asc")
    proj4string(biomes) <- CRS("+init=epsg:4326")
    coordinates(sites) <- ~Xadj+Yadj
    i <- 1
    params <- list()
    for (a in seq(0.02,0.04,0.005)) {
        for (delta in seq(40,60,5)) {
            for (cost in 1:10) {
                params[[i]] <- c(a, delta, cost)
                i <- i + 1
            }
        }
    }
    ncores <- detectCores() - 1
    cl <- makeCluster(ncores)
    clusterEvalQ(cl, library("gdistance"))
    clusterEvalQ(cl, library("rcarbon"))
    clusterExport(cl, varlist=c("biomes", "sites", "frontSpeed", "simulateDispersal", "sampleDates", "bootstrapDates"), envir=environment())
    res <- parLapply(cl, params, function(x) {
        a <- x[1]
        delta <- x[2]
        cost <- x[3]
        costSurface <- biomes
        costSurface[values(costSurface) > 1] <- cost
        speed <- frontSpeed(a, delta, 30)
        isochrones <- simulateDispersal(costSurface, c(-61.96, -10.96), 5000, speed)
        score <- sampleDates(isochrones, sites, verbose=FALSE)
        gc()
        return(c(a, delta, cost, score))
    })
    stopCluster(cl)
    res.df <- as.data.frame(matrix(unlist(res), nrow=length(res), byrow=TRUE))
    colnames(res.df) <- c("a", "delta", "cost", "score")
    res.df <- res.df[order(res.df$score),]
    return(res.df)
}
