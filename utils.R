library(sp)
library(rcarbon)
library(dplyr)
library(data.table)
library(gdistance)
library(parallel)

# projections
wgs <- CRS("+init=epsg:4326")
albers <- CRS("+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42
               +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs")

#' Models arrival times (in years BP) given a geographical centre of origin
#' and front speed. A cost surface determines acceleration or slowing down of
#' the expansion in different locations.
#'
#' @param cost A cost surface of class Raster*.
#' @param origin A pair of (x, y) coordinates for origin of the expansion.
#' @param date The start date in years BP.
#' @param speed The speed of propagation in km/yr.
#' @return a raster with simulated arrival times.
disperse <- function(cost, origin, date, speed) {
    tr <- transition(cost, function(x) 1 / mean(x), 16)
    tr <- geoCorrection(tr)
    ac <- accCost(tr, origin)
    arrivalTimes <- date - (ac / 1000 / speed)
    arrivalTimes[values(arrivalTimes) < 500] <- NA
    return(arrivalTimes)
}

#' Compares simulated arrival times with radiocarbon dates and returns the root
#' mean square error (RMSE).
#'
#' @param arrivalTimes A raster with simulated arrival times.
#' @param sites A SpatialPoint* object with archaeo sites and a column "bp".
#' @return the root mean square error.
getScore <- function(arrivalTimes, sites) {
    sites$simBP <- extract(arrivalTimes, sites)
    if ((sum(is.na(sites$simBP)) / nrow(sites)) > 0.1) {return(Inf)}
    sites <- sites[!is.na(sites$simBP),]
    rmse <- sqrt(sum((sites$bp - sites$simBP)^2) / nrow(sites))
    return(rmse)
}

#' Calculates the speed of a demic diffusion using the time-lag adjusted
#' Fisher's equation (Fort 2015) <doi:10.1098/rsif.2015.0166>.
#'
#' @param growth The annual growth rate.
#' @param displacement The mean individual displacement in km per generation.
#' @param generation The generation time in years.
#' @return a speed in km/yr.
frontSpeed <- function(growth, displacement, generation) {
    D <- (displacement^2)/(4*generation)
    return((2*sqrt(growth*D)) / (1 + (growth*generation/2)))
}

#' Tests different combinations of parameters.
#' @return a data frame with the scores for the various tested models.
testModels <- function() {
    sites <- read.csv("sites/tupi_dates.csv")
    biomes <- raster("layers/biomes.asc")
    proj4string(biomes) <- CRS("+init=epsg:4326")
    coordinates(sites) <- ~Xadj+Yadj
    growthRates <- seq(0.02,0.04,0.005)
    displacements <- seq(40,60,5)
    costValues <- 1:10
    params.df <- expand.grid(growthRates, displacements, costValues)
    params <- split(params.df, seq(nrow(params.df)))
    #i <- 1
    #params <- list()
    #for (a in seq(0.02,0.04,0.005)) {
    #    for (delta in seq(40,60,5)) {
    #        for (cost in 1:10) {
    #            params[[i]] <- c(a, delta, cost)
    #            i <- i + 1
    #        }
    #    }
    #}
    ncores <- detectCores() - 1
    cl <- makeCluster(ncores)
    clusterEvalQ(cl, library("gdistance"))
    clusterEvalQ(cl, library("rcarbon"))
    clusterExport(cl, varlist=c("biomes", "sites", "frontSpeed",
                                "disperse", "getScore"), envir=environment())
    cat(paste("Running", length(params), "dispersal models. This may take a while ...\n"))
    res <- parLapply(cl, params, function(x) {
        growth <- x[1,1]
        displacement <- x[1,2]
        cost <- x[1,3]
        costSurface <- biomes
        costSurface[values(costSurface) > 1] <- cost
        speed <- frontSpeed(growth, displacement, 30)
        arrivalTimes <- disperse(costSurface, c(-61.96, -10.96), 5000, speed)
        score <- getScore(arrivalTimes, sites)
        gc()
        return(c(growth, displacement, cost, score))
    })
    stopCluster(cl)
    res.df <- as.data.frame(matrix(unlist(res), nrow=length(res), byrow=TRUE))
    colnames(res.df) <- c("growth", "displacement", "cost", "score")
    res.df <- res.df[order(res.df$score),]
    return(res.df)
}
