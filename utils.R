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
    tupi <- read.csv("sites/tupi_filtered.csv")
    cal <- calibrate(tupi$C14Age, tupi$C14SD, calCurves=tupi$calCurves,
                     resOffsets=tupi$resOffsets, resErrors=tupi$resErrors)
    res <- matrix(ncol=3, nrow=max(sites$id))
    colnames(res) <- c("model", "null", "moist")
    for (i in 1:max(sites$id)) {
        simDates <- sites[sites$id == i,]
        for (j in c("null", "forest")) {
            veg_model <- simDates[simDates$forest %like% j,]
            cat("\nbootstrapping dates...\n")
            pb <- txtProgressBar(min = 0, max = num_iter, style = 3)
            errors <- vector(mode="numeric", length=num_iter)

            for (k in 1:num_iter) {
                sampled <- bootstrapDates(cal)
                rmse <- sqrt(sum((sampled - veg_model$sim_dates)^2) / length(sampled))
                errors[k] <- rmse
                setTxtProgressBar(pb, k)
            }
            res[i, "model"] <- i
            res[i, j] <- paste(mean(errors), sd(errors))
        }
    }
    res.df <- as.data.frame(res)
    write.csv(res.df, "img/res.csv")
    return(res.df)

    if(FALSE) {
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
    return(res.df) }
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
    sites$simBP[is.na(sites$simBP)] <- 0
    cal <- calibrate(sites$C14Age, sites$C14SD, calCurves=sites$calCurves, resOffsets=sites$resOffsets, resErrors=sites$resErrors, verbose=FALSE)
    errors <- vector(mode="numeric", length=num_iter)
    if (verbose) {
        pb <- txtProgressBar(min = 0, max = num_iter, style = 3)
    }
    for (k in 1:num_iter) {
        #sampled <- bootstrapDates(cal)
        #rmse <- sqrt(sum((sampled - sites$simBP)^2) / length(sampled))
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
    sites <- read.csv("sites/tupi_filtered_100.csv")
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
    clusterExport(cl, varlist=c("biomes", "sites", "frontSpeed", "simulateDispersal", "sampleDates", "bootstrapDates"), envir = environment())
    res <- parLapply(cl, params, function(x) {
        a <- x[1]
        delta <- x[2]
        cost <- x[3]
        costSurface <- biomes
        costSurface[values(costSurface) > 1] <- cost
        speed <- frontSpeed(a, delta, 30)
        isochrones <- simulateDispersal(costSurface, c(-61.9574, -10.9557), 5000, speed)
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

plotSim <- function(iso, title) {
    plot(iso, col=viridis(9), breaks=500*1:10, axes=F, box=F, main=title)
    plot(coast, add=T, col="black")
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
write.csv(tupi.df, "sites/tupi_filtered.csv")
}
