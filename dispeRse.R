library(gdistance)

simulateDispersal <- function(costRaster, origin, date, speed) {
    tr <- transition(costRaster, function(x) 1 / mean(x), 16)
    tr <- geoCorrection(tr)
    ac <- accCost(tr, origin)
    isochrones <- date - (ac / 1000 / speed)
    isochrones[values(isochrones) < 500] <- NA
    return(isochrones)
}

sampleDates <- function(isochrones, sites) {
    sites$simBP <- extract(isochrones, sites)
    sites$simBP[is.na(sites$simBP)] <- 0
    rmse <- sqrt(sum((sites$simBP - sites$bp)^2) / length(sites$bp))
    return(rmse)
}

frontSpeed <- function(a, delta, T) {
    D <- (delta^2)/(4*T)
    return((2*sqrt(a*D)) / (1 + (a*T/2)))
}
