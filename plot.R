library(rcarbon)
library(sp)

tp <- read.csv("tupi.csv")
cal <- calibrate(tp$C14Age, tp$C14SD)

coordinates(tp) <- ~Longitude+Latitude
proj4string(tp) <- "+init=epsg:4326"

tp$dist <- spDistsN1(tp, tp[1,], longlat=TRUE)

for (i in 1:nrow(tp)) {
    cal$grid[[i]]$PrDens <- cal$grid[[i]]$PrDens / max(cal$grid[[i]]$PrDens) * 50 + tp$dist[i]
}

plot(0, xlim=c(4000, 2000), ylim=c(1200, 2800))

dists <- c(1579, 1611, 2459, 1523, 2623)

for (i in 1:ncol(final)) {
    d <- density(final[,i])
    d$y <- d$y / max(d$y) * 50 + dists[i]
    polygon(d, col="yellow")
}

for (i in 1:ncol(final_veg)) {
    final_veg[final_veg == 0] <- NA
    d <- density(final_veg[,i], na.rm=TRUE)
    d$y <- d$y / max(d$y) * 50 + dists[i]
    polygon(d, col="green")
}

for (i in 1:nrow(tp)) {
    polygon(cal$grid[[i]]$calBP, cal$grid[[i]]$PrDens, col=rgb(0.5,0.5,0.5,0.25), border=rgb(0,0,0,0))
}
