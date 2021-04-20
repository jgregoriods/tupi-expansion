library(rcarbon)
library(sp)

tp <- read.csv("tupi.csv")
cal <- calibrate(tp$C14Age, tp$C14SD, calCurves="shcal20")

for (i in 1:nrow(tp)) {
    cal$grid[[i]]$PrDens <- cal$grid[[i]]$PrDens / max(cal$grid[[i]]$PrDens) * 0.8 + i
}

plot(0, xlim=c(3200, 1000), ylim=c(0.5, nrow(tp) + 1))

for (i in 1:nrow(tp)) {
    polygon(cal$grid[[i]]$calBP, cal$grid[[i]]$PrDens, col=rgb(0,0,0))
}

for (i in 1:nrow(tp)) {
    d <- density(res[,i])
    d_v <- density(res_v[,i])

    d$y <- d$y / max(d$y) * 0.8 + i
    d_v$y <- d_v$y / max(d_v$y) * 0.8 + i

    polygon(d, col=rgb(1, 1, 0))
    polygon(d_v, col=rgb(0, 1, 0))
}
