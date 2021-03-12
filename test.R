start_time <- Sys.time() # remove later

# Coordinates of Encontro site
coords <- c(-167889, 2409569)

NUM_ITER <- 100
NUM_DATES <- 6

date <- 4000
k <- 1
r <- 0.04
fission <- 0.75
leap <- c(0, 150)
forest_threshold <- c(0, 1)


res <- matrix(nrow=NUM_ITER, ncol=NUM_DATES)

pb = txtProgressBar(min=0, max=NUM_ITER, initial=0, style=3)
for (i in 1:NUM_ITER) {
    cmd <- paste("./expand", date, coords[1], coords[2], k, r, fission, leap[1], forest_threshold[1])
    a <- system(cmd, intern=TRUE)
    for (j in 1:NUM_DATES) {
        b <- strsplit(a[[j]], " ")
        res[i, j] <- as.numeric(b[[1]][2])
    }
    setTxtProgressBar(pb, i)
}
close(pb)


res_v <- matrix(nrow=NUM_ITER, ncol=NUM_DATES)

pb = txtProgressBar(min=0, max=NUM_ITER, initial=0, style=3)
for (i in 1:NUM_ITER) {
    cmd <- paste("./expand", date, coords[1], coords[2], k, r, fission, leap[2], forest_threshold[2])
    a <- system(cmd, intern=TRUE)
    for (j in 1:NUM_DATES) {
        b <- strsplit(a[[j]], " ")
        res_v[i, j] <- as.numeric(b[[1]][2])
    }
    setTxtProgressBar(pb, i)
}
close(pb)


print(Sys.time() - start_time) # remove later
