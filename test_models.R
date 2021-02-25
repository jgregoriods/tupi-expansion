library(parallel)

start_time <- Sys.time() # remove later

r_values <- c(0.03)
fission_values <- c(0.75)
leap_values <- c(0, 150)
start_dates <- c(4100)

param_df <- expand.grid(r_values, fission_values, leap_values, start_dates)
param_list <- split(param_df, seq(nrow(param_df)))

test_model <- function(params, coords, forest_threshold) {
    r <- params[1]
    fission <- params[2]
    leap <- params[3]
    date <- params[4]
    res <- matrix(nrow=10, ncol=5)
    for (i in 1:10) {
        cmd <- paste("./expand", date, coords[1], coords[2], 1, r, fission, leap, forest_threshold)
        a <- system(cmd, intern=TRUE)
        for (j in 1:5) {
            b <- strsplit(a[[j]], " ")
            res[i, j] <- as.numeric(b[[1]][2])
        }
    }
    return(res)
}

encontro_xy <- c(-167889, 2409569)
urupa_xy <- c(-208434, 2438371)

num_cores <- detectCores() - 1
final <- do.call(rbind, mclapply(param_list, test_model, encontro_xy, 0.0, mc.cores=num_cores))
final_veg <- do.call(rbind, mclapply(param_list, test_model, encontro_xy, 0.4, mc.cores=num_cores))

print(Sys.time() - start_time) # remove later