library(raster)

get_ecotone <- function(r) {
  eco <- focal(r, w=matrix(1,3,3), fun=mean)
  eco[values(eco)==1] <- 0
  eco[values(eco)!=0] <- 1
  return(eco)
}

run_eco_files <- function(){
  for(i in 1000*1:6) {
    veg <- raster(paste("veg_", i, ".asc", sep=""))
    eco <- get_ecotone(veg)
    writeRaster(eco, paste("eco_", i, ".asc", sep=""))
  }
}

eco_forest_cells <- function() {
  for(i in 1000*1:6) {
    veg <- raster(paste("veg_", i, ".asc", sep=""))
    eco <- raster(paste("eco_", i, ".asc", sep=""))
    veg[values(veg)==1 & values(eco)==1] <- 2
    writeRaster(veg, paste("veg2_", i, ".asc", sep=""))
  }
}