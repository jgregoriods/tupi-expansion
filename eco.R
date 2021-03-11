library(raster)

get_ecotone <- function(r) {
  eco <- focal(r, w=matrix(1,3,3), fun=mean)
  eco[values(eco)==1] <- 0
  eco[values(eco)!=0] <- 1
  return(eco)
}

for(i in 1000*1:6) {
  veg <- raster(paste("veg_", i, ".asc", sep=""))
  eco <- get_ecotone(veg)
  writeRaster(eco, paste("eco_", i, ".asc", sep=""))
}
