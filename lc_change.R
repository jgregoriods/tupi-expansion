library(motif)
library(stars)
library(tmap)

f <- NULL
for (i in 1:6) {
    a <- read_stars(paste("rasters/S", i, "biomes_", i - 1, "k.asc", sep=""))
    b <- read_stars(paste("rasters/S", i + 1, "biomes_", i, "k.asc", sep=""))
    c <- lsp_compare(a, b, window=2, type="cove", dist_fun="jensen-shannon")
    if (is.null(f)) {
        f <- c
    } else {
        f <- f + c
    }
}
f <- f / 6
r <- as(f, "Raster")