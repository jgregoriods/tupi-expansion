library(motif)
library(stars)
library(tmap)
library(ggplot2)
library(raster)
library(rnaturalearth)

coast <- ne_coastline(scale=50, returnclass="sf")
coast <- st_crop(coast, xmin=-82, xmax=-33, ymin=-60, ymax=12)

f <- NULL
for (i in 1:6) {
    a <- read_stars(paste("rasters/S", i, "biomes_", i - 1, "k.asc", sep=""))
    b <- read_stars(paste("rasters/S", i + 1, "biomes_", i, "k.asc", sep=""))
    c <- lsp_compare(a, b, window=10, type="cove", dist_fun="jensen-shannon")
    if (is.null(f)) {
        f <- c
    } else {
        f <- f + c
    }
}
f <- f / 6
r <- as(f, "Raster")

g <- ggplot() +
        geom_raster(data=as.data.frame(r[[4]], xy=T), aes(x=x, y=y, fill=dist)) +
        geom_sf(data=coast, fill=NA, color="black", size=0.5) +
        coord_sf() +
        scale_fill_viridis_c(na.value='#FFFFFF00', name='Distance (JSD)') +
        theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
            panel.background = element_rect(fill = "white", colour = NA),
            panel.border=element_blank(), axis.ticks=element_blank(),
            axis.text.x=element_blank(), axis.text.y=element_blank(),
            legend.position=c(0.8, 0.2))
plot(g)

dev.print(pdf, "img/lc_change.pdf")
dev.off()
