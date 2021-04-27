library(ggplot2)
library(raster)
library(tidyr)
library(sf)
library(gridExtra)
library(rnaturalearth)

bio <- read_sf("shp/forests.shp")
bio <- st_crop(bio, xmin=-85, xmax=-33, ymin=-60, ymax=12)
world <- ne_countries(scale=50, returnclass="sf")
coast <- ne_coastline(scale=50, returnclass="sf")
coast <- st_crop(coast, xmin=-85, xmax=-33, ymin=-60, ymax=12)
sam <- st_crop(world, xmin=-85, xmax=-33, ymin=-60, ymax=12)
tp <- read_sf("shp/tupi_languages.shp")

sites <- read.csv("sites/tupi_all.csv")
sites$size <- 1

plt_a <- ggplot() +
    geom_sf(data=sam, fill="white", color=NA) +
    geom_sf(data=bio, aes(fill="darkolivegreen3"), show.legend=T, color=NA) +
    geom_sf(data=coast, fill=NA, color="black", size=0.2) +
    geom_point(data=sites, aes(x=Longitude, y=Latitude, color="black"), size=0.25) +
    scale_fill_identity(labels=c("Tropical moist forests"), guide = "legend") +
    scale_color_identity(labels=c("Dated sites"), guide = "legend") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
          panel.background = element_rect(fill = "aliceblue"),
          legend.position="bottom", legend.title = element_blank())

plt_b <- ggplot() +
    geom_sf(data=sam, fill="white", color=NA) +
    geom_sf(data=tp, aes(fill="coral"), color=NA) +
    geom_sf(data=coast, fill=NA, color="black", size=0.2) +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
          panel.background = element_rect(fill = "aliceblue"),
          legend.position="bottom", legend.title = element_blank()) +
    scale_fill_identity(labels=c("Tupi languages"), guide = "legend")

grid.arrange(plt_a, plt_b, ncol = 2)
#dev.print(jpeg, "img/tupi.jpeg", res=300, width=2048)
dev.print(pdf, "img/tupi.pdf")
dev.off()
