#install.packages("ggOceanMaps")
library(ggOceanMaps); library(sf); library(stars); 
library(ggnewscale); library(here); library(tidyverse)

.ggOceanMapsenv <- new.env()
.ggOceanMapsenv$datapath <- 'C:/Users/Natasha.Howe/Work/Mapping/ggOceanMaps'

etopoPath <- "C:/Users/Natasha.Howe/Work/Mapping/ggOceanMaps"
outPath <- "./figures/map"

dt <- data.frame(lon = c(-170, -170, -35, -35), 
                 lat = c(29, 70, 70, 29))

popDF <- read.csv(paste0(here(),"./data/R/PopLatLong.csv"), header = T, row.names = NULL)
color_df <- read.delim2(paste0(here(),"./data/R/color_metadata_allpops_flip.txt"), 
                        header = T, row.names = NULL, sep = "\t")

colorPops <- popDF %>%
  left_join(color_df, by = "Pop")

mypalette <- colorPops$Color2
names(mypalette) <- colorPops$Population
colorPops$Population <- factor(colorPops$Population, levels = unique(colorPops$Population))

map1 <- basemap(data = dt, bathymetry = TRUE) + 
  geom_polygon(data = transform_coord(dt), aes(x = lon, y = lat), 
               color = NA, fill = NA) +
  ggnewscale::new_scale_fill() +
  geom_point(data = colorPops, aes(x = Longitude, y = Latitude, fill = Population), color = "black",
             size = 4, pch = 21) +
  scale_fill_manual(values = mypalette)
ggsave(paste0("./figures/map/arctic_cod_map_",format(Sys.Date(), format = "%Y%m%d"),".jpeg"), map1, width = 10, height = 6)

arcticCRS <- sf::st_crs(shapefile_list("Arctic")$crs)

# Binned raster bathymetry
dd_rbathy <- raster_bathymetry(
  bathy = file.path(etopoPath, "ETOPO_2022_v1_60s_N90W180_surface.nc"),
  depths = c(50, 300, 500, 1000, 1500, 2000, 4000, 6000, 10000), 
  proj.out = 4326,
  downsample = 1
)

save(dd_rbathy, file = file.path(outPath, "ggOceanMapsData/dd_rbathy.rda"), compress = "bzip2")

# Vector bathymetry
dd_bathy <- vector_bathymetry(dd_rbathy, drop.crumbs = 25, smooth = TRUE) 
save(dd_bathy, file = file.path(outPath, "ggOceanMapsData/dd_bathy.rda"), compress = "xz")

# Continuous raster bathymetry
dd_rbathy_cont <- raster_bathymetry(
  bathy = file.path(etopoPath, "ETOPO_2022_v1_60s_N90W180_surface.nc"),
  depths = NULL, 
  proj.out = 4326,
  downsample = 1
)

save(dd_rbathy_cont, file = file.path(outPath, "ggOceanMapsData/dd_rbathy_cont.rda"), compress = "bzip2")

########## ADD ICELAND ##########################

dt <- data.frame(lon = c(-170, -170, -10, -10), 
                 lat = c(29, 70, 70, 29))

popDF <- read.csv(paste0(here(),"./data/R/PopLatLong.csv"), header = T, row.names = NULL)
color_df <- read.delim2(paste0(here(),"./data/R/color_metadata_allpops_flip.txt"), 
                        header = T, row.names = NULL, sep = "\t")

colorPops <- popDF %>%
  left_join(color_df, by = "Pop")

mypalette <- colorPops$Color2
names(mypalette) <- colorPops$Population
colorPops$Population <- factor(colorPops$Population, levels = unique(colorPops$Population))
mypalette

map2 <- basemap(data = dt, bathymetry = TRUE) + 
  geom_polygon(data = transform_coord(dt), aes(x = lon, y = lat), 
               color = NA, fill = NA) +
  ggnewscale::new_scale_fill() +
  geom_point(data = colorPops, aes(x = Longitude, y = Latitude, fill = Population), color = "black",
             size = 4, pch = 21) +
  scale_fill_manual(values = mypalette)

ggsave(paste0("./figures/map/boreogadus_map_Iceland_",format(Sys.Date(), format = "%Y%m%d"),".jpeg"), 
       map2, width = 10, height = 6, dpi = 300)

