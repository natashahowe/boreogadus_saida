#install.packages("ggOceanMaps")
library(ggOceanMaps); library(sf); library(stars); library(ggnewscale)

.ggOceanMapsenv <- new.env()
.ggOceanMapsenv$datapath <- 'C:/Users/Natasha.Howe/Work/Mapping/ggOceanMaps'

etopoPath <- "C:/Users/Natasha.Howe/Work/Mapping/ggOceanMaps"
NEDPath <- ""
gebcoPath <- ""
outPath <- "C:/Users/Natasha.Howe/Work/Mapping/test_outputs"
EEAPath <- ""


dt <- data.frame(lon = c(-170, -170, -35, -35), 
                 lat = c(29, 70, 70, 29))

popDF <- data.frame(lon = c(-160, -55, -113), 
                    lat = c(71, 53, 68),
                    Pop = c("Chukchi", "Labrador", "Coronation Gulf"))

testmap1 <- basemap(data = dt, bathymetry = TRUE) + 
  geom_polygon(data = transform_coord(dt), aes(x = lon, y = lat), 
               color = NA, fill = NA) +
  ggnewscale::new_scale_fill() +
  geom_point(data = popDF, aes(x = lon, y = lat, fill = Pop), color = "black",
             size = 4, pch = 21) 
ggsave(paste0(outPath,"/testmap1_10162024.jpeg"), testmap1, width = 10, height = 6)

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