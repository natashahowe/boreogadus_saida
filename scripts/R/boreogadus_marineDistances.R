# find marine distance

# FOLLOW THIS FUNCTION
#https://jorgemfa.medium.com/how-to-calculate-minimum-marine-distances-in-r-49e897b7de0a

## Source the main function 
source("https://raw.githubusercontent.com/jorgeassis/marineDistances/master/Script.R")
## Read the landmass polygon
global.polygon <- "C:/Users/Natasha.Howe/Work/Mapping/ggOceanMaps/NOAA_GSHHS_Shoreline/GSHHS_shp/h/GSHHS_h_L1.shp"

## Run the function
contour(global.polygon = global.polygon, file= "./data/R/PopLatLong_noChukchiCoronation.txt", 
        file.sep = "\t", file.dec = ".", file.strucutre = 2, file.header = FALSE,
        resolution = 0.01, buffer = c(4,4,1,4), export.file = TRUE)

## file : the main file with the locations; should be text delimited
## global.polygon: the path of the polygon
## file.strucutre: the main file structure: 1 to “Name Lon Lat” or 2 to “Name Lat Lon”
## file.header: define if the text file has a header with the column names (TRUE or FALSE)
## resolution: the resolution of the study area and the buffer to use around the sites. 
## buffer: the buffer can be a simple value or a vector such as c(xmin,xmax,ymin,ymax). 
## export.file: file to export the results as a text delimited file (TRUE or FALSE)