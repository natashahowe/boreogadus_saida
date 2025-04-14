# find marine distances

# FOLLOW THIS FUNCTION
#https://jorgemfa.medium.com/how-to-calculate-minimum-marine-distances-in-r-49e897b7de0a

## Source the main function 
source("https://raw.githubusercontent.com/jorgeassis/marineDistances/master/Script.R")

# resolution of GSHHS file. can be c-l, with c being the least res and l being the most
res="h"

## Read the landmass polygon
global.polygon <- paste0("../Mapping/GSHHS_shp/",res,"/GSHHS_",res,"_L1.shp")

## Run the function
contour(global.polygon = global.polygon, file= "./data/R/PopLatLong_noChukchiCoronation.txt", 
        file.sep = "\t", file.dec = ".", file.strucutre = 2, file.header = FALSE,
        resolution = 0.01, buffer = c(4,6,1,4), export.file = TRUE)

    ## file : the main file with the locations; should be text delimited
    ## global.polygon: the path of the polygon
    ## file.strucutre: the main file structure: 1 to “Name Lon Lat” or 2 to “Name Lat Lon”
    ## file.header: define if the text file has a header with the column names (TRUE or FALSE)
    ## resolution: the resolution of the study area and the buffer to use around the sites. 
    ## buffer: the buffer can be a simple value or a vector such as c(xmin,xmax,ymin,ymax). 
    ## export.file: file to export the results as a text delimited file (TRUE or FALSE)

# move pdf
file.rename("Contour _ Study Region.pdf", 
            paste0("./figures/ibd/marineDistances/boreogadus_StudyRegion_noChukchiCoronation_res-",res,".pdf"))

# after it runs, read in the dataframe and rewrite it out
mat.marineDist <- as.matrix(read.delim2("Contour _ Pairwise Marine Distances.txt",sep="\t",
                           row.names = 1, header = T))
  # store values as.numeric
  mat.marineDist <- apply(mat.marineDist, 2, as.numeric)
  row.names(mat.marineDist) <- colnames(mat.marineDist)
  write.csv(mat.marineDist,
            paste0("./results/ibd/boreogadus_noChukchiCoronation_marinedistances_res-",res,".csv"),
            row.names = T, quote = F)

  mat.geoDist <- as.matrix(read.csv("./results/ibd/geodistances_matrix_20250224.csv",
                                    row.names = 1, header = T))
  simp.geoDist <- mat.geoDist[3:nrow(mat.geoDist),3:nrow(mat.geoDist)]
  simp.marineDist <- mat.marineDist[-5,-5]
  colnames(simp.geoDist)
  colnames(simp.marineDist)
  
  ## random check 
  mat.diff <- simp.geoDist - simp.marineDist[rownames(simp.geoDist), row.names(simp.geoDist)]
  View(mat.diff)

  
  
######## TRY TO RUN IT WITH CHUKCHI/CORONATION ##########################
  ## Run the function
  contour(global.polygon = global.polygon, file= "./data/R/PopLatLong.txt", 
          file.sep = "\t", file.dec = ".", file.strucutre = 2, file.header = FALSE,
          resolution = 0.01, buffer = c(4,22,1,12), export.file = TRUE)

  # move pdf
  file.rename("Contour _ Study Region.pdf", 
              paste0("./figures/ibd/marineDistances/boreogadus_StudyRegion_allPops_res-",res,".pdf"))
  
  # after it runs, read in the dataframe and rewrite it out
  mat.marineDist <- as.matrix(read.delim2("Contour _ Pairwise Marine Distances.txt",sep="\t",
                                          row.names = 1, header = T))
  # store values as.numeric
  row.names(mat.marineDist)[row.names(mat.marineDist) == "SE Baffin Island"] <- "SEBaffin"
    row.names(mat.marineDist) <- gsub(" .*","",row.names(mat.marineDist))
    colnames(mat.marineDist) <- row.names(mat.marineDist)
  mat.marineDist <- apply(mat.marineDist, 2, as.numeric)
    row.names(mat.marineDist) <- colnames(mat.marineDist)
  
  write.csv(mat.marineDist,
            paste0("./results/ibd/boreogadus_allPops_marinedistances_res-",res,"2.csv"),
            row.names = T, quote = F)
  
  
  mat.geoDist <- as.matrix(read.csv("./results/ibd/geodistances_matrix_20250224.csv",
                                    row.names = 1, header = T))
  simp.marineDist <- mat.marineDist[-5,-5]
  colnames(simp.geoDist)
  colnames(simp.marineDist)
  
  simp.marineDist <- simp.marineDist[rownames(mat.geoDist), row.names(mat.geoDist)]
  write.csv(simp.marineDist,
            paste0("./results/ibd/boreogadus_noLabrador_marinedistances_res-",res,"2.csv"),
            row.names = T, quote = F)
  
  ## random check 
  mat.diff <- mat.geoDist - simp.marineDist
  View(mat.diff)
  