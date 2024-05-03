# Gastropod SDM and Elaphostrongylus DD model
# Authors: Hannah Rose Vineer and Anna Ciezarek, University of Liverpool
# Contact: Hannah.Vineer@liverpool.ac.uk
# Last modified: 25/07/2022

# NOTE: It is our intention that each section can be run independently.
# You may therefore find a lot of repetition in this code.
# Being able to run each session independently is convenient for time management
# (because some sections are very slow), but also to allow you to restart R
# occassionally, which we have found frees up memory. 

# First, make space in your R environment for the large spatial datasets
require(usethis) 
usethis::edit_r_environ()
# In the new tab that opens, enter: R_MAX_VSIZE=100Gb (or whatever is appropriate for your computer), save, then restrat R Studio

# Packages required
require(raster)
require(rnaturalearth)
require(usdm)
require(maptools)
require(spatstat)
require(rJava)
require(dismo)
require(grid)
require(rgbif)
require(CoordinateCleaner)
require(dplyr)
require(countrycode)
require(beeswarm)
require(rgdal)
require(viridis)
require(tidyr)
require(ggplot2)


# Useful sources:
# https://mltconsecol.github.io/TU_LandscapeAnalysis_Documents/Assignments_web/Assignment08_SpeciesDistributionModeling_Pt2.pdf page 7
# http://modata.ceoe.udel.edu/dev/dhaulsee/class_rcode/r_pkgmanuals/maxentusingdismo.pdf 
# https://rdrr.io/cran/dismo/man/mess.html 
# https://groups.google.com/g/maxent/c/yRBlvZ1_9rQ 
# https://datasquad.at.sites.carleton.edu/data/storage-design/dealing-with-a-vector-memory-exhausted-error-in-r/
# https://cran.r-project.org/web/packages/CoordinateCleaner/vignettes/Cleaning_GBIF_data_with_CoordinateCleaner.html
# 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### 1. Raster preprocessing ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Specify the directory where you have saved your raster layers and where you want to save your Bioclim data (downloaded in next step)
setwd("ENTER THE DIRECTORY PATH HERE")
setwd("/Users/hannahvineer/OneDrive - The University of Liverpool/GIS_data/Bioclim/Bioclim_30s/")

# 1.1 Cropping ---------------------------------------------------------------

# Load your rasters and rename layers
elev = raster("wc2.1_30s_elev.tif")
bio.hist = stack(paste0(getwd(), "/wc2.1_30s_bio_historic/", list.files(path = "wc2.1_30s_bio_historic/")), elev)
names(bio.hist)
names(bio.hist) = c("bio1", "bio10", "bio11", "bio12", "bio13", "bio14", "bio15", 
                    "bio16", "bio17", "bio18", "bio19", "bio2", "bio3", "bio4",
                    "bio5", "bio6", "bio7", "bio8", "bio9", "elev")

bio.CN.ssp245.2030s = stack("wc2.1_30s_bioc_CNRM-ESM2-1_ssp245_2021-2040.tif", elev)
names(bio.CN.ssp245.2030s)
names(bio.CN.ssp245.2030s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")
bio.CN.ssp370.2030s = stack("wc2.1_30s_bioc_CNRM-ESM2-1_ssp370_2021-2040.tif", elev)
names(bio.CN.ssp370.2030s)
names(bio.CN.ssp370.2030s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                             "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                             "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.CN.ssp585.2030s = stack("wc2.1_30s_bioc_CNRM-ESM2-1_ssp585_2021-2040.tif", elev)
names(bio.CN.ssp585.2030s)
names(bio.CN.ssp585.2030s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.IP.ssp585.2030s = stack("wc2.1_30s_bioc_IPSL-CM6A-LR_ssp585_2021-2040.tif", elev)
names(bio.IP.ssp585.2030s)
names(bio.IP.ssp585.2030s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.CN.ssp245.2050s = stack("wc2.1_30s_bioc_CNRM-ESM2-1_ssp245_2041-2060.tif", elev)
names(bio.CN.ssp245.2050s)
names(bio.CN.ssp245.2050s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.CN.ssp370.2050s = stack("wc2.1_30s_bioc_CNRM-ESM2-1_ssp370_2041-2060.tif", elev)
names(bio.CN.ssp370.2050s)
names(bio.CN.ssp370.2050s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.CN.ssp585.2050s = stack("wc2.1_30s_bioc_CNRM-ESM2-1_ssp585_2041-2060.tif", elev)
names(bio.CN.ssp585.2050s)
names(bio.CN.ssp585.2050s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")


bio.IP.ssp585.2050s = stack("wc2.1_30s_bioc_IPSL-CM6A-LR_ssp585_2041-2060.tif", elev)
names(bio.IP.ssp585.2050s)
names(bio.IP.ssp585.2050s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.CN.ssp245.2070s = stack("wc2.1_30s_bioc_CNRM-ESM2-1_ssp245_2061-2080.tif", elev)
names(bio.CN.ssp245.2070s)
names(bio.CN.ssp245.2070s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.CN.ssp370.2070s = stack("wc2.1_30s_bioc_CNRM-ESM2-1_ssp370_2061-2080.tif", elev)
names(bio.CN.ssp370.2070s)
names(bio.CN.ssp370.2070s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.CN.ssp585.2070s = stack("wc2.1_30s_bioc_CNRM-ESM2-1_ssp585_2061-2080.tif", elev)
names(bio.CN.ssp585.2070s)
names(bio.CN.ssp585.2070s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.IP.ssp585.2070s = stack("wc2.1_30s_bioc_IPSL-CM6A-LR_ssp585_2061-2080.tif", elev)
names(bio.IP.ssp585.2070s)
names(bio.IP.ssp585.2070s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.CN.ssp245.2090s = stack("wc2.1_30s_bioc_CNRM-ESM2-1_ssp245_2081-2100.tif", elev)
names(bio.CN.ssp245.2090s)
names(bio.CN.ssp245.2090s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.CN.ssp370.2090s = stack("wc2.1_30s_bioc_CNRM-ESM2-1_ssp370_2081-2100.tif", elev)
names(bio.CN.ssp370.2090s)
names(bio.CN.ssp370.2090s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.CN.ssp585.2090s = stack("wc2.1_30s_bioc_CNRM-ESM2-1_ssp585_2081-2100.tif", elev)
names(bio.CN.ssp585.2090s)
names(bio.CN.ssp585.2090s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

bio.IP.ssp585.2090s = stack("wc2.1_30s_bioc_IPSL-CM6A-LR_ssp585_2081-2100.tif", elev)
names(bio.IP.ssp585.2090s)
names(bio.IP.ssp585.2090s) = c("bio1", "bio2", "bio3", "bio4", "bio5", "bio6", "bio7", 
                               "bio8", "bio9", "bio10", "bio11", "bio12", "bio13", "bio14",
                               "bio15", "bio16", "bio17", "bio18", "bio19", "elev")

# 1.1.1 Crop historic data =====================================================
# Get Europe boundary shapefile
europe = ne_countries(scale = "medium", type = "countries", continent = c("Europe", "Asia"))
europe = crop(europe, c(-25, 45, 35, 85))
plot(europe)

# Crop the historic data first
bio.hist = crop(bio.hist, extent(europe), progress = "text") # historic 1970-2000
plot(bio.hist$bio1)
# 1.1.2 Check and select environmental variables ===============================

# Need to check environmental variables to remove any highly correlated vars and
# remove any with a variance inflation factor of more than 10

# Stepwise removal of variables with VIF>10. Save as an object for subsetting
vif.selection = vifstep(bio.hist, th = 10)

# Use the object created in the previous line to subset the stack so that only 
# variables with no collinearity issues remain
bio.hist = exclude(bio.hist, vif.selection)
bio.CN.ssp245.2030s = exclude(bio.CN.ssp245.2030s, vif.selection)
bio.CN.ssp370.2030s = exclude(bio.CN.ssp370.2030s, vif.selection)
bio.CN.ssp585.2030s = exclude(bio.CN.ssp585.2030s, vif.selection)
bio.IP.ssp585.2030s = exclude(bio.IP.ssp585.2030s, vif.selection)
bio.CN.ssp245.2050s = exclude(bio.CN.ssp245.2050s, vif.selection)
bio.CN.ssp370.2050s = exclude(bio.CN.ssp370.2050s, vif.selection)
bio.CN.ssp585.2050s = exclude(bio.CN.ssp585.2050s, vif.selection)
bio.IP.ssp585.2050s = exclude(bio.IP.ssp585.2050s, vif.selection)
bio.CN.ssp245.2070s = exclude(bio.CN.ssp245.2070s, vif.selection)
bio.CN.ssp370.2070s = exclude(bio.CN.ssp370.2070s, vif.selection)
bio.CN.ssp585.2070s = exclude(bio.CN.ssp585.2070s, vif.selection)
bio.IP.ssp585.2070s = exclude(bio.IP.ssp585.2070s, vif.selection)
bio.CN.ssp245.2090s = exclude(bio.CN.ssp245.2090s, vif.selection)
bio.CN.ssp370.2090s = exclude(bio.CN.ssp370.2090s, vif.selection)
bio.CN.ssp585.2090s = exclude(bio.CN.ssp585.2090s, vif.selection)
bio.IP.ssp585.2090s = exclude(bio.IP.ssp585.2090s, vif.selection)

# 1.1.3 Crop all other raster stacks ================================

bio.CN.ssp245.2030s = crop(bio.CN.ssp245.2030s, extent(europe), progress = "text")
bio.CN.ssp370.2030s = crop(bio.CN.ssp370.2030s, extent(europe), progress = "text")
bio.CN.ssp585.2030s = crop(bio.CN.ssp585.2030s, extent(europe), progress = "text")
bio.IP.ssp585.2030s = crop(bio.IP.ssp585.2030s, extent(europe), progress = "text")
bio.CN.ssp245.2050s = crop(bio.CN.ssp245.2050s, extent(europe), progress = "text")
bio.CN.ssp370.2050s = crop(bio.CN.ssp370.2050s, extent(europe), progress = "text")
bio.CN.ssp585.2050s = crop(bio.CN.ssp585.2050s, extent(europe), progress = "text")
bio.IP.ssp585.2050s = crop(bio.IP.ssp585.2050s, extent(europe), progress = "text")
bio.CN.ssp245.2070s = crop(bio.CN.ssp245.2070s, extent(europe), progress = "text")
bio.CN.ssp370.2070s = crop(bio.CN.ssp370.2070s, extent(europe), progress = "text")
bio.CN.ssp585.2070s = crop(bio.CN.ssp585.2070s, extent(europe), progress = "text")
bio.IP.ssp585.2070s = crop(bio.IP.ssp585.2070s, extent(europe), progress = "text")
bio.CN.ssp245.2090s = crop(bio.CN.ssp245.2090s, extent(europe), progress = "text")
bio.CN.ssp370.2090s = crop(bio.CN.ssp370.2090s, extent(europe), progress = "text")
bio.CN.ssp585.2090s = crop(bio.CN.ssp585.2090s, extent(europe), progress = "text")
bio.IP.ssp585.2090s = crop(bio.IP.ssp585.2090s, extent(europe), progress = "text")

# 1.2 Masking ------------------------------------------------------------------
bio.hist = mask(bio.hist, europe, progress = "text") # historic 1970-2000
bio.CN.ssp245.2030s = mask(bio.CN.ssp245.2030s, europe, progress = "text")
bio.CN.ssp370.2030s = mask(bio.CN.ssp370.2030s, europe, progress = "text")
bio.CN.ssp585.2030s = mask(bio.CN.ssp585.2030s, europe, progress = "text")
bio.IP.ssp585.2030s = mask(bio.IP.ssp585.2030s, europe, progress = "text")
bio.CN.ssp245.2050s = mask(bio.CN.ssp245.2050s, europe, progress = "text")
bio.CN.ssp370.2050s = mask(bio.CN.ssp370.2050s, europe, progress = "text")
bio.CN.ssp585.2050s = mask(bio.CN.ssp585.2050s, europe, progress = "text")
bio.IP.ssp585.2050s = mask(bio.IP.ssp585.2050s, europe, progress = "text")
bio.CN.ssp245.2070s = mask(bio.CN.ssp245.2070s, europe, progress = "text")
bio.CN.ssp370.2070s = mask(bio.CN.ssp370.2070s, europe, progress = "text")
bio.CN.ssp585.2070s = mask(bio.CN.ssp585.2070s, europe, progress = "text")
bio.IP.ssp585.2070s = mask(bio.IP.ssp585.2070s, europe, progress = "text")
bio.CN.ssp245.2090s = mask(bio.CN.ssp245.2090s, europe, progress = "text")
bio.CN.ssp370.2090s = mask(bio.CN.ssp370.2090s, europe, progress = "text")
bio.CN.ssp585.2090s = mask(bio.CN.ssp585.2090s, europe, progress = "text")
bio.IP.ssp585.2090s = mask(bio.IP.ssp585.2090s, europe, progress = "text")

# Save rasters
setwd("ENTER THE DIRECTORY WHERE YOU WANT TO SAVE RASTERS HERE")
setwd("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Bioclim_Europe/")
lapply(X = ls(pattern = "bio."), FUN = function (x) writeRaster(get(x), filename = paste0(x, ".tif"), format="GTiff"))

for (i in ls(pattern="bio.")){
  tiff(filename = paste0("preview_",i,".tif"), width = 8, height = 8, units = "in", res = 300)
  par(mfrow = c(3,3))
  plot(get(i))
  dev.off()
}

# Clean your environment to avoid conflict and errors in the next section
rm(list = ls())
gc()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### 2. Occurrence data preprocessing and pseudoabsence selection ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Set working directory for sourcing input and saving occurrence data
setwd("ENTER YOUR DIRECTORY PATH HERE")
setwd("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Occurrence_data/")

# 2.1 Prepare occurrence data (dependent variable) for Maxent modelling --------

# 2.1.1 Get occurrence data ====================================================

gastropodspecies <- c("Arianta arbustorum", "Punctum pygmaeum", 
                      "Vitrina pellucida", "Euconulus fulvus",
                      "Deroceras laeve", "Arion subfuscus")
gastropod.occurrence = occ_search(scientificName = gastropodspecies, hasCoordinate = T, 
                 limit = 50000, decimalLongitude = "-25, 45",
                 decimalLatitude = "30, 85")

gastropod.copy = gastropod.occurrence # save a copy in case you overwrite in error below

gbif.citation = gbif_citation(gastropod.occurrence) # sources to cite in publications

# Check the data
names(gastropod.occurrence)
names(gastropod.occurrence[[gastropodspecies[1]]])
names(gastropod.occurrence[[gastropodspecies[1]]]$meta)
names(gastropod.occurrence[[gastropodspecies[1]]]$data)

# Extract just the data
Aa.occurrence = gastropod.occurrence$`Arianta arbustorum`$data
Pp.occurrence = gastropod.occurrence$`Punctum pygmaeum`$data
Vp.occurrence = gastropod.occurrence$`Vitrina pellucida`$data
Ef.occurrence = gastropod.occurrence$`Euconulus fulvus`$data
Dl.occurrence = gastropod.occurrence$`Deroceras laeve`$data
As.occurrence = gastropod.occurrence$`Arion subfuscus`$data

write.csv(x = dplyr::bind_rows(Aa.occurrence,Pp.occurrence, Vp.occurrence,Ef.occurrence,Dl.occurrence,As.occurrence),
            file = "GBIF_all_records.csv")

gastropod.occurrence = dplyr::bind_rows(Aa.occurrence,Pp.occurrence, Vp.occurrence,Ef.occurrence,Dl.occurrence,As.occurrence)

# Remove incomplete records
gastropod.occurrence = gastropod.occurrence%>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

# View on map
plot(europe)
plot(gastropod.occurrence$decimalLongitude, gastropod.occurrence$decimalLatitude)

# Use coordinate cleaner to remove problematic records
#convert country code from ISO2c to ISO3c
gastropod.occurrence$countryCode <-  countrycode::countrycode(gastropod.occurrence$countryCode, origin =  'iso2c', destination = 'iso3c')
#flag problems
gastropod.occurrence <- data.frame(gastropod.occurrence)
flags <- clean_coordinates(x = gastropod.occurrence, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                     "zeros", "countries", "seas", "validity")) # most test are on by default

summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

gastropods.cleaned = gastropod.occurrence[flags$.summary,]
unique(gastropods.cleaned$species)
# Remove records before 1970 to match the historic environmental data
gastropods.cleaned = gastropods.cleaned[which(gastropods.cleaned$year>=1970),]

# Remove any absences
gastropods.cleaned = gastropods.cleaned[which(gastropods.cleaned$occurrenceStatus=="PRESENT"),]
unique(gastropods.cleaned$species)

# Remove records with uncertainty <900m (approximately the resolution of the environmental data)
gastropods.cleaned = gastropods.cleaned[which(gastropods.cleaned$coordinateUncertaintyInMeters<=900),]
unique(gastropods.cleaned$species)

# Extract only human observations
gastropods.cleaned = gastropods.cleaned[which(gastropods.cleaned$basisOfRecord == "HUMAN_OBSERVATION"),]
unique(gastropods.cleaned$species)
table(gastropods.cleaned$species)

write.csv(x = gastropods.cleaned, file = "GBIF_Gastropods_occurrence_CLEANED.csv")


# 2.1.2 Spatial thinning =======================================================

# Create a new folder in the "Occurrence_data" folder
# This will fail if the directory already exists
# If this is the case, please go into the existing directory and make sure that
# there are no files that risk being overwritten in the loop below
dir.create("Occurrence_subsets/")
dir.create("Occurrence_subsets/thinned_final/")
dir.create("Occurrence_subsets/thinned_final/KDE/")
dir.create("Occurrence_subsets/thinned_final/SWD/")

# Thin the occurrence data to remove spatial clustering
# WARNING: THIS IS SLOW
# this has to be done for each species pres/abs subset separately
# definitions for loop
distance = 0.9 #thinning distances in km, to reflect the resolution of the environmental data

table(gastropods.cleaned$species)

Aa.pres_thinned =
  spThin::thin( loc.data = gastropods.cleaned[which(gastropods.cleaned$species=="Arianta arbustorum"),],
                lat.col = "decimalLatitude", long.col = "decimalLongitude",
                spec.col = "species",
                thin.par = as.numeric(distance), #the distance (in kilometers) that you want records to be separated by.
                reps = 1,
                locs.thinned.list.return = FALSE,
                write.files = TRUE,
                max.files = 1,
                out.dir = paste0(getwd(), "/Occurrence_subsets/thinned_final/"),
                out.base = paste0("pres_thinned_Aarbustorum"), #make sure nothing else called this or it has a funny five minutes and puts _new at the end
                verbose = TRUE)

As.pres_thinned =
  spThin::thin( loc.data = gastropods.cleaned[which(gastropods.cleaned$species=="Arion subfuscus"),],
                lat.col = "decimalLatitude", long.col = "decimalLongitude",
                spec.col = "species",
                thin.par = as.numeric(distance), #the distance (in kilometers) that you want records to be separated by.
                reps = 1,
                locs.thinned.list.return = FALSE,
                write.files = TRUE,
                max.files = 1,
                out.dir = paste0(getwd(), "/Occurrence_subsets/thinned_final/"),
                out.base = paste0("pres_thinned_Asubfuscus"), #make sure nothing else called this or it has a funny five minutes and puts _new at the end
                verbose = TRUE)

Dl.pres_thinned =
  spThin::thin( loc.data = gastropods.cleaned[which(gastropods.cleaned$species=="Deroceras laeve"),],
                lat.col = "decimalLatitude", long.col = "decimalLongitude",
                spec.col = "species",
                thin.par = as.numeric(distance), #the distance (in kilometers) that you want records to be separated by.
                reps = 1,
                locs.thinned.list.return = FALSE,
                write.files = TRUE,
                max.files = 1,
                out.dir = paste0(getwd(), "/Occurrence_subsets/thinned_final/"),
                out.base = paste0("pres_thinned_Dlaeve"), #make sure nothing else called this or it has a funny five minutes and puts _new at the end
                verbose = TRUE)

Ef.pres_thinned =
  spThin::thin( loc.data = gastropods.cleaned[which(gastropods.cleaned$species=="Euconulus fulvus"),],
                lat.col = "decimalLatitude", long.col = "decimalLongitude",
                spec.col = "species",
                thin.par = as.numeric(distance), #the distance (in kilometers) that you want records to be separated by.
                reps = 1,
                locs.thinned.list.return = FALSE,
                write.files = TRUE,
                max.files = 1,
                out.dir = paste0(getwd(), "/Occurrence_subsets/thinned_final/"),
                out.base = paste0("pres_thinned_Efulvus"), #make sure nothing else called this or it has a funny five minutes and puts _new at the end
                verbose = TRUE)

Pp.pres_thinned =
  spThin::thin( loc.data = gastropods.cleaned[which(gastropods.cleaned$species=="Punctum pygmaeum"),],
                lat.col = "decimalLatitude", long.col = "decimalLongitude",
                spec.col = "species",
                thin.par = as.numeric(distance), #the distance (in kilometers) that you want records to be separated by.
                reps = 1,
                locs.thinned.list.return = FALSE,
                write.files = TRUE,
                max.files = 1,
                out.dir = paste0(getwd(), "/Occurrence_subsets/thinned_final/"),
                out.base = paste0("pres_thinned_Ppygmaeum"), #make sure nothing else called this or it has a funny five minutes and puts _new at the end
                verbose = TRUE)

Vp.pres_thinned =
  spThin::thin( loc.data = gastropods.cleaned[which(gastropods.cleaned$species=="Vitrina pellucida"),],
                lat.col = "decimalLatitude", long.col = "decimalLongitude",
                spec.col = "species",
                thin.par = as.numeric(distance), #the distance (in kilometers) that you want records to be separated by.
                reps = 1,
                locs.thinned.list.return = FALSE,
                write.files = TRUE,
                max.files = 1,
                out.dir = paste0(getwd(), "/Occurrence_subsets/thinned_final/"),
                out.base = paste0("pres_thinned_Vpellucida"), #make sure nothing else called this or it has a funny five minutes and puts _new at the end
                verbose = TRUE)

# 2.2 Generate bias raster and select pseudoabsences ==========================

setwd("ENTER PATH TO OCCURRENCE DATA FOLDER HERE")
setwd("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Occurrence_data/")
GIS_dir = ("ENTER THE DIRECTORY WHERE YOU HAVE SAVED GIS LAYER HERE")
GIS_dir = ("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Bioclim_Europe/")

# First convert the boundary shapefile to owin format so that we can make a
# spatial point pattern dataset (ppp)
europe = ne_countries(scale = "medium", type = "countries", continent = "Europe")
europe = crop(europe, c(-25, 45, 30, 85))
proj4string(europe) <- ""
europe = as.owin.SpatialPolygons(europe)
plot(europe)

# Read in the thinned data
Aa = read.csv("Occurrence_subsets/thinned_final/pres_thinned_Aarbustorum_thin1.csv")
As = read.csv("Occurrence_subsets/thinned_final/pres_thinned_Asubfuscus_thin1.csv")
Dl = read.csv("Occurrence_subsets/thinned_final/pres_thinned_Dlaeve_thin1.csv")
Ef = read.csv("Occurrence_subsets/thinned_final/pres_thinned_Efulvus_thin1.csv")
Pp = read.csv("Occurrence_subsets/thinned_final/pres_thinned_Ppygmaeum_thin1.csv")
Vp = read.csv("Occurrence_subsets/thinned_final/pres_thinned_Vpellucida_thin1.csv")

# Combine all datasets
all_thinned = dplyr::bind_rows(Aa,As,Dl,Ef,Pp,Vp)

# Use the new europe owin and the background data (all records) to generate 
# a point pattern dataset
all.ppp = ppp(x = all_thinned$decimalLongitude, y = all_thinned$decimalLatitude, window=europe, check=F)
plot(all.ppp)

# Generate a KDE (kernel density estimate) based on the entire database of occurrence records
# set a severe bandwidth adjustment to account for the fact that study sites
# can be highly localised
bias = stats::density(all.ppp, adjust = 0.5)
# Normalise to values between 0 and 1 because we want to use these values as
# weights later, & convert to raster
normalize <- function(x)
{
  return((x- min(x)) /(max(x)-min(x)))
}
bias = raster(normalize(bias))
plot(bias)
plot(all.ppp, pch = 21, cex = 0.1, add = T) # Check that the bias layer is a good representation of the input data

# Resample so that it's the same resolution as the environmental data
# Get one of your environmental layers
bio.hist = stack(paste0(GIS_dir, "bio.hist.tif"))
bias = resample(bias,bio.hist$bio.hist.1, method = "ngb")
# Write to GTiff in case needed later
writeRaster(x = bias, filename = paste0(getwd(), "/Occurrence_subsets/thinned_final/KDE/KDE_RASTER.tif"), format = "GTiff")

# Convert the raster grid to points and select 10000 of these at random without replacement
# These are the background points that will be used for modelling
potential.bg = rasterToPoints(bias, spatial=F)
BG = potential.bg[sample(seq(1:nrow(potential.bg)), size = 10000, replace = F, prob = potential.bg[,"layer"]), 1:2]
points(BG, pch = 16, cex = 0.1) # check that the selection was correctly weighted using the bias layer

# 2.3 Combine presence & pseudoabsence data ------------------------------------

list.species = c("Aa", "As", "Dl", "Ef", "Pp", "Vp") # create a list of the species occurrence datasets. This should be the variable names used in section 2.2 above

for (i in list.species) {
dat = get(i) # get the dataset for species i
myResp = c(rep(1,length(dat[,1])), rep(0, 10000)) # Create a vector of 0 and 1 to represent presences and background points
colnames(BG) = colnames(dat[,2:3]) # change the column names in the background data to match the occurrence
myRespXY = rbind(dat[,2:3],BG)  # combine occurrence and background data
write.csv(cbind(myResp, myRespXY), file = paste0(getwd(), "/Occurrence_subsets/thinned_final/", i, "_pres_bg.csv")) # save as .csv
}
# Adapted from Guisan et al., 2017. Habitat suitability and distribution models with applications in R. Cambridge University Press, Cambridge, UK
# Chapter 14 (but with reference to other chapters too)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### 3. Extract environmental data for the presence & pseudoabsence points ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Set working directory for sourcing input and saving output
setwd("ENTER PATH TO OCCURRENCE DATA FOLDER HERE")
setwd("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Occurrence_data/")
GIS_dir = ("ENTER THE DIRECTORY WHERE YOU HAVE SAVED GIS LAYER HERE")
GIS_dir = setwd("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Bioclim_Europe/")

# Load your rasters
bio.hist = stack(paste0(GIS_dir, "/bio.hist.tif"))
par(mfrow = c(3,3))
plot(bio.hist) # compare this against the "preview" tiff generated above to be sure that the layers are correct
# Rename each layer to avoid long names
names(bio.hist) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

# Read in the presence-background data for each species
list.species = c("Aa", "As", "Dl", "Ef", "Pp", "Vp") # create a list of the species occurrence datasets. This should be the variable names used in section 2.2 above

for (i in list.species) {
  assign(i, read.csv(paste0("Occurrence_subsets/thinned_final/", i, "_pres_bg.csv")))
}

# For each species, extract the environmental data for model development
for (i in list.species) {
  dat = get(i)
  env = raster::extract(bio.hist, dat[,3:4])
  dat = cbind(dat,env)
  colnames(dat)[1] = i
  dat = na.omit(dat) # Remove any rows with missing data
  write.csv(dat, file = paste0(getwd(), "/Occurrence_subsets/thinned_final/SWD/", i, "_SWD.csv"))
}

# Clear the environment ready for SDM
rm(list = ls())
gc()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### 4. Maxent species distribution modelling ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Set working directory for sourcing input
setwd("ENTER PATH TO OCCURRENCE DATA FOLDER HERE")
setwd("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Occurrence_data/")

# Set the GIS directory where all of your environmental layers are held (the rasters produced in section 1 above)
GIS_dir = ("ENTER THE DIRECTORY WHERE YOU HAVE SAVED GIS LAYER HERE")
GIS_dir = ("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Bioclim_Europe/")

# Set the directory for saving the model output
mod_out_dir = ("ENTER YOUR DIRECTORY PATH HERE")
mod_out_dir = ("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Maxent_output_model_development/")

# 4.1 rJava troubleshooting for Mac users ----------------------------------------
# If you are using a Mac machine you may have trouble with rJava...
# If rJava package won't load, try the following solution from https://stackoverflow.com/questions/44081227/trouble-installing-and-loading-rjava-on-mac-el-capitan:
# In Terminal:
# > cd /Library/Frameworks/R.framework/Versions/3.4/Resources/lib
# > rm libjvm.dylib
# > ln -s /Library/Java/JavaVirtualMachines/jdk1.8.0_151.jdk/Contents/Home/jre/lib/server/libjvm.dylib libjvm.dylib
# now go to R or RStudio and try loading rJava again
# If you still can't run Maxent using dismo because of Java issues, try this 
# ... solution from https://community.rstudio.com/t/rstudio-crashing-with-tabulizer-need-to-install-the-legacy-java-se-6-runtime/87937/2 
# Sys.setenv(JAVA_HOME='/Library/Java/JavaVirtualMachines/jdk-11.0.14.jdk/Contents/Home/')
# detach("package:rJava", unload = TRUE)
# require(rJava)
# See also: https://github.com/wallaceEcoMod/wallace/blob/master/README.md 
# See also: https://stackoverflow.com/questions/30738974/rjava-load-error-in-rstudio-r-after-upgrading-to-osx-yosemite 

# If the above fails...
# Try uninstalling all Java and JDK from your Mac
# Download https://support.apple.com/kb/dl1572?locale=en_US
# Install the above by running the following in Terminal:
# > t=${TMPDIR:-/tmp}/java
# > hdiutil mount /path/to/javaforosx.dmg
# > pkgutil --expand /Volumes/Java\ for\ macOS\ 2017-001/JavaForOSX.pkg "$t"
# > hdiutil unmount /Volumes/Java\ for\ macOS\ 2017-001
# > sed -i '' 's/return false/return true/g' "$t"/Distribution
# > pkgutil --flatten "$t" ~/Desktop/Java.pkg
# > rm -rf "$t"
# > open ~/Desktop/Java.pkg
# Above commands taken from: https://apple.stackexchange.com/questions/375973/java-uninstalled-but-still-cannot-install-java-6-macos
# Download https://www.oracle.com/java/technologies/downloads/#java8 
# Install the above manually
# In Terminal, run the following:
# > sudo R CMD javareconf

# If you get a prefs root node error for Java, try this fix: https://github.com/julienvollering/MIAmaxent/issues/1#issuecomment-278985033

# 4.2 Save Maxent software on your computer --------------------------------------
# Download from https://biodiversityinformatics.amnh.org/open_source/maxent/ 
# Then paste the contents of the zip folder into the "java" folder under dismo 
# ... in your R library
# You can find out where this folder is using...
.Library

# Save the path to the jar file for later
jar = paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

# 4.3 Import environmental data ------------------------------------------------
# Load your rasters
bio.hist = stack(paste0(GIS_dir, "bio.hist.tif"))
# Rename each layer to avoid long names
names(bio.hist) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

# Get Europe boundary
europe = ne_countries(scale = "medium", type = "countries", continent = c("Europe", "Asia"))
europe = crop(europe, c(-25, 45, 35, 85))
plot(europe)


# 4.4. Run Maxent model --------------------------------------------------------

# test whether your system is correctly set up to run Maxent
# if so, a message stating "This is MaxEnt version x.x.x" will appear
# If you receive any errors, check you have downloaded MaxEnt, saved the .jar 
# file in the correct place (described in 5.1), and 
# resolved any rJava and Java issues.
maxent()

# 4.4.1 Maxent replicates ======================================================

spp = c("Aa", "As", "Dl", "Ef", "Pp", "Vp")
for (j in 1:length(spp)) {
  # Read in the presence-background dataset with appended environmental data
  DataSpecies = read.csv(paste0(getwd(), "/Occurrence_subsets/thinned_final/SWD/", spp[j],"_SWD.csv"))

  nV <- 10 # Number of validation replicates
  prop_test = 0.25 # Proportion of data to randomly set aside for validation in each replicate
  nRow <- nrow(DataSpecies)
  
  # Create an array to store the predicted probabilities for each replicate
  Pred_results <- array(0, c(nRow, 1, nV), dimnames=list(seq(1:nRow), "MAXENT", seq(1:nV)))
  
  # Create an array to store the model evaluation
  Eval_results = array(0, c(9, 1, nV), dimnames=list(c("n presences", "n absences", "AUC", "correlation coefficient", "P value correlation", "max TPR+TNR at", "TSS", "TPR", "TNR"), "MAXENT", seq(1:nV)))
  
  # Loop through the validation runs: 
  for(i in 1:nV){ # for each validation run...
    #separate the original data in one sub set for calibration and the other for evaluation. 
    a <- biomod2::SampleMat2(ref=DataSpecies[,3], ratio=1-prop_test) # function from the biomod2 package
    calib <- DataSpecies[a$calibration,] # specifies which data to use for calibration (training)...
    eval <- DataSpecies[a$evaluation,] # ...and for evaluation (testing)
    
    # Run Maxent
    # WARNING: YOU MUST MAKE SURE THE INDICES FOR THE CALIB AND EVAL OBJECTS 
    # MATCH THE INDICES FOR THE ENVIRONMENTAL DATA COLUMNS
    if (file.exists(jar) & require(rJava)) {
      maxent_mod = dismo::maxent(x = calib[, 6:14], p = calib[,3], path = paste0(mod_out_dir, spp[j], "_v_", i), args = c("responsecurves=TRUE", "pictures=TRUE", "jackknife=TRUE"))
      Pred_test =  dismo::predict(maxent_mod, eval[, 6:14], type="response")
      Pred_results[,"MAXENT",i] = dismo::predict(maxent_mod, DataSpecies, type="prob")
      
      # run evaluation to get threshold for mapping presence/absence map and
      # evaluation data
      # WARNING: YOU MUST MAKE SURE THE INDICES REPRESENT THE PRES/ABS COLUMN 
      # AND THE X AND Y COORDS
      e <- dismo::evaluate(model = maxent_mod, p = eval[which(eval[,3]==1),4:5], a = eval[which(eval[,3]==0),4:5], x = bio.hist)
      Eval_results[,"MAXENT",i] = rbind(e@np, e@na, e@auc, e@cor, e@pcor, e@t[which.max(e@TPR + e@TNR)], max(e@TPR+e@TNR)-1, e@TPR[which.max(e@TPR + e@TNR)], e@TNR[which.max(e@TPR + e@TNR)])

      # Predict probability of presence based on entire dataset
      px = dismo::predict(bio.hist, maxent_mod, progress='', filename= paste0(mod_out_dir, spp[j], "_v_", i,"/prediction"))
      writeRaster(px, filename = paste0(mod_out_dir, spp[j], "_v_", i,"/prediction_map.tif"), format = "GTiff") # Save the model prediction map
      threshold <- e@t[which.max(e@TPR + e@TNR)] # Find the threshold that maximises true positive rate and true negative rate
      binary_px = px>=threshold # generate a binary map where values above this threshold = 1
      writeRaster(binary_px, filename = paste0(mod_out_dir, spp[j], "_v_", i,"/prediction_map_BINARY.tif"), format = "GTiff") # Save binary maps
      
      tiff(paste0(mod_out_dir, spp[j], "_v_", i,"/Prediction_boxplot.tif"), width = 6, height = 4, units = "in", res = 300)
      boxplot(e, names = c("Pseudoabsence", "Presence"), ylab = "Predicted probability", xlab = "") # Plot the predicted probabilities for occurrence and pseudoabsence points
      dev.off()
      
    } else {
      print("cannot run this example because MaxEnt is not available")
    }
    print(paste("Finished MAXENT", spp[j], "v", i, sep = " ")) # Shows progress (because this loop can be slow)
  }
  # Save training and test results
  write.csv(Pred_results, file = paste0(mod_out_dir, spp[j], "_Validation_Pred_Results.csv"))
  write.csv(Eval_results, file = paste0(mod_out_dir, spp[j], "_Validation_Eval_Results.csv"))  
  }

# Note that the training AUC for each validation replicate is within the folder for
# that run (in the maxent results csv file). The test AUC is in the evaluation 
# csv file that contains the evaluation results for all of the validation
# runs.

# 4.5 Analyse replicate model output -------------------------------------------

# Set the directory for saving the model output
mod_out_dir = ("ENTER YOUR DIRECTORY PATH HERE")
mod_out_dir = ("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Maxent_output_model_development/")

# Each of the folders generated in the previous step contains a map of model predictions

# Get evaluation stats
# Transposing the data frame puts one variable per column (instead of one per row)
# but generates a matrix, which must be converted to a data frame to be able to use $

list.species = c("Aa", "As", "Dl", "Ef", "Pp", "Vp") # create a list of the species occurrence datasets. This should be the variable names used in section 2.2 above

for (i in list.species) {
  assign(paste0(i, ".eval"), as.data.frame(t(read.csv(paste0(mod_out_dir, i, "_Validation_Eval_Results.csv"), header = T, row.names = 1))))
  }

# Get training stats - these are all in separate folders for the cross-validations
# List all directories for the species
dirs <- list.dirs(path = mod_out_dir, recursive = F)

for (i in list.species){
  assign(paste0(i,".dirs"), grep(i, dirs, value = T))
  assign(paste0(i, ".train"), read.csv(paste0(get(paste0(i, ".dirs"))[1], "/maxentResults.csv")))
  for (k in 2:length(get(paste0(i, ".dirs")))) {
    assign(paste0(i, ".train"), rbind(get(paste0(i, ".train")), read.csv(paste0(get(paste0(i, ".dirs"))[k], "/maxentResults.csv"))))
  }
  write.csv(get(paste0(i, ".train")), paste0(mod_out_dir, i, "_Validation_Train_Results.csv"))
  }

# Summarise

tiff(paste0(mod_out_dir, "AUC_boxplots.tif"), width = 6, height = 4, units = "in", res = 300)
AUC_comp = data.frame(matrix(ncol = length(list.species), nrow = 10))
colnames(AUC_comp) = list.species
for (i in 1:length(list.species)) {
  AUC_comp[,i] = get(paste0(list.species[i], ".eval"))$AUC
  }
boxplot(AUC_comp, las = 1, names = list.species, ylab = "AUC", outline = F, col = "gray95")
beeswarm(AUC_comp, add = TRUE, pch = 16, cex = 0.8)
dev.off()

tiff(paste0(mod_out_dir, "TSS_boxplots.tif"), width = 6, height = 4, units = "in", res = 300)
TSS_comp = data.frame(matrix(ncol = length(list.species), nrow = 10))
colnames(TSS_comp) = list.species
for (i in 1:length(list.species)) {
  TSS_comp[,i] = get(paste0(list.species[i], ".eval"))$TSS
}
boxplot(TSS_comp, las = 1, names = list.species, ylab = "TSS", outline = F, col = "gray95")
beeswarm(TSS_comp, add = TRUE, pch = 16, cex = 0.8)
dev.off()

for (i in list.species) {
tiff(paste0(mod_out_dir, i, "_contributions.tif"), width = 6, height = 6, units = "in", res = 300)
boxplot(get(paste0(i, ".train"))[,c(11:15, 8:10, 16)], outline = F, col = "gray95", las = 2, names = c("bio2", "bio3", "bio4", "bio8", "bio9", "bio15", "bio18", "bio19", "elevation"), ylab = "Contribution to model performance (%)")
beeswarm(get(paste0(i, ".train"))[,c(11:15, 8:10, 16)], add = TRUE, pch = 16, cex = 0.8)
dev.off()
}

for (i in list.species) {
  tiff(paste0(mod_out_dir, i, "_permutation_importance.tif"), width = 6, height = 6, units = "in", res = 300)
  boxplot(get(paste0(i, ".train"))[,c(20:24, 17:19, 25)], outline = F, col = "gray95", las = 2, names = c("bio2", "bio3", "bio4", "bio8", "bio9", "bio15", "bio18", "bio19", "elevation"), ylab = "Permutation importance (%)")
  beeswarm(get(paste0(i, ".train"))[,c(20:24, 17:19, 25)], add = TRUE, pch = 16, cex = 0.8)
  dev.off()
}


# Read the maps in to R
for (i in list.species) {
  for (k in 1:length(get(paste0(i, ".dirs")))) {
    assign(paste0("prob.pres.", i, k), raster(paste0(get(paste0(i, ".dirs"))[k], "/prediction_map.tif")))
    assign(paste0("binary.", i, k), raster(paste0(get(paste0(i, ".dirs"))[k], "/prediction_map_BINARY.tif")))
    }

  # Create a stack
  # I couldn't find a way to use assign and ls to stack the rasters, so another loop is necessary!
  list.rasters = ls(pattern = paste0("prob.pres.", i))
  assign(paste0(i, ".prob.stk"), stack(get(list.rasters[1]), get(list.rasters[2])))
  for (ras in 3:length(list.rasters)) {
    assign(paste0(i, ".prob.stk"), stack(get(paste0(i, ".prob.stk")), get(list.rasters[ras])))
    }
  
  list.rasters = ls(pattern = paste0("binary.", i))
  assign(paste0(i, ".binary.stk"), stack(get(list.rasters[1]), get(list.rasters[2])))
  for (ras in 3:length(list.rasters)) {
    assign(paste0(i, ".binary.stk"), stack(get(paste0(i, ".binary.stk")), get(list.rasters[ras])))
    }

  # Calculate mean, standard deviation and model agreement  
  assign(paste0(i, ".mean"), calc(get(paste0(i, ".prob.stk")), mean))
  assign(paste0(i, ".sd"), calc(get(paste0(i, ".prob.stk")), sd))
  assign(paste0(i, ".agreement"), calc(get(paste0(i, ".binary.stk")), sum))

  writeRaster(get(paste0(i, ".mean")), filename = paste0(mod_out_dir, i, "_mean_probability_of_presence.tif"), format = "GTiff")
  writeRaster(get(paste0(i, ".sd")), filename = paste0(mod_out_dir, i, "_sd_probability_of_presence.tif"), format = "GTiff")
  writeRaster(get(paste0(i, ".agreement")), filename = paste0(mod_out_dir, i, "i_binary_model_agreement.tif"), format = "GTiff")
}


# Read sumary maps in (if you are starting from this point in a new R session)
for (i in list.species) {
  assign(paste0(i, ".mean"), raster(paste0(mod_out_dir, i, "_mean_probability_of_presence.tif")))
  assign(paste0(i, ".sd"), raster(paste0(mod_out_dir, i, "_sd_probability_of_presence.tif")))
  assign(paste0(i, ".agreement"), raster(paste0(mod_out_dir, i, "i_binary_model_agreement.tif")))
  }

tiff(paste0(mod_out_dir,"Summary_maps.tif"), width = 2*length(list.species), height = 6, units = "in", res = 300)
  par(mfrow = c(3,length(list.species)), oma = c(1,1,1,1))
  plot(get(paste0(list.species[1],".mean")), ylab = "Latitude", main = list.species[1], xaxt = "n")
  mtext(text = "Mean", side = 3, line = 0, cex = 0.8)
  for (i in list.species[2:length(list.species)]) {
    plot(get(paste0(i,".mean")), main = i, xaxt = "n")
    mtext(text = "Mean", side = 3, line = 0, cex = 0.8)
    }
  plot(get(paste0(list.species[1],".sd")), ylab = "Latitude", xaxt = "n")
  mtext(text = "S.D.", side = 3, line = 0, cex = 0.8)
  for (i in list.species[2:length(list.species)]) {
    plot(get(paste0(i,".sd")), xaxt = "n")
    mtext(text = "S.D.", side = 3, line = 0, cex = 0.8)
    }
  plot(get(paste0(list.species[1],".agreement")), xlab = "Longitude", ylab = "Latitude")
  mtext(text = "Agreement", side = 3, line = 0, cex = 0.8)
  for (i in list.species[2:length(list.species)]) {
    plot(get(paste0(i,".agreement")), xlab = "Longitude")
    mtext(text = "Agreement", side = 3, line = 0, cex = 0.8)
    }
  dev.off()
  

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### 5. Projections onto future climate ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Set working directory for sourcing input
setwd("ENTER PATH TO OCCURRENCE DATA FOLDER HERE")
setwd("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Occurrence_data/")

# Set the GIS directory where all of your environmental layers are held (the rasters produced in section 1 above)
GIS_dir = ("ENTER THE DIRECTORY WHERE YOU HAVE SAVED GIS LAYER HERE")
GIS_dir = ("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Bioclim_Europe/")

# Set the directory for saving the model output
mod_out_dir = ("ENTER YOUR DIRECTORY PATH HERE")
mod_out_dir = ("/Users/hannahvineer/OneDrive - The University of Liverpool/Ciezarek_models/Maxent_output_model_projection/")

# Save the path to the jar file for later
jar = paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

# 5.1 Import environmental data ------------------------------------------------
# Load your rasters
bio.hist = stack(paste0(GIS_dir, "/bio.hist.tif"))
par(mfrow = c(3,3))
plot(bio.hist) # compare this against the "preview" tiff generated above to be sure that the layers are correct
# Rename each layer to avoid long names
names(bio.hist) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.CN.ssp245.2030s = stack(paste0(GIS_dir, "bio.CN.ssp245.2030s.tif"))
names(bio.CN.ssp245.2030s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.CN.ssp370.2030s = stack(paste0(GIS_dir, "bio.CN.ssp370.2030s.tif"))
names(bio.CN.ssp370.2030s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.CN.ssp585.2030s = stack(paste0(GIS_dir, "bio.CN.ssp585.2030s.tif"))
names(bio.CN.ssp585.2030s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.IP.ssp585.2030s = stack(paste0(GIS_dir, "bio.IP.ssp585.2030s.tif"))
names(bio.IP.ssp585.2030s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.CN.ssp245.2050s = stack(paste0(GIS_dir, "bio.CN.ssp245.2050s.tif"))
names(bio.CN.ssp245.2050s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.CN.ssp370.2050s = stack(paste0(GIS_dir, "bio.CN.ssp370.2050s.tif"))
names(bio.CN.ssp370.2050s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.CN.ssp585.2050s = stack(paste0(GIS_dir, "bio.CN.ssp585.2050s.tif"))
names(bio.CN.ssp585.2050s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.IP.ssp585.2050s = stack(paste0(GIS_dir, "bio.IP.ssp585.2050s.tif"))
names(bio.IP.ssp585.2050s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.CN.ssp245.2070s = stack(paste0(GIS_dir, "bio.CN.ssp245.2070s.tif"))
names(bio.CN.ssp245.2070s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.CN.ssp370.2070s = stack(paste0(GIS_dir, "bio.CN.ssp370.2070s.tif"))
names(bio.CN.ssp370.2070s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.CN.ssp585.2070s = stack(paste0(GIS_dir, "bio.CN.ssp585.2070s.tif"))
names(bio.CN.ssp585.2070s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.IP.ssp585.2070s = stack(paste0(GIS_dir, "bio.IP.ssp585.2070s.tif"))
names(bio.IP.ssp585.2070s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.CN.ssp245.2090s = stack(paste0(GIS_dir, "bio.CN.ssp245.2090s.tif"))
names(bio.CN.ssp245.2090s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.CN.ssp370.2090s = stack(paste0(GIS_dir, "bio.CN.ssp370.2090s.tif"))
names(bio.CN.ssp370.2090s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.CN.ssp585.2090s = stack(paste0(GIS_dir, "bio.CN.ssp585.2090s.tif"))
names(bio.CN.ssp585.2090s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above

bio.IP.ssp585.2090s = stack(paste0(GIS_dir, "bio.IP.ssp585.2090s.tif"))
names(bio.IP.ssp585.2090s) = c("bio15", "bio18", "bio19", "bio2", "bio3", "bio4","bio8", "bio9","elev") # corresponding to the "preview" tiff generated above


# # Get Europe boundary
# europe = ne_countries(scale = "medium", type = "countries", continent = "Europe")
# europe = crop(europe, c(-25, 45, 30, 85))

# 5.2. Run Maxent model --------------------------------------------------------

# In section 4 we validated the models
# In this section we will run Maxent using 100% of the occurrence records, and the historic data,
# then project onto future climate scenarios

# 5.2.1 Maxent =================================================================

spp = c("Aa", "As", "Dl", "Ef", "Pp", "Vp")
for (j in 1:length(spp)) {
  # Read in the presence-background dataset with appended environmental data
  DataSpecies = read.csv(paste0(getwd(), "/Occurrence_subsets/thinned_final/SWD/", spp[j],"_SWD.csv"))
  
  nRow <- nrow(DataSpecies)

    # Run Maxent
    # WARNING: YOU MUST MAKE SURE THE INDICES FOR THE CALIB AND EVAL OBJECTS 
    # MATCH THE INDICES FOR THE ENVIRONMENTAL DATA COLUMNS
    if (file.exists(jar) & require(rJava)) {
      maxent_mod = dismo::maxent(x = DataSpecies[, 6:14], p = DataSpecies[,3], path = paste0(mod_out_dir, spp[j]), args = c("responsecurves=TRUE", "pictures=TRUE", "jackknife=TRUE"))
      # Predict probability of presence for all scenarios
      for (k in ls(pattern = "bio.")) {
      px = dismo::predict(get(k), maxent_mod, progress='text')
      writeRaster(px, filename = paste0(mod_out_dir, spp[j],"/prediction_", k, ".tif"), format = "GTiff") # Save the model prediction map
      }
    } else {
      print("cannot run this example because MaxEnt is not available")
    }
    print(paste("Finished MAXENT", spp[j], sep = " ")) # Shows progress (because this loop can be slow)
  }

# Note that the training AUC for each validation replicate is within the folder for
# that run (in the maxent results csv file). The test AUC is in the evaluation 
# csv file that contains the evaluation results for all of the validation
# runs.



# # Qualitative (visual) validation using the lower accuracy records
# tiff("Maxent_output/Maxent_output_v2/Overlay_all_data_maps_Dmarginatus.tif", width = 8, height = 6, units = "in", res = 300)
# par(mfrow = c(1,2))
# plot(dm.mean, ylab = "Latitude", xlab = "Longitude", main = expression(italic("D. marginatus")))
# mtext(text = "Mean", side = 3, line = 0, cex = 0.8)
# drehmen.dat = read.csv("Occurrence_data/Data from Drehmann et al 2020.csv", sep = ";")
# names(drehmen.dat)
# head(drehmen.dat)
# new.dat = read.csv("Occurrence_data/Dermacentor data Germany unpublished.csv", sep = ";")
# names(new.dat)
# head(new.dat)
# # Add column describing the data source
# new.dat$source = "unpublished"
# drehmen.dat$source = "drehmen"
# dat = rbind(drehmen.dat[,c(1,2,3,9,10,11,31)], new.dat[,c(1,2,3,9,10,11,18)])
# points(dat$Latitude[which(dat$Species=="Dermacentor marginatus")]~dat$Longitude[which(dat$Species=="Dermacentor marginatus")], add = T, pch = 20, cex = 0.2)
# plot(dm.agreement, ylab = "Latitude", xlab = "Longitude", main = expression(italic("D. marginatus")))
# mtext(text = "Model agreement", side = 3, line = 0, cex = 0.8)
# points(dat$Latitude[which(dat$Species=="Dermacentor marginatus")]~dat$Longitude[which(dat$Species=="Dermacentor marginatus")], add = T, pch = 20, cex = 0.2)
# dev.off()
# 
# # Qualitative (visual) validation using the lower accuracy records
# tiff("Maxent_output/Maxent_output_v2/Overlay_all_data_maps_Dreticulatus.tif", width = 8, height = 6, units = "in", res = 300)
# par(mfrow = c(1,2))
# plot(dr.mean, ylab = "Latitude", xlab = "Longitude", main = expression(italic("D. reticulatus")))
# mtext(text = "Mean", side = 3, line = 0, cex = 0.8)
# drehmen.dat = read.csv("Occurrence_data/Data from Drehmann et al 2020.csv", sep = ";")
# names(drehmen.dat)
# head(drehmen.dat)
# new.dat = read.csv("Occurrence_data/Dermacentor data Germany unpublished.csv", sep = ";")
# names(new.dat)
# head(new.dat)
# # Add column describing the data source
# new.dat$source = "unpublished"
# drehmen.dat$source = "drehmen"
# dat = rbind(drehmen.dat[,c(1,2,3,9,10,11,31)], new.dat[,c(1,2,3,9,10,11,18)])
# points(dat$Latitude[which(dat$Species=="Dermacentor reticulatus")]~dat$Longitude[which(dat$Species=="Dermacentor reticulatus")], add = T, pch = 20, cex = 0.2)
# plot(dr.agreement, ylab = "Latitude", xlab = "Longitude", main = expression(italic("D. reticulatus")))
# mtext(text = "Model agreement", side = 3, line = 0, cex = 0.8)
# points(dat$Latitude[which(dat$Species=="Dermacentor reticulatus")]~dat$Longitude[which(dat$Species=="Dermacentor reticulatus")], add = T, pch = 20, cex = 0.2)
# dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### 6. Collate SDM Output ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 6.1 Make SDM output into shapefiles --------------------------------------------------------
setwd("/Users/pippin/Documents/Reindeer_2022/")
As.results <- read.csv("As_Validation_Train_Results.csv")
Aa.results <- read.csv("Aa_Validation_Train_Results.csv")
Ef.results <- read.csv("Ef_Validation_Train_Results.csv")
Dl.results <- read.csv("Dl_Validation_Train_Results.csv")
Pp.results <- read.csv("Pp_Validation_Train_Results.csv")
Vp.results <- read.csv("Vp_Validation_Train_Results.csv")

As.thresh <- mean(As.results$Equate.entropy.of.thresholded.and.original.distributions.Cloglog.threshold)
Aa.thresh <- mean(Aa.results$Equate.entropy.of.thresholded.and.original.distributions.Cloglog.threshold)
Ef.thresh <- mean(Ef.results$Equate.entropy.of.thresholded.and.original.distributions.Cloglog.threshold)
Dl.thresh <- mean(Dl.results$Equate.entropy.of.thresholded.and.original.distributions.Cloglog.threshold)
Pp.thresh <- mean(Pp.results$Equate.entropy.of.thresholded.and.original.distributions.Cloglog.threshold)
Vp.thresh <- mean(Vp.results$Equate.entropy.of.thresholded.and.original.distributions.Cloglog.threshold)

As.mat <- matrix(c(0, As.thresh, 0, As.thresh, 1, 1), ncol = 3, byrow = TRUE)
Aa.mat <- matrix(c(0, Aa.thresh, 0, Aa.thresh, 1, 1), ncol = 3, byrow = TRUE)
Ef.mat <- matrix(c(0, Ef.thresh, 0, Ef.thresh, 1, 1), ncol = 3, byrow = TRUE)
Dl.mat <- matrix(c(0, Dl.thresh, 0, Dl.thresh, 1, 1), ncol = 3, byrow = TRUE)
Pp.mat <- matrix(c(0, Pp.thresh, 0, Pp.thresh, 1, 1), ncol = 3, byrow = TRUE)
Vp.mat <- matrix(c(0, Vp.thresh, 0, Vp.thresh, 1, 1), ncol = 3, byrow = TRUE)

setwd("/Users/pippin/Documents/Reindeer_2022/Figures/")
As.present <- raster("As_mean_probability_of_presence.tif")
Aa.present <- raster("Aa_mean_probability_of_presence.tif")
Ef.present <- raster("Ef_mean_probability_of_presence.tif")
Dl.present <- raster("Dl_mean_probability_of_presence.tif")
Pp.present <- raster("Pp_mean_probability_of_presence.tif")
Vp.present <- raster("Vp_mean_probability_of_presence.tif")

As.present <- reclassify(As.present, As.mat)
Aa.present <- reclassify(Aa.present, Aa.mat)
Ef.present <- reclassify(Ef.present, Ef.mat)
Dl.present <- reclassify(Dl.present, Dl.mat)
Pp.present <- reclassify(Pp.present, Pp.mat)
Vp.present <- reclassify(Vp.present, Vp.mat)

Allspp.present <- As.present + As.present + Ef.present + Dl.present + Pp.present + Vp.present
writeRaster(Allspp.present, filename = "Allspp.present.tif", format="GTiff")

Allspp.present.shp <- rasterToPolygons(As.present, fun = function(x){x>=1})

writeOGR(Allspp.present.shp, "Allspp.present", driver = "ESRI Shapefile", layer = "Allspp.present.shp")

# Combine individual species by climate model ----------------------------------
setwd("/Users/pippin/Documents/Reindeer_2022/")

Aa.future = stack(file.path(paste0(getwd(), "/Aa"), list.files(paste0(getwd(), "/Aa/"))))
Aa.future <- reclassify(Aa.future, Aa.mat)
As.future = stack(file.path(paste0(getwd(), "/As"), list.files(paste0(getwd(), "/As/"))))
As.future <- reclassify(As.future, As.mat)
Dl.future = stack(file.path(paste0(getwd(), "/Dl"), list.files(paste0(getwd(), "/Dl/"))))
Dl.future <- reclassify(Dl.future, Dl.mat)
Ef.future = stack(file.path(paste0(getwd(), "/Ef"), list.files(paste0(getwd(), "/Ef/"))))
Ef.future <- reclassify(Ef.future, Ef.mat)
Pp.future = stack(file.path(paste0(getwd(), "/Pp"), list.files(paste0(getwd(), "/Pp/"))))
Pp.future <- reclassify(Pp.future, Pp.mat)
Vp.future = stack(file.path(paste0(getwd(), "/Vp"), list.files(paste0(getwd(), "/Vp/"))))
Vp.future <- reclassify(Vp.future, Vp.mat)


ssp245.2030s <- Aa.future[[1]] + As.future[[1]] + Dl.future[[1]] + Ef.future[[1]] + Pp.future[[1]] + Vp.future[[1]]
ssp245.2050s <- Aa.future[[2]] + As.future[[2]] + Dl.future[[2]] + Ef.future[[2]] + Pp.future[[2]] + Vp.future[[2]]
ssp245.2070s <- Aa.future[[3]] + As.future[[3]] + Dl.future[[3]] + Ef.future[[3]] + Pp.future[[3]] + Vp.future[[3]]
ssp245.2090s <- Aa.future[[4]] + As.future[[4]] + Dl.future[[4]] + Ef.future[[4]] + Pp.future[[4]] + Vp.future[[4]]

ssp370.2030s <- Aa.future[[5]] + As.future[[5]] + Dl.future[[5]] + Ef.future[[5]] + Pp.future[[5]] + Vp.future[[5]]
ssp370.2050s <- Aa.future[[6]] + As.future[[6]] + Dl.future[[6]] + Ef.future[[6]] + Pp.future[[6]] + Vp.future[[6]]
ssp370.2070s <- Aa.future[[7]] + As.future[[7]] + Dl.future[[7]] + Ef.future[[7]] + Pp.future[[7]] + Vp.future[[7]]
ssp370.2090s <- Aa.future[[8]] + As.future[[8]] + Dl.future[[8]] + Ef.future[[8]] + Pp.future[[8]] + Vp.future[[8]]

ssp585.2030s <- Aa.future[[9]] + As.future[[9]] + Dl.future[[9]] + Ef.future[[9]] + Pp.future[[9]] + Vp.future[[9]]
ssp585.2050s <- Aa.future[[10]] + As.future[[10]] + Dl.future[[10]] + Ef.future[[10]] + Pp.future[[10]] + Vp.future[[10]]
ssp585.2070s <- Aa.future[[11]] + As.future[[11]] + Dl.future[[11]] + Ef.future[[11]] + Pp.future[[11]] + Vp.future[[11]]
ssp585.2090s <- Aa.future[[12]] + As.future[[12]] + Dl.future[[12]] + Ef.future[[12]] + Pp.future[[12]] + Vp.future[[12]]

lapply(X = ls(pattern = "ssp"), FUN = function (x) writeRaster(get(x), filename = paste0(x, ".tif"), format="GTiff"))

# 6.2 Compare SDM outputs  ----------------------------------
# 6.2.1 Species richness
# Europe
spec.rich.df <- cbind((as.data.frame(table(as.vector(ssp245.2030s)))),
                      (as.data.frame(table(as.vector(ssp245.2050s))))[,2],
                      (as.data.frame(table(as.vector(ssp245.2070s))))[,2],
                      (as.data.frame(table(as.vector(ssp245.2090s))))[,2],
                      (as.data.frame(table(as.vector(ssp370.2030s))))[,2],
                      (as.data.frame(table(as.vector(ssp370.2050s))))[,2],
                      (as.data.frame(table(as.vector(ssp370.2070s))))[,2],
                      (as.data.frame(table(as.vector(ssp370.2090s))))[,2],
                      (as.data.frame(table(as.vector(ssp585.2030s))))[,2],
                      (as.data.frame(table(as.vector(ssp585.2050s))))[,2],
                      (as.data.frame(table(as.vector(ssp585.2070s))))[,2],
                      (as.data.frame(table(as.vector(ssp585.2090s))))[,2])


present <- raster("Allspp.present.tif")
spec.rich.df <- cbind(spec.rich.df, (as.data.frame(table(as.vector(present))))[,2])

colnames(spec.rich.df) <- c("richness", "ssp245.2030s", "ssp245.2050s", "ssp245.2070s", "ssp245.2090s", 
                            "ssp370.2030s", "ssp370.2050s", "ssp370.2070s", "ssp370.2090s",
                            "ssp585.2030s", "ssp585.2050s", "ssp585.2070s", "ssp585.2090s", "historic")


spec.rich.df <- pivot_longer(spec.rich.df, ssp245.2030s:historic, names_to = "projection")
spec.rich.df <- as.data.frame(spec.rich.df)
spec.rich.df$projection <- as.factor(spec.rich.df$projection)

png("Species_richness_Europe.png")
ggplot(spec.rich.df, aes(projection, value)) + geom_col(aes(fill = richness)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Projection", y = "Number of pixels")+
  scale_fill_brewer("Number of species")
dev.off()

# Crop by reindeer distribution
g.dist <- stack(present, ssp245.2030s, ssp245.2050s, ssp245.2070s, ssp245.2090s,
                ssp370.2030s, ssp370.2050s, ssp370.2070s, ssp370.2090s,
                ssp585.2030s, ssp585.2050s, ssp585.2070s, ssp585.2090s)

r.dist <- mask(g.dist, reindeer)

r.spec.rich.df <- cbind((as.data.frame(table(as.vector(r.dist[[1]])))),
                        (as.data.frame(table(as.vector(r.dist[[2]]))))[,2],
                        (as.data.frame(table(as.vector(r.dist[[3]]))))[,2],
                        (as.data.frame(table(as.vector(r.dist[[4]]))))[,2],
                        (as.data.frame(table(as.vector(r.dist[[5]]))))[,2],
                        (as.data.frame(table(as.vector(r.dist[[6]]))))[,2],
                        (as.data.frame(table(as.vector(r.dist[[7]]))))[,2],
                        (as.data.frame(table(as.vector(r.dist[[8]]))))[,2],
                        c(0, (as.data.frame(table(as.vector(r.dist[[9]]))))[,2]), # rasters 9 and 13 have no pixels with 0 species
                        (as.data.frame(table(as.vector(r.dist[[10]]))))[,2],
                        (as.data.frame(table(as.vector(r.dist[[11]]))))[,2],
                        (as.data.frame(table(as.vector(r.dist[[12]]))))[,2],
                        c(0, (as.data.frame(table(as.vector(r.dist[[13]]))))[,2]))

colnames(r.spec.rich.df) <- c("richness", "historic", "ssp245.2030s", "ssp245.2050s", "ssp245.2070s", "ssp245.2090s", 
                              "ssp370.2030s", "ssp370.2050s", "ssp370.2070s", "ssp370.2090s",
                              "ssp585.2030s", "ssp585.2050s", "ssp585.2070s", "ssp585.2090s")

r.spec.rich.df <- pivot_longer(r.spec.rich.df, historic:ssp585.2090s, names_to = "projection")
spec.rich.df <- as.data.frame(spec.rich.df)
spec.rich.df$projection <- as.factor(spec.rich.df$projection)

png("Species_richness_Reindeer.png")
ggplot(r.spec.rich.df, aes(projection, value)) + geom_col(aes(fill = richness)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Projection", y = "Number of pixels")+
  scale_fill_brewer("Number of species")
dev.off()

# 6.2.2 Individual species ranges  ---------------------------------------------
# Europe
Aa.future.df <- data.frame(matrix(nrow = 2, ncol = 12))
for (i in 1:nlayers(Aa.future)){
  Aa.future.df[,i] <- as.data.frame(table(as.vector(Aa.future[[i]])))[,2]
}
Aa.future.df <- cbind(Aa.future.df, as.data.frame(table(as.vector(Aa.present)))[,2])
Aa.area <- Aa.future.df[2,]/(Aa.future.df[1,]+Aa.future.df[2,])*100 
Aa.area <- cbind("Arianta arbustorum", Aa.area)

As.future.df <- data.frame(matrix(nrow = 2, ncol = 12))
for (i in 1:nlayers(As.future)){
  As.future.df[,i] <- as.data.frame(table(as.vector(As.future[[i]])))[,2]
}
As.future.df <- cbind(As.future.df, as.data.frame(table(as.vector(As.present)))[,2])
As.area <- As.future.df[2,]/(As.future.df[1,]+As.future.df[2,])*100 
As.area <- cbind("Arion subfuscus", As.area)

Dl.future.df <- data.frame(matrix(nrow = 2, ncol = 12))
for (i in 1:nlayers(Dl.future)){
  Dl.future.df[,i] <- as.data.frame(table(as.vector(Dl.future[[i]])))[,2]
}
Dl.future.df <- cbind(Dl.future.df, as.data.frame(table(as.vector(Dl.present)))[,2])
Dl.area <- Dl.future.df[2,]/(Dl.future.df[1,]+Dl.future.df[2,])*100 
Dl.area <- cbind("Deroceras laeve", Dl.area)


Ef.future.df <- data.frame(matrix(nrow = 2, ncol = 12))
for (i in 1:nlayers(Ef.future)){
  Ef.future.df[,i] <- as.data.frame(table(as.vector(Ef.future[[i]])))[,2]
}
Ef.future.df <- cbind(Ef.future.df, as.data.frame(table(as.vector(Ef.present)))[,2])
Ef.area <- Ef.future.df[2,]/(Ef.future.df[1,]+Ef.future.df[2,])*100 
Ef.area <- cbind("Euconulus fulvus", Ef.area)


Pp.future.df <- data.frame(matrix(nrow = 2, ncol = 12))
for (i in 1:nlayers(Pp.future)){
  Pp.future.df[,i] <- as.data.frame(table(as.vector(Pp.future[[i]])))[,2]
}
Pp.future.df <- cbind(Pp.future.df, as.data.frame(table(as.vector(Pp.present)))[,2])
Pp.area <- Pp.future.df[2,]/(Pp.future.df[1,]+Pp.future.df[2,])*100 
Pp.area <- cbind("Punctum pygmaeum", Pp.area)


Vp.future.df <- data.frame(matrix(nrow = 2, ncol = 12))
for (i in 1:nlayers(Vp.future)){
  Vp.future.df[,i] <- as.data.frame(table(as.vector(Vp.future[[i]])))[,2]
}
Vp.future.df <- cbind(Vp.future.df, as.data.frame(table(as.vector(Vp.present)))[,2])
Vp.area <- Vp.future.df[2,]/(Vp.future.df[1,]+Vp.future.df[2,])*100 
Vp.area <- cbind("Vitrina pellucida", Vp.area)


column.names <- c("species","ssp245.2030s", "ssp245.2050s", "ssp245.2070s", "ssp245.2090s", 
                  "ssp370.2030s", "ssp370.2050s", "ssp370.2070s", "ssp370.2090s",
                  "ssp585.2030s", "ssp585.2050s", "ssp585.2070s", "ssp585.2090s", "historic")

colnames(Aa.area) <- column.names
colnames(As.area) <- column.names
colnames(Dl.area) <- column.names
colnames(Ef.area) <- column.names
colnames(Pp.area) <- column.names
colnames(Vp.area) <- column.names

spec.area <- rbind(Aa.area, As.area, Dl.area, Ef.area, Pp.area, Vp.area)

spec.area <- pivot_longer(spec.area, ssp245.2030s:historic, names_to = "projection")
spec.area$species <- as.factor(spec.area$species)

png("Species_area_Europe.png")
ggplot(spec.area, aes(projection, value)) + geom_col() + facet_wrap(spec.area$species) +
  theme(strip.text = element_text(face = "italic")) +
  theme(axis.text.x = element_text(angle = 90)) + labs(x = "Projection", y = "% Area") 
dev.off()  

# Reindeer area
Aa.future.r <- mask(Aa.future, reindeer)
As.future.r <- mask(As.future, reindeer)
Dl.future.r <- mask(Dl.future, reindeer)
Ef.future.r <- mask(Ef.future, reindeer)
Pp.future.r <- mask(Pp.future, reindeer)
Vp.future.r <- mask(Vp.future, reindeer)

Aa.future.r.df <- data.frame(matrix(nrow = 2, ncol = 12))
for (i in 1:nlayers(Aa.future.r)){
  Aa.future.r.df[,i] <- as.data.frame(table(as.vector(Aa.future.r[[i]])))[,2]
}
Aa.future.r.df <- cbind(Aa.future.r.df, as.data.frame(table(as.vector(Aa.present)))[,2])
Aa.area.r <- Aa.future.r.df[2,]/(Aa.future.r.df[1,]+Aa.future.r.df[2,])*100 
Aa.area.r <- cbind("Arianta arbustorum", Aa.area.r)

As.future.r.df <- data.frame(matrix(nrow = 2, ncol = 12))
for (i in 1:nlayers(As.future.r)){
  As.future.r.df[,i] <- as.data.frame(table(as.vector(As.future.r[[i]])))[,2]
}
As.future.r.df <- cbind(As.future.r.df, as.data.frame(table(as.vector(As.present)))[,2])
As.area.r <- As.future.r.df[2,]/(As.future.r.df[1,]+As.future.r.df[2,])*100 
As.area.r <- cbind("Arion subfuscus", As.area.r)

Dl.future.r.df <- data.frame(matrix(nrow = 2, ncol = 12))
for (i in 1:nlayers(Dl.future.r)){
  Dl.future.r.df[,i] <- as.data.frame(table(as.vector(Dl.future.r[[i]])))[,2]
}
Dl.future.r.df <- cbind(Dl.future.r.df, as.data.frame(table(as.vector(Dl.present)))[,2])
Dl.area.r <- Dl.future.r.df[2,]/(Dl.future.r.df[1,]+Dl.future.r.df[2,])*100 
Dl.area.r <- cbind("Deroceras laeve", Dl.area.r)


Ef.future.r.df <- data.frame(matrix(nrow = 2, ncol = 12))
for (i in 1:nlayers(Ef.future.r)){
  Ef.future.r.df[,i] <- as.data.frame(table(as.vector(Ef.future.r[[i]])))[,2]
}
Ef.future.r.df <- cbind(Ef.future.r.df, as.data.frame(table(as.vector(Ef.present)))[,2])
Ef.area.r <- Ef.future.r.df[2,]/(Ef.future.r.df[1,]+Ef.future.r.df[2,])*100 
Ef.area.r <- cbind("Euconulus fulvus", Ef.area.r)


Pp.future.r.df <- data.frame(matrix(nrow = 2, ncol = 12))
for (i in 1:nlayers(Pp.future.r)){
  Pp.future.r.df[,i] <- as.data.frame(table(as.vector(Pp.future.r[[i]])))[,2]
}
Pp.future.r.df <- cbind(Pp.future.r.df, as.data.frame(table(as.vector(Pp.present)))[,2])
Pp.area.r <- Pp.future.r.df[2,]/(Pp.future.r.df[1,]+Pp.future.r.df[2,])*100 
Pp.area.r <- cbind("Punctum pygmaeum", Pp.area.r)


Vp.future.r.df <- data.frame(matrix(nrow = 2, ncol = 12))
for (i in 1:nlayers(Vp.future.r)){
  Vp.future.r.df[,i] <- as.data.frame(table(as.vector(Vp.future.r[[i]])))[,2]
}
Vp.future.r.df <- cbind(Vp.future.r.df, as.data.frame(table(as.vector(Vp.present)))[,2])
Vp.area.r <- Vp.future.r.df[2,]/(Vp.future.r.df[1,]+Vp.future.r.df[2,])*100 
Vp.area.r <- cbind("Vitrina pellucida", Vp.area.r)

colnames(Aa.area.r) <- column.names
colnames(As.area.r) <- column.names
colnames(Dl.area.r) <- column.names
colnames(Ef.area.r) <- column.names
colnames(Pp.area.r) <- column.names
colnames(Vp.area.r) <- column.names

spec.area.r <- rbind(Aa.area.r, As.area.r, Dl.area.r, Ef.area.r, Pp.area.r, Vp.area.r)

spec.area.r <- pivot_longer(spec.area.r, ssp245.2030s:historic, names_to = "projection")
spec.area.r$species <- as.factor(spec.area.r$species)

png("Species_area_Reindeer.png")
ggplot(spec.area.r, aes(projection, value)) + geom_col() + facet_wrap(spec.area$species) +
  theme(strip.text = element_text(face = "italic")) +
  theme(axis.text.x = element_text(angle = 90)) + labs(x = "Projection", y = "% Area") 
dev.off() 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#### 7. Degree-day model ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 7.1 Degree-day calculation from Rose Vineer et al. 2021 doi:10.3389/fvets.2020.603990
dev.dat = data.frame(gastropod = c("Euconulus fulvus",
                                   "Euconulus fulvus",
                                   "Euconulus fulvus",
                                   "Euconulus fulvus",
                                   "Euconulus fulvus",
                                   "Arianta arbustorum",
                                   "Arianta arbustorum",
                                   "Arianta arbustorum",
                                   "Arianta arbustorum",
                                   "Arianta arbustorum",
                                   "Succinea pfeifferi",
                                   "Cochlicopa lubrica",
                                   "Discus ruderatus",
                                   "Arion subfuscus",
                                   "Arion hortensis",
                                   "Deroceras reticulatum",
                                   "Deroceras laeve",
                                   "Euconulus fulvus",
                                   "Clausilia bidentata",
                                   "Arianta arbustorum",
                                   "Arianta arbustorum",
                                   "Arianta arbustorum",
                                   "Arianta arbustorum",
                                   "Arianta arbustorum"),
                     temperature = c(12, 16, 20, 24, 28, 12, 16, 20, 24, 28, rep(20, 14)),
                     days_to_L3 = c(60, 30, 20, 15, 15, 75, 35, 20, 15, 12, 18, 24, 18, 18, 18, 12, 18, 18, 24, 18, 21, 21, 21, 28),
                     source = c(rep("H&S", 10), rep("S&H", 10), rep("S", 4)))

dev.dat$DD = (dev.dat$temperature-8)*dev.dat$days_to_L3

# 7.2 Read in files and pre-process ---------------------------------------
worldmap <- ne_countries(scale = 'medium', type = 'map_units')
europe = ne_countries(scale = "medium", type = "countries", continent = c("Europe", "Asia"))
europe = crop(europe, c(-25, 45, 35, 85))
justnorway = worldmap[176,]
reindeer <- readOGR(dsn = "/Users/pippin/Documents/MPhil/Norwegian_ReindeerAreas/AllReindeerAreas3.shp" , p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")

setwd("/Volumes/Backup Plus/Future_worldclim")

max.CN.ssp245.2030s <- stack("wc2.1_30s_tmax_CNRM-ESM2-1_ssp245_2021-2040.tif")
min.CN.ssp245.2030s <- stack("wc2.1_30s_tmin_CNRM-ESM2-1_ssp245_2021-2040.tif")
max.CN.ssp245.2050s <- stack("wc2.1_30s_tmax_CNRM-ESM2-1_ssp245_2041-2060.tif")
min.CN.ssp245.2050s <- stack("wc2.1_30s_tmin_CNRM-ESM2-1_ssp245_2041-2060.tif")
max.CN.ssp245.2070s <- stack("wc2.1_30s_tmax_CNRM-ESM2-1_ssp245_2061-2080.tif")
min.CN.ssp245.2070s <- stack("wc2.1_30s_tmin_CNRM-ESM2-1_ssp245_2061-2080.tif")
max.CN.ssp245.2090s <- stack("wc2.1_30s_tmax_CNRM-ESM2-1_ssp245_2081-2100.tif")
min.CN.ssp245.2090s <- stack("wc2.1_30s_tmin_CNRM-ESM2-1_ssp245_2081-2100.tif")
max.CN.ssp370.2030s <- stack("wc2.1_30s_tmax_CNRM-ESM2-1_ssp370_2021-2040.tif")
min.CN.ssp370.2030s <- stack("wc2.1_30s_tmin_CNRM-ESM2-1_ssp370_2021-2040.tif")
max.CN.ssp370.2050s <- stack("wc2.1_30s_tmax_CNRM-ESM2-1_ssp370_2041-2060.tif")
min.CN.ssp370.2050s <- stack("wc2.1_30s_tmin_CNRM-ESM2-1_ssp370_2041-2060.tif")
max.CN.ssp370.2070s <- stack("wc2.1_30s_tmax_CNRM-ESM2-1_ssp370_2061-2080.tif")
min.CN.ssp370.2070s <- stack("wc2.1_30s_tmin_CNRM-ESM2-1_ssp370_2061-2080.tif")
max.CN.ssp370.2090s <- stack("wc2.1_30s_tmax_CNRM-ESM2-1_ssp370_2081-2100.tif")
min.CN.ssp370.2090s <- stack("wc2.1_30s_tmin_CNRM-ESM2-1_ssp370_2081-2100.tif")
max.CN.ssp585.2030s <- stack("wc2.1_30s_tmax_CNRM-ESM2-1_ssp585_2021-2040.tif")
min.CN.ssp585.2030s <- stack("wc2.1_30s_tmin_CNRM-ESM2-1_ssp585_2021-2040.tif")
max.CN.ssp585.2050s <- stack("wc2.1_30s_tmax_CNRM-ESM2-1_ssp585_2041-2060.tif")
min.CN.ssp585.2050s <- stack("wc2.1_30s_tmin_CNRM-ESM2-1_ssp585_2041-2060.tif")
max.CN.ssp585.2070s <- stack("wc2.1_30s_tmax_CNRM-ESM2-1_ssp585_2061-2080.tif")
min.CN.ssp585.2070s <- stack("wc2.1_30s_tmin_CNRM-ESM2-1_ssp585_2061-2080.tif")
max.CN.ssp585.2090s <- stack("wc2.1_30s_tmax_CNRM-ESM2-1_ssp585_2081-2100.tif")
min.CN.ssp585.2090s <- stack("wc2.1_30s_tmin_CNRM-ESM2-1_ssp585_2081-2100.tif")

setwd("Users/pippin/documents/Reindeer_2022/")
# Crop and mask before processing to speed things up
max.CN.ssp245.2030s <- crop(max.CN.ssp245.2030s, extent(europe), progress = "text")
min.CN.ssp245.2030s <- crop(min.CN.ssp245.2030s, extent(europe), progress = "text")
max.CN.ssp245.2050s <- crop(max.CN.ssp245.2050s, extent(europe), progress = "text")
min.CN.ssp245.2050s <- crop(min.CN.ssp245.2050s, extent(europe), progress = "text")
max.CN.ssp245.2070s <- crop(max.CN.ssp245.2070s, extent(europe), progress = "text")
min.CN.ssp245.2070s <- crop(min.CN.ssp245.2070s, extent(europe), progress = "text")
max.CN.ssp245.2090s <- crop(max.CN.ssp245.2090s, extent(europe), progress = "text")
min.CN.ssp245.2090s <- crop(min.CN.ssp245.2090s, extent(europe), progress = "text")
max.CN.ssp370.2030s <- crop(max.CN.ssp370.2030s, extent(europe), progress = "text")
min.CN.ssp370.2030s <- crop(min.CN.ssp370.2030s, extent(europe), progress = "text")
max.CN.ssp370.2050s <- crop(max.CN.ssp370.2050s, extent(europe), progress = "text")
min.CN.ssp370.2050s <- crop(min.CN.ssp370.2050s, extent(europe), progress = "text")
max.CN.ssp370.2070s <- crop(max.CN.ssp370.2070s, extent(europe), progress = "text")
min.CN.ssp370.2070s <- crop(min.CN.ssp370.2070s, extent(europe), progress = "text")
max.CN.ssp370.2090s <- crop(max.CN.ssp370.2090s, extent(europe), progress = "text")
min.CN.ssp370.2090s <- crop(min.CN.ssp370.2090s, extent(europe), progress = "text")
max.CN.ssp585.2030s <- crop(max.CN.ssp585.2030s, extent(europe), progress = "text")
min.CN.ssp585.2030s <- crop(min.CN.ssp585.2030s, extent(europe), progress = "text")
max.CN.ssp585.2050s <- crop(max.CN.ssp585.2050s, extent(europe), progress = "text")
min.CN.ssp585.2050s <- crop(min.CN.ssp585.2050s, extent(europe), progress = "text")
max.CN.ssp585.2070s <- crop(max.CN.ssp585.2070s, extent(europe), progress = "text")
min.CN.ssp585.2070s <- crop(min.CN.ssp585.2070s, extent(europe), progress = "text")
max.CN.ssp585.2090s <- crop(max.CN.ssp585.2090s, extent(europe), progress = "text")
min.CN.ssp585.2090s <- crop(min.CN.ssp585.2090s, extent(europe), progress = "text")

max.CN.ssp245.2030s <- mask(max.CN.ssp245.2030s, europe, progress = "text")
min.CN.ssp245.2030s <- mask(min.CN.ssp245.2030s, europe, progress = "text")
max.CN.ssp245.2050s <- mask(max.CN.ssp245.2050s, europe, progress = "text")
min.CN.ssp245.2050s <- mask(min.CN.ssp245.2050s, europe, progress = "text")
max.CN.ssp245.2070s <- mask(max.CN.ssp245.2070s, europe, progress = "text")
min.CN.ssp245.2070s <- mask(min.CN.ssp245.2070s, europe, progress = "text")
max.CN.ssp245.2090s <- mask(max.CN.ssp245.2090s, europe, progress = "text")
min.CN.ssp245.2090s <- mask(min.CN.ssp245.2090s, europe, progress = "text")
max.CN.ssp370.2030s <- mask(max.CN.ssp370.2030s, europe, progress = "text")
min.CN.ssp370.2030s <- mask(min.CN.ssp370.2030s, europe, progress = "text")
max.CN.ssp370.2050s <- mask(max.CN.ssp370.2050s, europe, progress = "text")
min.CN.ssp370.2050s <- mask(min.CN.ssp370.2050s, europe, progress = "text")
max.CN.ssp370.2070s <- mask(max.CN.ssp370.2070s, europe, progress = "text")
min.CN.ssp370.2070s <- mask(min.CN.ssp370.2070s, europe, progress = "text")
max.CN.ssp370.2090s <- mask(max.CN.ssp370.2090s, europe, progress = "text")
min.CN.ssp370.2090s <- mask(min.CN.ssp370.2090s, europe, progress = "text")
max.CN.ssp585.2030s <- mask(max.CN.ssp585.2030s, europe, progress = "text")
min.CN.ssp585.2030s <- mask(min.CN.ssp585.2030s, europe, progress = "text")
max.CN.ssp585.2050s <- mask(max.CN.ssp585.2050s, europe, progress = "text")
min.CN.ssp585.2050s <- mask(min.CN.ssp585.2050s, europe, progress = "text")
max.CN.ssp585.2070s <- mask(max.CN.ssp585.2070s, europe, progress = "text")
min.CN.ssp585.2070s <- mask(min.CN.ssp585.2070s, europe, progress = "text")
max.CN.ssp585.2090s <- mask(max.CN.ssp585.2090s, europe, progress = "text")
min.CN.ssp585.2090s <- mask(min.CN.ssp585.2090s, europe, progress = "text")

# Calculate mean monthly temperatures
mean.CN.ssp245.2030s <- (max.CN.ssp245.2030s + min.CN.ssp245.2030s)/2
mean.CN.ssp245.2050s <- (max.CN.ssp245.2050s + min.CN.ssp245.2050s)/2
mean.CN.ssp245.2070s <- (max.CN.ssp245.2070s + min.CN.ssp245.2070s)/2
mean.CN.ssp245.2090s <- (max.CN.ssp245.2090s + min.CN.ssp245.2090s)/2
mean.CN.ssp370.2030s <- (max.CN.ssp370.2030s + min.CN.ssp370.2030s)/2
mean.CN.ssp370.2050s <- (max.CN.ssp370.2050s + min.CN.ssp370.2050s)/2
mean.CN.ssp370.2070s <- (max.CN.ssp370.2070s + min.CN.ssp370.2070s)/2
mean.CN.ssp370.2090s <- (max.CN.ssp370.2090s + min.CN.ssp370.2090s)/2
mean.CN.ssp585.2030s <- (max.CN.ssp585.2030s + min.CN.ssp585.2030s)/2
mean.CN.ssp585.2050s <- (max.CN.ssp585.2050s + min.CN.ssp585.2050s)/2
mean.CN.ssp585.2070s <- (max.CN.ssp585.2070s + min.CN.ssp585.2070s)/2
mean.CN.ssp585.2090s <- (max.CN.ssp585.2090s + min.CN.ssp585.2090s)/2

# Save rasters
setwd("ENTER THE DIRECTORY WHERE YOU WANT TO SAVE RASTERS HERE")
setwd("/Users/pippin/Reindeer_2022/DD_Europe/")
lapply(X = ls(pattern = "mean."), FUN = function (x) writeRaster(get(x), filename = paste0(x, ".tif"), format="GTiff"))


# Read in mean rasters that were saved
setwd("/Volumes/Backup Plus/")

mean.CN.ssp245.2030s <- stack("mean.CN.ssp245.2030s.tif")
mean.CN.ssp245.2050s <- stack("mean.CN.ssp245.2050s.tif")
mean.CN.ssp245.2070s <- stack("mean.CN.ssp245.2070s.tif")
mean.CN.ssp245.2090s <- stack("mean.CN.ssp245.2090s.tif")
mean.CN.ssp370.2030s <- stack("mean.CN.ssp370.2030s.tif")
mean.CN.ssp370.2050s <- stack("mean.CN.ssp370.2050s.tif")
mean.CN.ssp370.2070s <- stack("mean.CN.ssp370.2070s.tif")
mean.CN.ssp370.2090s <- stack("mean.CN.ssp370.2090s.tif")
mean.CN.ssp585.2030s <- stack("mean.CN.ssp585.2030s.tif")
mean.CN.ssp585.2050s <- stack("mean.CN.ssp585.2050s.tif")
mean.CN.ssp585.2070s <- stack("mean.CN.ssp585.2070s.tif")
mean.CN.ssp585.2090s <- stack("mean.CN.ssp585.2090s.tif")


# Change temperatures above 21 to 21
mean.CN.ssp245.2030s[mean.CN.ssp245.2030s > 21] <- 21
mean.CN.ssp245.2050s[mean.CN.ssp245.2050s > 21] <- 21
mean.CN.ssp245.2070s[mean.CN.ssp245.2070s > 21] <- 21
mean.CN.ssp245.2090s[mean.CN.ssp245.2090s > 21] <- 21
mean.CN.ssp370.2030s[mean.CN.ssp370.2030s > 21] <- 21
mean.CN.ssp370.2050s[mean.CN.ssp370.2050s > 21] <- 21
mean.CN.ssp370.2070s[mean.CN.ssp370.2070s > 21] <- 21
mean.CN.ssp370.2090s[mean.CN.ssp370.2090s > 21] <- 21
mean.CN.ssp585.2030s[mean.CN.ssp585.2030s > 21] <- 21
mean.CN.ssp585.2050s[mean.CN.ssp585.2050s > 21] <- 21
mean.CN.ssp585.2070s[mean.CN.ssp585.2070s > 21] <- 21
mean.CN.ssp585.2090s[mean.CN.ssp585.2090s > 21] <- 21

# 7.3 Degree-day calculation for all projections -------------------------
# calculate daily degree days by subtracting min temp for development 
# then removing negative values
# then multiplying monthly figures by number of days in the month
days.per.month <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

dd.CN.ssp245.2030s <- calc(mean.CN.ssp245.2030s, fun=function(x){x-8})
dd.CN.ssp245.2030s[dd.CN.ssp245.2030s<0] <- 0
dd.CN.ssp245.2030s <- dd.CN.ssp245.2030s*days.per.month
dd.CN.ssp245.2050s <- calc(mean.CN.ssp245.2050s, fun=function(x){x-8})
dd.CN.ssp245.2050s[dd.CN.ssp245.2050s<0] <- 0
dd.CN.ssp245.2050s <- dd.CN.ssp245.2050s*days.per.month
dd.CN.ssp245.2070s <- calc(mean.CN.ssp245.2070s, fun=function(x){x-8})
dd.CN.ssp245.2070s[dd.CN.ssp245.2070s<0] <- 0
dd.CN.ssp245.2070s <- dd.CN.ssp245.2070s*days.per.month
dd.CN.ssp245.2090s <- calc(mean.CN.ssp245.2090s, fun=function(x){x-8})
dd.CN.ssp245.2090s[dd.CN.ssp245.2090s<0] <- 0
dd.CN.ssp245.2090s <- dd.CN.ssp245.2090s*days.per.month

dd.CN.ssp370.2030s <- calc(mean.CN.ssp370.2030s, fun=function(x){x-8})
dd.CN.ssp370.2030s[dd.CN.ssp370.2030s<0] <- 0
dd.CN.ssp370.2030s <- dd.CN.ssp370.2030s*days.per.month
dd.CN.ssp370.2050s <- calc(mean.CN.ssp370.2050s, fun=function(x){x-8})
dd.CN.ssp370.2050s[dd.CN.ssp370.2050s<0] <- 0
dd.CN.ssp370.2050s <- dd.CN.ssp370.2050s*days.per.month
dd.CN.ssp370.2070s <- calc(mean.CN.ssp370.2070s, fun=function(x){x-8})
dd.CN.ssp370.2070s[dd.CN.ssp370.2070s<0] <- 0
dd.CN.ssp370.2070s <- dd.CN.ssp370.2070s*days.per.month
dd.CN.ssp370.2090s <- calc(mean.CN.ssp370.2090s, fun=function(x){x-8})
dd.CN.ssp370.2090s[dd.CN.ssp370.2090s<0] <- 0
dd.CN.ssp370.2090s <- dd.CN.ssp370.2090s*days.per.month

dd.CN.ssp585.2030s <- calc(mean.CN.ssp585.2030s, fun=function(x){x-8})
dd.CN.ssp585.2030s[dd.CN.ssp585.2030s<0] <- 0
dd.CN.ssp585.2030s <- dd.CN.ssp585.2030s*days.per.month
dd.CN.ssp585.2050s <- calc(mean.CN.ssp585.2050s, fun=function(x){x-8})
dd.CN.ssp585.2050s[dd.CN.ssp585.2050s<0] <- 0
dd.CN.ssp585.2050s <- dd.CN.ssp585.2050s*days.per.month
dd.CN.ssp585.2070s <- calc(mean.CN.ssp585.2070s, fun=function(x){x-8})
dd.CN.ssp585.2070s[dd.CN.ssp585.2070s<0] <- 0
dd.CN.ssp585.2070s <- dd.CN.ssp585.2070s*days.per.month
dd.CN.ssp585.2090s <- calc(mean.CN.ssp585.2090s, fun=function(x){x-8})
dd.CN.ssp585.2090s[dd.CN.ssp585.2090s<0] <- 0
dd.CN.ssp585.2090s <- dd.CN.ssp585.2090s*days.per.month

# Save rasters
setwd("ENTER THE DIRECTORY WHERE YOU WANT TO SAVE RASTERS HERE")
setwd("/Users/pippin/Reindeer_2022/DD_Europe/")
lapply(X = ls(pattern = "dd."), FUN = function (x) writeRaster(get(x), filename = paste0(x, ".tif"), format="GTiff"))

# 7.4 Degree-day index projections -------------------------------------------
# Read in dd rasters
setwd("/Users/pippin/Documents/Reindeer_2022/DD_Europe/")
dd.CN.ssp245.2030s <- stack("dd.CN.ssp245.2030s.tif")
dd.CN.ssp245.2050s <- stack("dd.CN.ssp245.2050s.tif")
dd.CN.ssp245.2070s <- stack("dd.CN.ssp245.2070s.tif")
dd.CN.ssp245.2090s <- stack("dd.CN.ssp245.2090s.tif")
dd.CN.ssp370.2030s <- stack("dd.CN.ssp370.2030s.tif")
dd.CN.ssp370.2050s <- stack("dd.CN.ssp370.2050s.tif")
dd.CN.ssp370.2070s <- stack("dd.CN.ssp370.2070s.tif")
dd.CN.ssp370.2090s <- stack("dd.CN.ssp370.2090s.tif")
dd.CN.ssp585.2030s <- stack("dd.CN.ssp585.2030s.tif")
dd.CN.ssp585.2050s <- stack("dd.CN.ssp585.2050s.tif")
dd.CN.ssp585.2070s <- stack("dd.CN.ssp585.2070s.tif")
dd.CN.ssp585.2090s <- stack("dd.CN.ssp585.2090s.tif")

# read in sdm rasters
ssp245.2030s <- raster("ssp245.2030s.tif")
ssp245.2050s <- raster("ssp245.2050s.tif")
ssp245.2070s <- raster("ssp245.2070s.tif")
ssp245.2090s <- raster("ssp245.2090s.tif")
ssp370.2030s <- raster("ssp370.2030s.tif")
ssp370.2050s <- raster("ssp370.2050s.tif")
ssp370.2070s <- raster("ssp370.2070s.tif")
ssp370.2090s <- raster("ssp370.2090s.tif")
ssp585.2030s <- raster("ssp585.2030s.tif")
ssp585.2050s <- raster("ssp585.2050s.tif")
ssp585.2070s <- raster("ssp585.2070s.tif")
ssp585.2090s <- raster("ssp585.2090s.tif")


# Convert gastropod distributions into shape files
#### This wouldnt run on my laptop and I dont know how to increase 
# memory capacity on windows to run it on my desktop so actually did this in QGIS
ssp245.2030s.shp <- rasterToPolygons(ssp245.2030s, fun = function(x){x>=1})
ssp245.2050s.shp <- rasterToPolygons(ssp245.2050s, fun = function(x){x>=1})
ssp245.2070s.shp <- rasterToPolygons(ssp245.2070s, fun = function(x){x>=1})
ssp245.2090s.shp <- rasterToPolygons(ssp245.2090s, fun = function(x){x>=1})
ssp370.2030s.shp <- rasterToPolygons(ssp370.2030s, fun = function(x){x>=1})
ssp370.2050s.shp <- rasterToPolygons(ssp370.2050s, fun = function(x){x>=1})
ssp370.2070s.shp <- rasterToPolygons(ssp370.2070s, fun = function(x){x>=1})
ssp370.2090s.shp <- rasterToPolygons(ssp370.2090s, fun = function(x){x>=1})
ssp585.2030s.shp <- rasterToPolygons(ssp585.2030s, fun = function(x){x>=1})
ssp585.2050s.shp <- rasterToPolygons(ssp585.2050s, fun = function(x){x>=1})
ssp585.2070s.shp <- rasterToPolygons(ssp585.2070s, fun = function(x){x>=1})
ssp585.2090s.shp <- rasterToPolygons(ssp585.2090s, fun = function(x){x>=1})

# Read in gastropod shapefiles
setwd("/Volumes/Backup Plus/")
setwd("D:/")
g.ssp245.2030s.shp <- readOGR("g.ssp245.2030s.diss.gpkg", p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
g.ssp245.2050s.shp <- readOGR("g.ssp245.2050s.diss.gpkg", p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
g.ssp245.2070s.shp <- readOGR("g.ssp245.2070s.diss.gpkg", p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
g.ssp245.2090s.shp <- readOGR("g.ssp245.2090s.diss.gpkg", p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
g.ssp370.2030s.shp <- readOGR("g.ssp370.2030s.diss.gpkg", p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
g.ssp370.2050s.shp <- readOGR("g.ssp370.2050s.diss.gpkg", p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
g.ssp370.2070s.shp <- readOGR("g.ssp370.2070s.diss.gpkg", p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
g.ssp370.2090s.shp <- readOGR("g.ssp370.2090s.diss.gpkg", p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
g.ssp585.2030s.shp <- readOGR("g.ssp585.2030s.diss.gpkg", p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
g.ssp585.2050s.shp <- readOGR("g.ssp585.2050s.diss.gpkg", p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
g.ssp585.2070s.shp <- readOGR("g.ssp585.2070s.diss.gpkg", p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")
g.ssp585.2090s.shp <- readOGR("g.ssp585.2090s.diss.gpkg", p4s = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 ")


# Calculate annual DD sums
dd.ssp245.2030s.tot <- sum(dd.CN.ssp245.2030s)
dd.ssp245.2050s.tot <- sum(dd.CN.ssp245.2050s)
dd.ssp245.2070s.tot <- sum(dd.CN.ssp245.2070s)
dd.ssp245.2090s.tot <- sum(dd.CN.ssp245.2090s)
dd.ssp370.2030s.tot <- sum(dd.CN.ssp370.2030s)
dd.ssp370.2050s.tot <- sum(dd.CN.ssp370.2050s)
dd.ssp370.2070s.tot <- sum(dd.CN.ssp370.2070s)
dd.ssp370.2090s.tot <- sum(dd.CN.ssp370.2090s)
dd.ssp585.2030s.tot <- sum(dd.CN.ssp585.2030s)
dd.ssp585.2050s.tot <- sum(dd.CN.ssp585.2050s)
dd.ssp585.2070s.tot <- sum(dd.CN.ssp585.2070s)
dd.ssp585.2090s.tot <- sum(dd.CN.ssp585.2090s)

# Mask by gastropod distribution

dd.ssp245.2030s.g <- mask(dd.ssp245.2030s.tot, g.ssp245.2030s.shp)
dd.ssp245.2050s.g <- mask(dd.ssp245.2050s.tot, g.ssp245.2050s.shp)
dd.ssp245.2070s.g <- mask(dd.ssp245.2070s.tot, g.ssp245.2070s.shp)
dd.ssp245.2090s.g <- mask(dd.ssp245.2090s.tot, g.ssp245.2090s.shp)
dd.ssp370.2030s.g <- mask(dd.ssp370.2030s.tot, g.ssp370.2030s.shp)
dd.ssp370.2050s.g <- mask(dd.ssp370.2050s.tot, g.ssp370.2050s.shp)
dd.ssp370.2070s.g <- mask(dd.ssp370.2070s.tot, g.ssp370.2070s.shp)
dd.ssp370.2090s.g <- mask(dd.ssp370.2090s.tot, g.ssp370.2090s.shp)
dd.ssp585.2030s.g <- mask(dd.ssp585.2030s.tot, g.ssp585.2030s.shp)
dd.ssp585.2050s.g <- mask(dd.ssp585.2050s.tot, g.ssp585.2050s.shp)
dd.ssp585.2070s.g <- mask(dd.ssp585.2070s.tot, g.ssp585.2070s.shp)
dd.ssp585.2090s.g <- mask(dd.ssp585.2090s.tot, g.ssp585.2090s.shp)

# crop by reindeer distribution
dd.ssp245.g <- stack(dd.ssp245.2030s.g, dd.ssp245.2050s.g, dd.ssp245.2070s.g, dd.ssp245.2090s.g)
dd.ssp370.g <- stack(dd.ssp370.2030s.g, dd.ssp370.2050s.g, dd.ssp370.2070s.g, dd.ssp370.2090s.g)
dd.ssp585.g <- stack(dd.ssp585.2030s.g, dd.ssp585.2050s.g, dd.ssp585.2070s.g, dd.ssp585.2090s.g)

dd.ssp245.r <- mask(dd.ssp245.g, reindeer)
dd.ssp370.r <- mask(dd.ssp370.g, reindeer)
dd.ssp585.r <- mask(dd.ssp585.g, reindeer)

# Calculate annual dd index
dd.ssp245.r.index = dd.ssp245.r/round(mean(dev.dat$DD))
dd.ssp370.r.index = dd.ssp370.r/round(mean(dev.dat$DD))
dd.ssp585.r.index = dd.ssp585.r/round(mean(dev.dat$DD))

# 7.5 Plot results --------------------------------------------------
worldmap <- ne_countries(scale = 'medium', type = 'map_units')
Europe <- worldmap[worldmap$continent == 'Europe',]
Europe_cropped = Europe[-c(19, 44, 56),] # remove Norway, Sweden and Finland
Europe_no_norway = Europe[-44,]
Norway = worldmap[c(77, 176, 225),] # Norway, Finland, Sweden
justnorway = worldmap[176,]


boundarylayer = list("sp.polygons", as(Europe_cropped, "SpatialPolygons"), col = "black", fill = "light grey")
norwaylayer = list("sp.polygons", Norway, col = "black", first = F, lwd = 2)
boundarylayer2 = list("sp.polygons", as(Europe_no_norway, "SpatialPolygons"), col = "black", fill = "light grey")
norwaylayer2 = list("sp.polygons", justnorway, col = "black", first = F, lwd = 2)


ssp245.plot <- crop(dd.ssp245.r.index, Norway)
png("future_decadal_average_ssp245.png", width = 12, height = 8,units = "in", res = 300)
spplot(ssp245.plot, names.attr = c("2021-2040", "2041-2060", "2061-2080", "2081-2100"),
       col.regions = rev(magma(10)[c(3:10)]), sp.layout = list(norwaylayer2, boundarylayer2),
       scales = list(draw = T), main = "Thermal Suitability Index", xlab = "Longitude",
       ylab = "Latitude", at = c(0,1,2,3,4,5,6,7,8))
dev.off()

ssp370.plot <- crop(dd.ssp370.r.index, Norway)
png("future_decadal_average_ssp370.png", width = 12, height = 8,units = "in", res = 300)
spplot(ssp370.plot, names.attr = c("2021-2040", "2041-2060", "2061-2080", "2081-2100"),
       col.regions = rev(magma(10)[c(3:10)]), sp.layout = list(norwaylayer2, boundarylayer2),
       scales = list(draw = T), main = "Thermal Suitability Index", xlab = "Longitude",
       ylab = "Latitude", at = c(0,1,2,3,4,5,6,7,8))
dev.off()

ssp585.plot <- crop(dd.ssp585.r.index, Norway)
png("future_decadal_average_ssp585.png", width = 12, height = 8,units = "in", res = 300)
spplot(ssp585.plot, names.attr = c("2021-2040", "2041-2060", "2061-2080", "2081-2100"),
       col.regions = rev(magma(10)[c(3:10)]), sp.layout = list(norwaylayer2, boundarylayer2),
       scales = list(draw = T), main = "Thermal Suitability Index", xlab = "Longitude",
       ylab = "Latitude", at = c(0,1,2,3,4,5,6,7,8))
dev.off()


