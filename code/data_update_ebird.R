# data update full ebird dataset
library(tidyverse)
library(raster)
library(ncdf4)
library(sf)
library(prioritizr)

# species list
setwd("/vast/palmer/home.mccleary/jc3893/30x30/")
list = read_csv("species_list_2024.csv")

# get files
setwd("/gpfs/gibbs/pi/jetz/data/species_datasets/occurrence/ebird/ebird_Apr2024/global/raw/")
lf = list.files("all_files/")

# series of loops to break up data
# process 1000 files at a time
for (h in 0:4){ 
#for (i in 124:length(lf)){
for (i in ((h*1000)+1):((h*1000)+1000)){
  try({
  print(i)
    # fuse file types to get all checklist information available
 eb2 = read_delim(paste0("all_files_newcol/", lf[i]))
 eb1 = read_delim(paste0("all_files/", lf[i])) %>%
   left_join(eb2, by="molid") %>%
   # month and year columns
   mutate(month = month(eventDate),
          year = year(eventDate)) %>%
   # filtering checklists by effort, location etc
   filter(all_species_reported == 1,
          effort_distance_km <= 2,
          duration_minutes <= 180,
          number_observers <= 5,
          samplingprotocol %in% c("Traveling", "Stationary"),
          between(latitude, 0, 87.11),
          between(longitude, -179.9, -42.68),
          month %in% c(1,2,6,7,8,12),
          year >= 2004,
          scientificname %in% list$Scientific) %>%
   # select one list per sampling event
   group_by(sampling_event_identifier) %>%
   sample_n(size = 1) %>%
   ungroup() %>%
   # column selection
   dplyr::select(sampling_event_identifier,  scientificname,
                 longitude, latitude, eventDate, number_observers,
                 effort_distance_km, duration_minutes,
                 time_observations_started)
 # combine
 if (i==((h*1000)+1)){eb = eb1}else{
   eb = rbind(eb, eb1)}
  })
}
  write_csv(eb, paste0("/vast/palmer/home.mccleary/jc3893/30x30/ebird_2024_",
            h,"_unfilled.csv"))
}

# rebuild all pieces to one file
for (h in 0:4){
eb = read_csv(paste0("/vast/palmer/home.mccleary/jc3893/30x30/ebird_2024_",
            h,"_unfilled.csv"))
if (h==0){ebird = eb}else{ebird = rbind(ebird, eb)}
}
# divide seasonally
ebird$month = month(ebird$eventDate)
summer = ebird %>% filter(month %in% 6:8)
write_csv(summer, "/vast/palmer/home.mccleary/jc3893/30x30/ebird_2024_summer_unfilled.csv")
winter = ebird %>% filter(month %in% c(1,2,12))
write_csv(winter, "/vast/palmer/home.mccleary/jc3893/30x30/ebird_2024_winter_unfilled.csv")

# zero-filling species absences for each checklist
for (season in c("summer","winter")){
  ebird = read_csv(paste0("/vast/palmer/home.mccleary/jc3893/30x30/ebird_2024_",
                          season,"_unfilled.csv"))
  # select one list per sampling event
  samps <- ebird %>%
    group_by(sampling_event_identifier) %>%
    sample_n(1) %>%
    ungroup()
  samps$scientificname = NULL
  # zero fill by species (loop species)
  for (j in 1:length(list$Scientific)){
    sp = list$Scientific[j]
    ebsub = ebird[ebird$scientificname==sp,]
    pres_rows <- which(samps$sampling_event_identifier %in% ebsub$sampling_event_identifier)
    samps$observation_count <- 0
    samps$observation_count[pres_rows] <- 1
    colnames(samps)[ncol(samps)] <- gsub(" ","_",sp)
  }
  write_csv(samps, paste0("/vast/palmer/home.mccleary/jc3893/30x30/ebird_2024_",
                          season,"_zerofilled.csv"))
  
}


# annotating
# Climate/Topography - checklists
# standard bioclim vars, load and stack
setwd("/gpfs/gibbs/pi/jetz/data/environmental_datasets/")
bio1 = raster("WC/wc2.1_30s_bio_1.tif")
bio4 = raster("WC/wc2.1_30s_bio_4.tif")
bio12 = raster("WC/wc2.1_30s_bio_12.tif")
bio15 = raster("WC/wc2.1_30s_bio_15.tif")
biostack <- stack(bio1, bio4, bio12, bio15) %>%
  round(1)
names(biostack) <- c("bio1","bio4","bio12","bio15")
# load other variables from earthENV, stack
elevr <- raster("gmted/elevation_1KMmd_GMTEDmd.tif") # elevation
evisum <- raster("gmted/EVI_Summer_6to8_resampled_NA30x30_cropped.tif") # summer EVI
eviwin <- raster("gmted/EVI_Winter_11to2_resampled_NA30x30_cropped.tif") # winter EVI
tri <- raster("gmted/tri_1KMmd_GMTEDmd_resampled_masked_NA30x30_cropped.tif") #topographic roughness index
gmtedstack <- stack(evisum, eviwin, tri)
names(gmtedstack) <- c("evisum", "eviwin", "tri")
# Percent Landcover (ESA CCI)
lc <- list.files("esacci/2019/", pattern="-00000000", full.names=T) %>%
  stack() %>%
  round(1)
names(lc) <- c("pland_10_cropland_rainfed", "pland_100_mosaic_tree_shrub",
               "pland_11_cropland_rainfed",  "pland_110_mosaic_herbacious",
               "pland_12_cropland_rainfed", "pland_120_shrubland",
               "pland_121_shrubland", "pland_122_shrubland",
               "pland_130_grassland", "pland_140_lichens_mosses",
               "pland_150_sparse","pland_152_sparse",
               "pland_153_sparse","pland_160_flooded_freshwater",
               "pland_170_flooded_saltwater","pland_180_flooded_shrub",
               "pland_190_urban", "pland_20_cropland_irrigated",
               "pland_200_barren", "pland_201_barren",
               "pland_202_barren", "pland_210_water",
               "pland_220_ice", "pland_30_mosaic_cropland",
               "pland_40_mosaic_natural_veg", "pland_50_evergreen_broadleaf",
               "pland_60_deciduous_broadleaf",  "pland_61_deciduous_broadleaf",
               "pland_62_deciduous_broadleaf",  "pland_70_evergreen_needleleaf",
               "pland_71_evergreen_needleleaf", "pland_72_evergreen_needleleaf",
               "pland_80_deciduous_needleleaf", "pland_81_deciduous_needleleaf",
               "pland_82_deciduous_needleleaf", "pland_90_mixed_forest")
# frank extreme variables, stack, project to our coords system
# ehe = extreme heat; ece = extreme cold
ehe_s = raster("frank_extremes/prop-ehe-avg-bs.tif")
ece_s = raster("frank_extremes/prop-ece-avg-bs.tif")
s_stack = stack(ehe_s, ece_s) %>%
  projectRaster(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
ehe_w = raster("frank_extremes/prop-ehe-avg-bw.tif")
ece_w = raster("frank_extremes/prop-ece-avg-bw.tif")
w_stack = stack(ehe_w, ece_w) %>%
  projectRaster(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")

# spei (already projected)
spei_s = raster("/vast/palmer/home.mccleary/jc3893/30x30/extremes/rasters/spei_summer.tif") 
spei_w = raster("/vast/palmer/home.mccleary/jc3893/30x30/extremes/rasters/spei_winter.tif") 

# ebird annotation
for (season in c("summer","winter")){ # seasonal loop
  # load data, make spatial
  data <- read_csv(paste0("/vast/palmer/home.mccleary/jc3893/30x30/ebird_2024_",
                          season,"_extracted.csv")) %>%
    st_as_sf(coords=c("longitude","latitude"))
  # extract all envi variables, sometimes seasonally
  bioex <- extract(biostack, data)
  elev <- extract(elevr, data)
  gmtedex <- extract(gmtedstack, data)
  lcex <- extract(lc, data)
  if(season=="summer"){
   pdsiex <- extract(pdsi_s, data)
   frankex <- extract(s_stack, data)
   }else{
   pdsiex <- extract(pdsi_w, data)  
   frankex <- extract(w_stack, data)}
  # transform to CEA projection
  data = data %>%
    st_set_crs("+proj=longlat +datum=WGS84 +no_defs ") %>%
    st_transform(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
  data <- data.frame(sf:::as_Spatial(data))
  # data organizing and naming
  data$optional <- NULL
  names(data)[697:698] <- c("x", "y")
  data = cbind(data, biox, elev, gmtedex, lecex, pdsiex, frankex)
  if(season=="summer"){names(data)[750] <- "pdsi_s"}else{names(data)[750] <- "pdsi_w"}
  if(season=="summer"){names(data)[751:752] <- c("ehe_s","ece_s")}else{
    names(data)[751:752] <- c("ehe_w","ece_w")}
  # hour of observation column
  data$time_observations_started = hour(data$time_observations_started)
  # day of year column
  data$day_of_year = yday(data$eventDate)
  # year column
  data$year = year(data$eventDate)
  data$eventDate = NULL
  # save
  write_csv(data, paste0("/vast/palmer/home.mccleary/jc3893/30x30/ebird_2024_",
                         season,"_extracted.csv"))
  gc()
}

# add in SPEI later (it was a later addition and I didn't rewrite code elegantly)
for (season in c("summer","winter")){
  data <- read_csv(paste0("/gpfs/gibbs/pi/jetz/data/temporary_datasets/ebird_extremes_jeremy/ebird_2024_",
                          season,"_extracted.csv")) %>%
    st_as_sf(coords=c("x","y"))
  # extract seasonal SPEI
  if(season=="summer"){
    speiex <- extract(spei_s, data)
  }else{
    speiex <- extract(spei_w, data)}
  # data organization
  data <- data.frame(sf:::as_Spatial(data))
  data$optional <- NULL
  names(data)[751:752] <- c("x", "y")
  data = cbind(data, speiex)
  if(season=="summer"){names(data)[753] <- "spei_s"}else{names(data)[753] <- "spei_w"}
  # save
  write_csv(data, paste0("/gpfs/gibbs/pi/jetz/data/temporary_datasets/ebird_extremes_jeremy/ebird_2024_",
                         season,"_extracted.csv"))
  gc()
}


### Prediction surface annotation using existing surface (add evi- winter)
surf = read_csv("/vast/palmer/home.mccleary/jc3893/30x30/prediction-surface_ex.csv")
names(surf)
surf = st_as_sf(surf, coords=c("x","y"))
eviwin <- raster("gmted/EVI_Winter_11to2_resampled_NA30x30_cropped.tif") %>%
  projectRaster(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
ex <- extract(eviwin, surf)
surf <- data.frame(sf:::as_Spatial(surf))
surf = cbind(surf, ex)
names(surf)
colnames(surf)[ncol(surf)] = "eviwin"
surf$optional = NULL
names(surf)[58:59] <- c("x", "y")
write_csv(surf, "/vast/palmer/home.mccleary/jc3893/30x30/prediction-surface_ex.csv")

# add in SPEI layer
surf = read_csv("/gpfs/gibbs/pi/jetz/data/temporary_datasets/ebird_extremes_jeremy/prediction-surface_ex2.csv")
surf = st_as_sf(surf, coords=c("x","y"))
ex1 <- extract(spei_s, surf)
ex2 <- extract(spei_w, surf)
surf <- data.frame(sf:::as_Spatial(surf))
surf = cbind(surf, ex1, ex2)
names(surf)
colnames(surf)[62:63] = c("spei_s","spei_w")
surf$optional = NULL
names(surf)[59:60] <- c("x", "y")
write_csv(surf, "/gpfs/gibbs/pi/jetz/data/temporary_datasets/ebird_extremes_jeremy/prediction-surface_ex.csv")

{ # Optional plot
  data2 = data[complete.cases(data),]
  world <- ne_countries(scale='medium', returnclass = 'sf')
  datasp = st_as_sf(data2, coords=c("x","y"))
  lim <- st_bbox(datasp)
  # plot
  map <- ggplot() +
    geom_sf(data=world, bg="gray95") +
    coord_sf(xlim=lim[c(1,3)], ylim=lim[c(2,4)],
             expand=F) +
    theme_bw() +
    theme(text = element_text(size=13),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text.x = element_text(color = "black", size=13),
          axis.text.y = element_text(color = "black", size=13),
          plot.title = element_text(size = 22, face = "bold", hjust=.5)) +
    geom_point(data2, mapping=aes(x=x, y=y), col="black", size=.2) +
    xlab("") +
    ylab("")
}



# generate and add dummy layers with similar autocorrelation structure
# first load extreme layers above, then variability from bioclim
# crop to our extent, mask oceans
setwd("/gpfs/gibbs/pi/jetz/data/environmental_datasets/")
# frank extreme layers
ehe_s = raster("frank_extremes/prop-ehe-avg-bs.tif")
ece_s = raster("frank_extremes/prop-ece-avg-bs.tif")
s_stack = stack(ehe_s, ece_s) %>%
  projectRaster(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
ehe_w = raster("frank_extremes/prop-ehe-avg-bw.tif")
ece_w = raster("frank_extremes/prop-ece-avg-bw.tif")
w_stack = stack(ehe_w, ece_w) %>%
  projectRaster(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
# spei (already projected)
spei_s = raster("/vast/palmer/home.mccleary/jc3893/30x30/extremes/rasters/spei_summer.tif") 
spei_w = raster("/vast/palmer/home.mccleary/jc3893/30x30/extremes/rasters/spei_winter.tif") 
# seasonality from bioclim
bio4 = raster("WC/wc2.1_30s_bio_4.tif") %>%
  crop(c(-150, -40, 10, 60)) %>%
  projectRaster(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
bio15 = raster("WC/wc2.1_30s_bio_15.tif") %>%
  crop(c(-150, -40, 10, 60)) %>%
  projectRaster(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")


# load country polygons, transform, crop everything
countries = st_read("/vast/palmer/home.mccleary/jc3893/global/shapefiles", "countries") %>%
  st_transform(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
bbox = c(-17284324,-4691257,1399609,6500000)
s_stack = crop(s_stack, bbox) %>% mask(countries)
w_stack = crop(w_stack, bbox) %>% mask(countries)
frank_stack = stack(s_stack, w_stack)

spei_s = crop(spei_s, bbox) # oceans already masked
spei_w = crop(spei_w, bbox) # oceans already masked
bio4 = crop(bio4, bbox)
bio15 = crop(bio4, bbox)

# simulate layers with similar autocorrelation structure
# s = summer, w = winter
# ehe = extreme heat, ece = extreme cold
s_ehe = frank_stack[[1]]
s_ece = frank_stack[[2]]
w_ehe = frank_stack[[3]]
w_ece = frank_stack[[4]]
s_ehe_sim = simulate_data(s_ehe)
s_ece_sim = simulate_data(s_ece)
w_ehe_sim = simulate_data(w_ehe)
w_ece_sim = simulate_data(w_ece)
frank_sim_stack = stack(s_ehe_sim, s_ece_sim, w_ehe_sim, w_ece_sim)
writeRaster(frank_sim_stack, "frank_extremes/frank_sims.tif", overwrite=TRUE)

spei_s_sim = simulate_data(spei_s)
spei_w_sim = simulate_data(spei_w)
writeRaster(spei_s_sim, "/vast/palmer/home.mccleary/jc3893/30x30/extremes/rasters/spei_sim_summer.tif")
writeRaster(spei_w_sim, "/vast/palmer/home.mccleary/jc3893/30x30/extremes/rasters/spei_sim_winter.tif")

bio4_sim = simulate_data(bio4)
bio15_sim = simulate_data(bio15)
writeRaster(bio4_sim, "/vast/palmer/home.mccleary/jc3893/30x30/extremes/rasters/bio4_sim.tif")
writeRaster(bio15_sim, "/vast/palmer/home.mccleary/jc3893/30x30/extremes/rasters/bio15_sim.tif")

# reloading
frank_sim_stack = stack("frank_extremes/frank_sims.tif")
s_ehe_sim = frank_sim_stack[[1]]
s_ece_sim = frank_sim_stack[[2]]
w_ehe_sim = frank_sim_stack[[3]]
w_ece_sim = frank_sim_stack[[4]]
spei_s_sim = raster("/vast/palmer/home.mccleary/jc3893/30x30/extremes/rasters/spei_sim_summer.tif")
spei_w_sim = raster("/vast/palmer/home.mccleary/jc3893/30x30/extremes/rasters/spei_sim_winter.tif")
bio4_sim = raster("/vast/palmer/home.mccleary/jc3893/30x30/extremes/rasters/bio4_sim.tif")
bio15_sim = raster("/vast/palmer/home.mccleary/jc3893/30x30/extremes/rasters/bio15_sim.tif")
biostack = stack(bio4_sim, bio15_sim)

# extract to ebird
for (season in c("summer","winter")){
  data <- read_csv(paste0("/gpfs/gibbs/pi/jetz/data/temporary_datasets/ebird_extremes_jeremy/ebird_2024_",
                          season,"_extracted.csv")) %>%
    st_as_sf(coords=c("x","y"))
  # extract
  if(season=="summer"){
    frankex1 <- extract(s_ehe_sim, data)
    frankex2 <- extract(s_ece_sim, data)
    speiex <- extract(spei_s_sim, data)
  }else{
    frankex1 <- extract(w_ehe_sim, data)
    frankex2 <- extract(w_ece_sim, data)
    speiex <- extract(spei_w_sim, data)}
  bioex = extract(biostack, data)
  data <- data.frame(sf:::as_Spatial(data))
  data$optional <- NULL
  names(data)[752:753] <- c("x", "y")
  data = cbind(data, frankex1, frankex2, speiex, bioex)
  if(season=="summer"){names(data)[754:758] <- c("sim_ehe_s","sim_ece_s","sim_spei_s",
                                                 "sim_bio4","sim_bio15")
  }else{
    names(data)[754:758] <- c("sim_ehe_w","sim_ece_w","sim_spei_w",
                              "sim_bio4","sim_bio15")}
  write_csv(data, paste0("/gpfs/gibbs/pi/jetz/data/temporary_datasets/ebird_extremes_jeremy/ebird_2024_",
                         season,"_extracted.csv"))
  gc()
}

# extract to prediction surface
s_sim = stack(s_ehe_sim, s_ece_sim)
w_sim = stack(w_ehe_sim, w_ece_sim)
spei_sim = stack(spei_s_sim, spei_w_sim)
data <- read_csv(paste0("/gpfs/gibbs/pi/jetz/data/temporary_datasets/ebird_extremes_jeremy/prediction-surface_ex.csv"))
data = st_as_sf(data, coords=c("x","y"))
ex1 = raster::extract(s_sim, data)
ex2 = raster::extract(w_sim, data)
ex3 = raster::extract(spei_sim, data)
ex4 = raster::extract(biostack, data)
data <- data.frame(sf:::as_Spatial(data))
data2 = cbind(data, ex1, ex2, ex3, ex4)
colnames(data2)[61:68] = c("sim_ehe_s","sim_ece_s","sim_ehe_w",
                           "sim_ece_w","sim_spei_s","sim_spei_w",
                           "sim_bio4","sim_bio15")
colnames(data2)[69:70] <- c("x", "y")
data2$optional = NULL
write_csv(data2, "/gpfs/gibbs/pi/jetz/data/temporary_datasets/ebird_extremes_jeremy/prediction-surface_ex.csv")




























