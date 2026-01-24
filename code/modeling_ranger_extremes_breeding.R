# SDM workflow (summer)

{library(maps)
  library(rgdal)
  library(raster)
  library(dismo)
  library(data.table)	
  library(scales)
  library(ranger)
  library(scam)
  library(fields)  
  library(PresenceAbsence) 
  library(pdp)
  library(tidyverse)
  library(gbm)
  library(sf)
  library(rgeos)
  library(plyr)
  library(ff)
  library(ffbase)
  library(foreach)
  library(parallel)
  library(doParallel)
  library(ggpubr)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(blockCV)
  library(jpeg)
  # library(devtools)
  # install_github('mcooper/moranfast')
  library(moranfast)  
  # library(terra)
  
  sessionInfo()
  capabilities() 
}

###############################
#MODELS (RF)
###############################
{ version <- "VE8"
setwd("/vast/palmer/home.mccleary/jc3893/")

# Species names
list <- read_csv("30x30/species_list_2024.csv") %>%
  subset(hawaii=="no" & seabird=="no")

# shapefile
countries <- readOGR("30x30/shapefiles", "countries")
countries_proj <- spTransform(countries, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ")
aba <- countries[c(9,38),] %>%
  spTransform(countries, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>%
  st_as_sf()

{# grid for thinning
  bb <- bbox(countries_proj)
  cs <- c(5000,5000)
  cc <- bb[, 1] + (cs/2)  # cell offset
  cd <- ceiling(diff(t(bb))/cs)  # number of cells per direction
  grd <- GridTopology(cellcentre.offset=cc, cellsize=cs, cells.dim=cd)
  grd
  sp_grd <- SpatialGridDataFrame(grd,
                                 data=data.frame(id=1:prod(cd)),
                                 proj4string=CRS(proj4string(countries)))
  }

# define season - change as needed
season="summer"; ss="breeding"; months=6:8; mindate=153; maxdate=243; scode=2
#season <- "winter"; ss="nonbreeding"; months <- c(12,1,2); mindate=335; maxdate=59; scode=3

# ebird data
ebird_full <- read.csv.ffdf(file=paste0("/gpfs/gibbs/pi/jetz/data/temporary_datasets/ebird_extremes_jeremy/ebird_2024_",
                                        season,"_extracted.csv"),
                            colClasses=c(rep("factor",2), rep("numeric",4),
                                         rep("integer",689), rep("numeric",58)))

# prediction surface
pred_surface_full <- read.csv.ffdf(file="/gpfs/gibbs/pi/jetz/data/temporary_datasets/ebird_extremes_jeremy/prediction-surface_ex.csv",
                                   colClasses=rep("numeric",62))

# standard grid
bbox <- readOGR("30x30/shapefiles", "30x30BBox_NAmerica") %>%
  spTransform(countries, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs ") %>%
  st_as_sf()
grid <- raster("30x30/grids/Global_reference_raster_1km_CEA.tif") %>%
  crop(bbox)

# Predictors
preds <- c(
  "year","day_of_year","time_observations_started","duration_minutes",
  "effort_distance_km","number_observers","bio1","bio12",
  "tri","elev","pland_10_cropland_rainfed", 
  "pland_100_mosaic_tree_shrub","pland_11_cropland_rainfed", 
  "pland_110_mosaic_herbacious","pland_12_cropland_rainfed", 
  "pland_120_shrubland","pland_121_shrubland","pland_122_shrubland",
  "pland_130_grassland","pland_140_lichens_mosses","pland_150_sparse",
  "pland_152_sparse","pland_153_sparse","pland_160_flooded_freshwater",
  "pland_170_flooded_saltwater","pland_180_flooded_shrub", 
  "pland_190_urban","pland_20_cropland_irrigated","pland_200_barren",     
  "pland_201_barren","pland_202_barren","pland_210_water",
  "pland_220_ice","pland_30_mosaic_cropland","pland_40_mosaic_natural_veg",
  "pland_50_evergreen_broadleaf","pland_60_deciduous_broadleaf", 
  "pland_61_deciduous_broadleaf","pland_62_deciduous_broadleaf", 
  "pland_70_evergreen_needleleaf","pland_71_evergreen_needleleaf",
  "pland_72_evergreen_needleleaf","pland_80_deciduous_needleleaf",
  "pland_81_deciduous_needleleaf","pland_82_deciduous_needleleaf",
  "pland_90_mixed_forest", "evisum", # change!
  "bio4", "bio15",
  "ehe_s", "ece_s", "spei_s",
  "sim_bio4","sim_bio15",
  "sim_ehe_s","sim_ece_s","sim_spei_s") # change!
# subset for models
preds_m = preds[c(1:47,53:57)]; preds_m = gsub("\\_s$","",preds_m)
preds_v = preds[c(1:49,55:57)]; preds_v = gsub("\\_s$","",preds_v)
preds_e = preds[1:52]; preds_e = gsub("\\_s$","",preds_e)

# mol ranges
mol.range <- readOGR("30x30/range_polygons", "jetz_ranges_sub")

# fix stat ellipse
source("fix_stat_ellipse.R")

# directories
dir <- paste0("/gpfs/gibbs/pi/jetz/from_loomis/data/results_30x30/birds/",version,"_30x30/")
dir.create(dir)
dir.create(paste0(dir, "maps/"))
dir.create(paste0(dir, "niche/"))
dir.create(paste0(dir, "pdp/"))
dir.create(paste0(dir, "reports/"))
dir.create(paste0(dir, "pointmaps/"))
setwd(dir)

# parallel computing
#n.cores <- parallel::detectCores()
n.cores <- 3
my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

}

# ----------------------------------------------------------------------
# RANGER
# ----------------------------------------------------------------------		
# find a species
# which(list[[1]]=="Hooded Oriole")
# list[451,]$Common
{# loop species
  series <- sample(1:nrow(list))
  foreach (i=series) %dopar% { 
    # for(i in series){
    try({
      {library(maps)
        library(rgdal)
        library(raster)
        library(dismo)
        library(data.table)	
        library(scales)
        library(ranger)
        library(scam)
        library(fields)  
        library(PresenceAbsence) 
        library(pdp)
        library(tidyverse)
        library(gbm)
        library(sf)
        library(rgeos)
        library(plyr)
        library(ff)
        library(ffbase)
        library(foreach)
        library(parallel)
        library(doParallel)
        library(ggpubr)
        library(rnaturalearth)
        library(rnaturalearthdata)
        library(blockCV)
        library(jpeg)
        library(moranfast)
        #library(terra)
      }
      #removeTmpFiles(0)
      gc()
      rm(sprange)
      # Establish names
      common <- as.character(list[i,1])
      sci.name <- as.character(list[i,2])
      sci_name <- gsub(" ","_",sci.name)
      code <- as.character(list[i,5])   
      print(c(i, common, sci_name))
      # check if species has been done
      setwd(dir)
      if(!file.exists(paste0(sci_name,"_complete_",season,".csv"))){
        rm(sprange)
        # Expert range = modeling domain
        setwd("/vast/palmer/home.mccleary/jc3893/")
        try({
          total.range <- readOGR(unzip(
            paste0("30x30/range_polygons/",code,"-range-2020.gpkg.zip"),
            paste0(code,"-range-mr-2020.gpkg")), "range")
          # Get seasonal range of species
          sprange <- total.range[total.range$season_name==ss |
                                   total.range$season_name=="resident" ,]
        })
        if(!exists("sprange")){
          try({
            total.range <- readOGR(unzip(
              paste0("30x30/range_polygons/",code,"-range-2021.gpkg.zip"),
              paste0(code,"-range-mr-2021.gpkg")), "range")
            # Get seasonal range of species
            sprange <- total.range[total.range$season_name==ss |
                                     total.range$season_name=="resident" ,]
          })
        }
        if(!exists("sprange")){
          sprange <- mol.range[mol.range$sciname==sci.name & 
                                 (mol.range$season==scode |
                                    mol.range$season==1),]
          sprange$season <- season
        }
        if(nrow(sprange)==0){STOP}
        try(file.remove(paste0(code,"-range-mr-2020.gpkg")))
        {# Unproj, box range
        sprange <- spTransform(sprange, CRS="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
        sprangebuff <- buffer(sprange, 200000)
        box <- bbox(sprangebuff)
        # plot(sprangebuff, col="bisque")
        # Back to SP dataframe
        range.df <- data.frame(ID=1, row.names="buffer")
        sprangebuff <- SpatialPolygonsDataFrame(sprangebuff, range.df, match.ID=F)
        # ST range
        rangest <- st_as_sf(sprangebuff)
        ograngest = st_as_sf(sprange)
        # check that range is in our bbox
        sf::sf_use_s2(FALSE)
        intersect <- st_intersects(aba, rangest, sparse=T)}
        # check that range intersects US/Canada, otherwise abort
        if(sum(lengths(intersect))>0){
          print(paste("START", i, common, season))
          # directories
          setwd(dir)
          dir.create(paste0(sci_name,"/"))
          dir.create(paste0(sci_name,"/predictions/"))
          dir.create(paste0(sci_name,"/predictions/ROR/"))
          dir.create(paste0(sci_name,"/predictions/PA/"))
          dir.create(paste0(sci_name,"/points/"))
          dir.create(paste0(sci_name,"/accuracy/"))
          dir.create(paste0(sci_name,"/rangemap/"))
          
          # Save range
          writeOGR(sprangebuff, paste0(sci_name,"/rangemap"),
                   paste0(sci_name,"_",season), driver="ESRI Shapefile",
                   overwrite_layer=T)
          
          # ebird - get cropped rows, columns needed
          # filter out arctic and western alaska
          # for shorebirds, cut data to June, for others check if they are early migrants and cut August if so
          if (season=="summer" & list[i,"order"]=="Charadriiformes"){
            ebird_index <- ffwhich(ebird_full, 
                                   y<box[2,2] &
                                     y>box[2,1] &
                                     x<box[1,2] &
                                     x>box[1,1] &
                                     day_of_year<183 &
                                     y <= 6005522 &
                                     x >= -13508079)
          }else{
          if (season=="summer" & list[i,"aug_migration"]=="migrating"){
            ebird_index <- ffwhich(ebird_full, 
                                   y<box[2,2] &
                                     y>box[2,1] &
                                     x<box[1,2] &
                                     x>box[1,1] &
                                     day_of_year<213 &
                                     y <= 6005522 &
                                     x >= -13508079)
          }else{
            ebird_index <- ffwhich(ebird_full, 
                                   y<box[2,2] &
                                     y>box[2,1] &
                                     x<box[1,2] &
                                     x>box[1,1] &
                                     y <= 6005522 &
                                     x >= -13508079)
          }}
          ebird <- data.frame(ebird_full[ebird_index,
                                         c(sci_name,preds,"y","x",
                                           "sampling_event_identifier")])
          # NA removal (only a few rows)
          ebird <- ebird[complete.cases(ebird),]
          # Fix extreme variable names (leave duplicates matching pred surface)
          if(season=="summer"){
            colnames(ebird) = gsub("\\_s$","",colnames(ebird))
          }else{
            colnames(ebird) = gsub("\\_w$","",colnames(ebird))
          }
          # Fix extreme variable class
          ebird$sim_ehe = as.numeric(as.character(ebird$sim_ehe))
          ebird$sim_ece = as.numeric(as.character(ebird$sim_ece))
          
          {# Cut ebird to range
          print("processing ebird data")
          ebird <- st_as_sf(ebird, coords = c("x", "y"), 
                            crs = crs(rangest))
          ebird.int <- st_intersects(ebird, rangest, sparse=F)
          ebird <- ebird[as.vector(ebird.int),]
          #plot(ebird[sample.int(nrow(ebird), nrow(ebird)/20),], max.plot=1)
          ebird <- data.frame(sf:::as_Spatial(ebird))
          ebird$optional <- NULL
          names(ebird)[(ncol(ebird)-1):ncol(ebird)] <- 
            c("x", "y")
          rm(ebird.int)}
          
          # Check if there's enough data
          if(nrow(ebird)<50){
            writeLines("stop", paste0(sci_name,"_",common,"_failed_too_few_points_",season,".txt"))
            STOP
          }
          
          {# Column for presence/absence
          so = ebird[,gsub(" ",".",sci_name)]
          ebird$species_observed = ifelse(so >= 1, 1, 0)
          pres = ebird[ebird$species_observed==1,]
          abs = ebird[ebird$species_observed==0,]
          pres_st = st_as_sf(pres, coords = c("x", "y"), 
                             crs = crs(rangest))}
          
          { # Optional plot
            world <- ne_countries(scale='medium', returnclass = 'sf') %>%
              st_transform(crs="+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
            lim <- st_bbox(pres_st)            
            # plot
            sprange2 <- st_as_sf(sprange)
            map <- ggplot() +
              geom_sf(data=world, bg="gray95") + 
              geom_sf(data=sprange2, bg="bisque1") +
              coord_sf(xlim=lim[c(1,3)], ylim=lim[c(2,4)],
                       expand=F) +
              theme_bw() +
              theme(text = element_text(size=10),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    axis.line = element_line(color = "black"),
                    axis.text.x = element_text(color = "black", size=9),
                    axis.text.y = element_text(color = "black", size=9),
                    plot.title = element_text(size = 22, face = "bold", hjust=.5)) +
              geom_point(ebird, mapping=aes(x=x, y=y), col="gray50", size=.2) +
              geom_point(pres, mapping=aes(x=x, y=y), col="black", bg="white", size=.5) +
              xlab("") +
              ylab("")
            ggsave(paste0("pointmaps/",sci_name,"_point_map_with_abs_",season,".jpeg"), map, "jpeg",
                   height=7, width=7, units="in")}
          
          {# Cut prediction surface to range
          print("processing prediction surface")
          surface_index <- ffwhich(pred_surface_full, 
                                   y<box[2,2] &
                                     y>box[2,1] &
                                     x<box[1,2] &
                                     x>box[1,1] &
                                     y <= 6005522 &
                                     x >= -13508079)
          pred_surface <- data.frame(pred_surface_full[surface_index,])
          pred_surface <- st_as_sf(pred_surface, coords = c("x", "y"), 
                                   crs = crs(sprange))
          surface.int <- st_intersects(pred_surface, rangest, sparse=F)
          pred_surface <- pred_surface[as.vector(surface.int),]
          pred_surface <- data.frame(sf:::as_Spatial(pred_surface))
          pred_surface$optional <- NULL
          names(pred_surface)[(ncol(pred_surface)-1):ncol(pred_surface)] <- 
            c("x", "y")
          pred_surface$twi <- NULL
          pred_surface <- pred_surface[complete.cases(pred_surface),]
          rm(surface.int)
          gc()}
          
          # Check if there's enough surface
          if(nrow(pred_surface)<10){
            writeLines("stop", paste0(sci_name,"_",common,"_failed_too_few_cells_",season,".txt"))
            STOP
          }
          
          { # Populate with values for non-spatially varying covariates
            pred_surface$duration_minutes <- 60
            pred_surface$effort_distance_km <- 1
            pred_surface$number_observers <- 1
            pred_surface$year <- 2020
            pred_surface$optional <- NULL
            # fix bio1
            pred_surface$bio1 = pred_surface$bio1/10
            
            # Find max hr of day for spatial predictions
            hrs <- ebird[,c("time_observations_started","species_observed")] %>%
              group_by(round(time_observations_started)) %>%
              summarize_at(vars(species_observed), mean)
            pred_surface$time_observations_started <- as.numeric(
              hrs[which.max(hrs$species_observed),1])
            # set day for spatial predictions
            if(season=="winter"){
              pred_surface$day_of_year <- sample(
                c(1:maxdate,mindate:365), nrow(pred_surface), replace=T)
            }else{
              pred_surface$day_of_year <- sample(
                mindate:maxdate, nrow(pred_surface), replace=T)
            }
            # Fix extreme variable names 
            if(season=="summer"){
              colnames(pred_surface) = gsub("\\_s$","",colnames(pred_surface))
            }else{
              colnames(pred_surface) = gsub("\\_w$","",colnames(pred_surface))
            }
          }
          
          {# Niche plots
          pres.means = data.frame(x=mean(pres$bio1), y=mean(pres$ehe))
          abs.means = data.frame(x=mean(abs$bio1), y=mean(abs$ehe))
          gg1 <- ggplot(pres, mapping=aes(x = bio1, y = ehe)) +
            geom_point(abs, mapping=aes(x = bio1, y = ehe),
                       col="gray80", alpha=0.2) +
            stat_clip_ellipse_niche(abs, mapping=aes(x = bio1, y = ehe),
                                    geom="polygon", level=.9, fill="gray80",
                                    color="black", alpha=0.5, size=.8) +
            geom_point(pres, mapping=aes(x = bio1, y = ehe),
                       col="orange", alpha=0.2) +
            stat_clip_ellipse_niche(geom="polygon", level=.9, fill="orange",
                                    color="black", alpha=0.5, size=.8) +
            theme_bw() +
            theme(text = element_text(size=12),
                  panel.grid.major = element_line(),
                  panel.grid.minor = element_line(),
                  axis.line = element_line(color = "black"),
                  axis.text.x = element_text(color = "black", size=9),
                  axis.text.y = element_text(color = "black", size=9),
                  plot.title = element_text(size = 15, face = "bold", hjust=.5)) +
            geom_point(data=pres.means, aes(x, y), size=3, col="black") +
            geom_point(data=abs.means, aes(x, y), size=3, col="gray50") +
            ggtitle("Temperature - heat") +
            xlab("Mean Annual Temp. (°C)") +
            ylab("Extreme heat index")

          pres.means = data.frame(x=mean(pres$bio1), y=mean(pres$ece))
          abs.means = data.frame(x=mean(abs$bio1), y=mean(abs$ece))
          gg2 <- ggplot(pres, mapping=aes(x = bio1, y = ece)) +
            geom_point(abs, mapping=aes(x = bio1, y = ece),
                       col="gray80", alpha=0.2) +
            stat_clip_ellipse_niche(abs, mapping=aes(x = bio1, y = ece),
                                    geom="polygon", level=.9, fill="gray80",
                                    color="black", alpha=0.5, size=.8) +
            geom_point(pres, mapping=aes(x = bio1, y = ece),
                       col="dodgerblue", alpha=0.2) +
            stat_clip_ellipse_niche(geom="polygon", level=.9, fill="dodgerblue",
                                    color="black", alpha=0.5, size=.8) +
            theme_bw() +
            theme(text = element_text(size=12),
                  panel.grid.major = element_line(),
                  panel.grid.minor = element_line(),
                  axis.line = element_line(color = "black"),
                  axis.text.x = element_text(color = "black", size=9),
                  axis.text.y = element_text(color = "black", size=9),
                  plot.title = element_text(size = 15, face = "bold", hjust=.5)) +
            geom_point(data=pres.means, aes(x, y), size=3, col="black") +
            geom_point(data=abs.means, aes(x, y), size=3, col="gray50") +
            ggtitle("Temperature - cold") +
            xlab("Mean Annual Temp. (°C)") +
            ylab("Extreme cold index")

          pres.means = data.frame(x=mean(pres$bio12), y=mean(pres$spei))
          abs.means = data.frame(x=mean(abs$bio12), y=mean(abs$spei))
          gg3 <- ggplot(pres, aes(x = bio12, y = spei)) +
            geom_point(abs, mapping=aes(x = bio12, y = spei),
                       col="gray80", alpha=0.2) +
            stat_clip_ellipse_niche(abs, mapping=aes(x = bio12, y = spei),
                                    geom="polygon", level=.9, fill="gray80",
                                    color="black", alpha=0.5, linewidth=.8) +
            geom_point(pres, mapping=aes(x = bio12, y = spei),
                       col="plum1", alpha=0.2) +
            stat_clip_ellipse_niche(geom="polygon", level=.9, fill="plum1",
                                    color="black", alpha=0.5, linewidth=.8) +
            theme_bw() +
            theme(text = element_text(size=12),
                  panel.grid.major = element_line(),
                  panel.grid.minor = element_line(),
                  axis.line = element_line(color = "black"),
                  axis.text.x = element_text(color = "black", size=9),
                  axis.text.y = element_text(color = "black", size=9),
                  plot.title = element_text(size = 15, face = "bold", hjust=.5)) +
            geom_point(data=pres.means, aes(x, y), size=3, col="black") +
            geom_point(data=abs.means, aes(x, y), size=3, col="gray50") +
            ggtitle("Precipitation") +
            xlab("Mean Annual Precip. (mm)") +
            ylab("SPEI")
          gga_niche <- ggarrange(gg1, gg2, gg3, nrow=1, ncol=3)
          ggsave(paste0("niche/",sci_name,"_niche_plots_",season,".jpeg"), gga_niche, "jpeg",
                 height=4, width=12, units="in")
          }
          
          {### Ellipse slopes
            
            # Presence
            pres_st = data.frame(scale(pres[,c("bio1","bio12","ehe","ece","spei")]))
            s1p <- coef(lm(pres_st[,c("bio1","ehe")]))[2]
            s2p <- coef(lm(pres_st[,c("bio1","ece")]))[2]
            s3p <- coef(lm(pres_st[,c("bio12","spei")]))[2]
            # Absence
            abs_st = data.frame(scale(abs[,c("bio1","bio12","ehe","ece","spei")]))
            s1a <- coef(lm(abs_st[,c("bio1","ehe")]))[2]
            s2a <- coef(lm(abs_st[,c("bio1","ece")]))[2]
            s3a <- coef(lm(abs_st[,c("bio12","spei")]))[2]
            
            slopes = c(s1p, s2p, s3p, s1a, s2a, s3a)
            slope.names = c("slope_pres_bio1_ehe", "slope_pres_bio1_ece",
                           "slope_pres_bio12_spei", "slope_abs_bio1_ehe", 
                           "slope_abs_bio1_ece", "slope_abs_bio12_spei")
          }
          
          # Splits
          print("splitting and filtering data")
          
          ebirdpts <- st_as_sf(ebird, coords = c("x", "y"), 
                               crs = crs(rangest))
          {#training and testing samples
            sb <- cv_spatial(x = ebirdpts,
                             column = "species_observed",
                             rows_cols = c(8, 8),
                             k = 4,
                             hexagon = FALSE,
                             selection = "systematic")
          }
          # limit folds to those with sufficient train and test records
          whk = which(sb$records$test_1>10 & sb$records$train_1>5)
          # loop folds
          for (k in whk){
            print(k)
            train1 = ebird[sb$biomod_table[,k]==T,]
            test = ebird[sb$biomod_table[,k]==F,]
            test = test[complete.cases(test),]
            # random split OOB from train data
            index = runif(nrow(train1))
            train = train1[which(index>.25),]
            train = train[complete.cases(train),]
            train.oob = train1[which(index<=.25),]
            train.oob = train.oob[complete.cases(train.oob),]
            
            {# Spatiotemporal thinning
              train$id <- NULL
              coordinates(train) <- c("x", "y")
              crs(train) <- crs(sp_grd)
              over <- over(train, sp_grd)
              train <- cbind(data.frame(train), over)
              # get cell id and week number for each checklist
              train$optional <- NULL
              checklist_cell <- train %>% 
                mutate(cell = id,
                       week = ceiling(day_of_year/7))
              # sample one checklist per grid cell per week
              # sample detection/non-detection independently 
              train <- checklist_cell %>% 
                group_by(species_observed, year, week, cell) %>% 
                sample_n(size = 1) %>% 
                ungroup()
            }
            
            {# Assess spatial autocorrelation (sample 20k rows for speed)
              outtrain = c(); outtest = c()
              for (aa in 1:3){
                if(nrow(train>20000)){
                trainsa = sample_n(train, 20000)}else{
                  trainsa = train}
                if(nrow(test>20000)){
                  testsa = sample_n(test, 20000)}else{
                    testsa = test}
                autotrain = moranfast(trainsa$species_observed, 
                                 trainsa$x, trainsa$y)
                autotest = moranfast(testsa$species_observed, 
                                     testsa$x, testsa$y)
                outtrain = c(outtrain, autotrain$observed)
                outtest = c(outtest, autotest$observed)
                print(aa)
              }
              # compile summary stats
              autocorr = c(mean(outtrain), sd(outtrain), mean(outtest), sd(outtest))
              assign(paste0("ac.",k), autocorr)
            }
            
            {# Assemble TrainOOB
              train.oob$id <- NULL
              coordinates(train.oob) <- c("x","y")
              crs(train.oob) <- crs(sp_grd)
              over <- over(train.oob, sp_grd)
              train.oob <- cbind(data.frame(train.oob), over)
              train$optional <- NULL
              checklist_cell_oob <- train.oob %>% 
                mutate(cell = id,
                       week = ceiling(day_of_year/7))
              train.oob <- checklist_cell_oob %>% 
                group_by(species_observed, year, week, cell) %>% 
                sample_n(size = 1) %>% 
                ungroup()
            }
            
            {# remove oversampling for occurrence model, limit size
              train <- train[!duplicated(train$sampling_event_identifier), ]
              # Binary Occurrence Response, coded as numeric
              train$pres_abs <- as.numeric(train[[which(names(train) == "species_observed")]] > 0)
              # balanced sampling
              pos_fraction <- mean(as.numeric(train$pres_abs))
              # Check that positives ARE the minority class
              if (pos_fraction > 0.5) pos_fraction <- 1 - pos_fraction
              # Ranger binary response model requires code response as factor
              train$pres_abs <- as.factor(as.numeric(train$pres_abs))}
            
            # Save points alone
            if(k==max(whk)){
              write_csv(train[,c("y","x","pres_abs")], 
                        paste0(sci_name,"/points/pts_binary_occurrences_",
                               sci_name,"_",season,"_00.csv"))
            }
            
            # Loop different temp models
            for (temp in c("means","variability","extremes")){
              print(paste("model", temp))
              if(temp=="means"){pred_names=preds_m}
              if(temp=="variability"){pred_names=preds_v}
              if(temp=="extremes"){pred_names=preds_e}
              
              {#MODEL
                # Model Formula, threads
                m_formula <- paste("pres_abs ~", paste(pred_names, collapse = "+"))
                n_threads <- 7
                # Balanced Ranger Occurrence Model
                rf_occ <- ranger::ranger(
                  formula =  m_formula, 
                  num.trees = 100, 
                  #max.depth = 0, 
                  importance = "impurity",
                  num.threads = n_threads,
                  respect.unordered.factors = "order",
                  always.split.variables = NULL,
                  probability = TRUE,
                  replace = TRUE, 
                  sample.fraction = c(pos_fraction, pos_fraction),
                  data = train)
              }
              
              # Test Set PPMS
              test$occ_pred <- predict(
                rf_occ, data = test, 
                type = "response",
                num.threads = n_threads)$predictions[,2]
              
              # Partial dependence plots of extremes
              if(temp == "extremes" & k == max(whk)){
                vars = c("ehe","ece","spei")
                varnames = c("Extreme heat index", "Extreme cold index",
                             "SPEI")
                # loop three extreme weather layers, get partial dependence, plot
                for(u in 1:3){
                  par = pdp::partial(rf_occ, pred.var=vars[u],
                                     trim.outliers = TRUE, chull = TRUE, parallel = TRUE,
                                     grid.resolution = 30,  paropts = list(.packages = "ranger"))
                  names(par)[1] = "var"
                  ggpd = ggplot(data=par, aes(x=var, y=yhat)) +
                    geom_point() +
                    geom_smooth() +
                    ylab("Predicted occurrence") +
                    xlab(varnames[u]) +
                    theme_light() +
                    theme(legend.position = "none",
                          text = element_text(size = 12))
                  assign(paste0("ggpd.",vars[u]), ggpd)
                }
                gga_pd <- ggarrange(ggpd.ehe, ggpd.ece, ggpd.spei, nrow=1, ncol=3)
                ggsave(paste0("pdp/",sci_name,"_pdp_",season,"_",temp,".jpeg"), gga_pd, "jpeg",
                       height=4, width=11, units="in")
              }
              
              
              {#Select Optimal Threshold on train OOB data
                trainOOB.pred <- predict(rf_occ, data=train.oob, type="response")
                pa.df <- data.frame("space.holder.tag",
                                    obs = train.oob$species_observed, 
                                    ppp = trainOOB.pred$predictions[,2] )
                pa.metrics <- presence.absence.accuracy(
                  pa.df,
                  threshold = quantile(
                    pa.df$ppp,
                    probs = seq(from=0, to=1, length=1000), na.rm =T),
                  na.rm = T, st.dev = F)
                pa.metrics$PCC <- pa.metrics$Kappa <- NULL
                pa.metrics$TSS <- (pa.metrics$sensitivity + pa.metrics$specificity) - 1
                optimal_thresh_position <- which.max(pa.metrics$sensitivity +
                                                       pa.metrics$specificity)
                
                # Save confusion matrix
                if(k==max(whk)){
                  write_csv(pa.metrics, 
                            paste0(sci_name,"/accuracy/accuracy_matrix_",
                                   sci_name,"_",season,"_",temp,"_00.csv"))
                }
                
              }
              {# Compute Test set PPMs
                test.pred <- predict(rf_occ, data=test, type="response")
                pa.df <- data.frame("space.holder.tag",
                                    obs = test$species_observed, 
                                    ppp = test.pred$predictions[,2] )
                pa.metrics <- presence.absence.accuracy(
                  pa.df,
                  threshold = pa.metrics$threshold[ optimal_thresh_position ] ,
                  na.rm = T, st.dev = F)
                thresh <- pa.metrics$threshold
                pa.metrics["PCC"] <- pa.metrics["Kappa"] <- NULL
                pa.metrics["TSS"] <- (pa.metrics$sensitivity + pa.metrics$specificity) - 1
                
                assign(paste0("pa.metrics.",temp,".",k), pa.metrics)
              }
              
              # Spatial predictions
              pred_surface$occ_pred <- predict(
                rf_occ, data = pred_surface, type = "response",
                num.threads = n_threads)$predictions[,2]
              pred_surface$occ_thresh <- ifelse(
                pred_surface$occ_pred>pa.metrics$threshold,1,0)
              
              {# How much of range overlaps expert range?
                new_surface = pred_surface[pred_surface$occ_thresh==1,]
                new_surface <- st_as_sf(new_surface, coords = c("x", "y"), 
                                        crs = crs(sprange))
                new.surface.int <- st_intersects(new_surface, ograngest, sparse=F)
                new_surface2 <- new_surface[as.vector(new.surface.int),]
                prop = nrow(new_surface2)/nrow(new_surface)
                assign(paste0("prop.",temp,".",k), prop)
                # also save range area
                assign(paste0("range.",temp,".",k), nrow(new_surface))}
              
              if(k==max(whk)){
                {# Save rasterized predictions
                  coordinates(pred_surface) <- c("x","y")
                  ras <- rasterize(pred_surface, grid, pred_surface$occ_pred, fun='first')
                  writeRaster(ras, paste0(sci_name,"/predictions/ROR/PO_predictions_",
                                          sci_name,"_",season,"_",temp,"_00.tif"),
                              overwrite=T)
                  assign(paste0("ras.",temp), ras) # need for edge/core assessment
                  
                  pred_ones <- pred_surface[pred_surface$occ_thresh==1,]
                  ras <- rasterize(pred_ones, grid, pred_ones$occ_thresh, fun='first')
                  
                  writeRaster(ras, paste0(sci_name,"/predictions/PA/Thresholded_predictions_",
                                          sci_name,"_",season,"_",temp,"_00.tif"),
                              overwrite=T)
                  pred_surface <- data.frame(pred_surface)
                  pred_ones <- data.frame(pred_ones)
                  
                  occ.probs <- pred_surface[,c("x","y","occ_pred")]
                  write_csv(occ.probs, paste0(sci_name,"/predictions/ROR/PO_predictions_",
                                              sci_name,"_",season,"_",temp,"_00.csv"))
                  
                  occ.probs <- pred_ones[,c("x","y","occ_thresh")]
                  write_csv(occ.probs, paste0(sci_name,"/predictions/PA/Thresholded_predictions_",
                                              sci_name,"_",season,"_",temp,"_00.csv"))
                  rm(pred_ones)}
                
                
                {# Plot maps
                  jpeg(paste0("maps/PO_spatial_predictions_",
                              sci_name,"_",season,"_",temp,"_00.jpeg"),
                       width=500, height=500, units="px")
                  # Plot prob of occurrence
                  plotres <- 500
                  quilt.plot(
                    x = pred_surface$x,
                    y = pred_surface$y,
                    z = pred_surface$occ_pred,
                    nrow = plotres,
                    ncol = plotres,
                    na.rm = T,
                    main = paste(common,"-",sci.name,"\n",season,"-",temp))
                  plot(countries_proj, col="transparent", border="black",
                       usePolypath=F, add=T)
                  dev.off()
                  # save ror
                  ror = pred_surface$occ_pred
                  assign(paste0("ror.",temp), ror)
                  # Plot thresholded occurrence
                  jpeg(paste0("maps/Thresholded_spatial_predictions_",
                              sci_name,"_",season,"_",temp,"_00.jpeg"),
                       width=500, height=500, units="px")
                  quilt.plot(
                    x = pred_surface$x,
                    y = pred_surface$y,
                    z = pred_surface$occ_thresh,
                    nrow = plotres,
                    ncol = plotres,
                    na.rm = T,
                    add.legend = F,
                    col = c("transparent", "darkgreen"),
                    main = paste(common,"-",sci.name,"\n",season,"-",temp))
                  plot(countries_proj, col="transparent", border="black",
                       usePolypath=F, add=T)
                  dev.off()
                  # save occ thresh
                  ot = pred_surface$occ_thresh
                  assign(paste0("ot.",temp), ot)
                  } # end plot
                
                # Averaging across folds
                for (w in whk){
                  pa.metrics.w = get(paste0("pa.metrics.",temp,".",w))[2:6]
                  if(w==min(whk)){pa.metrics.table=pa.metrics.w}else{
                    pa.metrics.table=rbind(pa.metrics.table, pa.metrics.w)}
                  
                  ac.w = get(paste0("ac.",w))
                  if(w==min(whk)){ac.table=ac.w}else{
                    ac.table=rbind(ac.table, ac.w)}
                  
                  prop.w = get(paste0("prop.",temp,".",w))
                  if(w==min(whk)){props=prop.w}else{props=c(props, prop.w)}
                  
                  range.w = get(paste0("range.",temp,".",w))
                  if(w==min(whk)){ras=range.w}else{ras=c(ras, range.w)}
                }
                pa.metrics.summary = colMeans(pa.metrics.table, na.rm=T)
                ac.summary = colMeans(ac.table, na.rm=T)
                prop.summary = mean(props, na.rm=T)
                range.summary = mean(ras, na.rm=T)
                
                {# Species values
                  modName <- paste0(sci_name,"_",season,"_",temp)
                  imp.scores <- rf_occ$variable.importance
                  # if(temp=="means"){imp.scores = c(imp.scores, rep(NA, 5))} # 7 originally
                  # if(temp=="variability"){imp.scores = c(imp.scores, rep(NA, 3))}
                  model.output <- unlist(c(common, sci_name, season, temp, version, 
                                           modName, nrow(train), range.summary, 
                                           pa.metrics.summary[4],
                                           pa.metrics.summary[c(1,5,2,3)],
                                           prop.summary, ac.summary, slopes, imp.scores))
                  names(model.output)<- c("common_name","sciname", 
                                          "season","temperature","version","modName",
                                          "noPts","range_area",
                                          "AUC","SPSthreshold","TSS",
                                          "Sensitivity","Specificity","prop_in_range",
                                          "train.ac.mean","train.sc.sd",
                                          "test.ac.mean","test.ac.sd", slope.names,
                                          paste0(preds_e,"_importance"))
                  if(temp=="means"){sp.output <- model.output
                  }else{sp.output <- data.frame(rbind(sp.output, model.output))}}
              } # end k=max(whk) only 
            } # end model loop
          } # end k loop
          
          {# pdf output
          jpeg(paste0("maps/",sci_name,"_",season,"_predictions.jpeg"),
               width=800, height=1000)
          par(mfrow=c(2,3))
          quilt.plot(
            x = pred_surface$x,
            y = pred_surface$y,
            z = ot.means,
            nrow = plotres,
            ncol = plotres,
            na.rm = T,
            add.legend = F,
            col = c("transparent", "darkgreen"),
            main = paste(common,"-",sci.name,"\n",season,"- Means"))
          plot(countries_proj, col="transparent", border="black",
               usePolypath=F, add=T)
          quilt.plot(
            x = pred_surface$x,
            y = pred_surface$y,
            z = ot.variability,
            nrow = plotres,
            ncol = plotres,
            na.rm = T,
            add.legend = F,
            col = c("transparent", "darkgreen"),
            main = paste(common,"-",sci.name,"\n",season,"- Variability"))
          plot(countries_proj, col="transparent", border="black",
               usePolypath=F, add=T)
          quilt.plot(
            x = pred_surface$x,
            y = pred_surface$y,
            z = ot.extremes,
            nrow = plotres,
            ncol = plotres,
            na.rm = T,
            add.legend = F,
            col = c("transparent", "darkgreen"),
            main = paste(common,"-",sci.name,"\n",season,"- Extremes"))
          plot(countries_proj, col="transparent", border="black",
               usePolypath=F, add=T)
          quilt.plot(
            x = pred_surface$x,
            y = pred_surface$y,
            z = ror.means,
            nrow = plotres,
            ncol = plotres,
            na.rm = T,
            main = paste(common,"-",sci.name,"\n",season,"-",temp))
          plot(countries_proj, col="transparent", border="black",
               usePolypath=F, add=T)
          quilt.plot(
            x = pred_surface$x,
            y = pred_surface$y,
            z = ror.variability,
            nrow = plotres,
            ncol = plotres,
            na.rm = T,
            main = paste(common,"-",sci.name,"\n",season,"-",temp))
          plot(countries_proj, col="transparent", border="black",
               usePolypath=F, add=T)
          quilt.plot(
            x = pred_surface$x,
            y = pred_surface$y,
            z = ror.extremes,
            nrow = plotres,
            ncol = plotres,
            na.rm = T,
            main = paste(common,"-",sci.name,"\n",season,"-",temp))
          plot(countries_proj, col="transparent", border="black",
               usePolypath=F, add=T)
          dev.off()

          ggblank = ggplot() + theme_void()
          ggarr = ggarrange(map, ggblank, ggblank, gg1, gg2, gg3,
                            ggpd.ehe, ggpd.ece, ggpd.spei,
                            nrow=3, ncol=3)
          ggsave(paste0("reports/",sci_name,"_",season,"_ggplots.jpeg"), ggarr,
                 width=10, height=6, units="in")

          jpeg(paste0("reports/",sci_name,"_",season,".jpeg"), width=1600, height=1600)
          par(mar=c(0,0,0,0), mfrow=c(2,1))
          im <- readJPEG(paste0("maps/",sci_name,"_",season,"_predictions.jpeg"))
          plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
          rasterImage(im,-0.03,-0.03,1.03,1.03)
          im <- readJPEG(paste0("reports/",sci_name,"_",season,"_ggplots.jpeg"))
          plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
          rasterImage(im,-0.03,-0.03,1.03,1.03)
          dev.off()

          pdf(paste0("reports/",sci_name,"_",season,".pdf"), width=25, height=25)
          par(mar=c(0,0,0,0), mfrow=c(2,1))
          im <- readJPEG(paste0("maps/",sci_name,"_",season,"_predictions.jpeg"))
          plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
          rasterImage(im,-0.03,-0.03,1.03,1.03)
          im <- readJPEG(paste0("reports/",sci_name,"_",season,"_ggplots.jpeg"))
          plot(0:1,0:1,type="n",ann=FALSE,axes=FALSE)
          rasterImage(im,-0.03,-0.03,1.03,1.03)
          dev.off()
          } # end pdf output
          
          # csv output
          write_csv(data.frame(sp.output), 
                    paste0(sci_name,"_complete_",season,".csv"))
        } # end range exists in bbox check
      } # end file exists check
    }) # end try
  }#end sp
} #end workflow



