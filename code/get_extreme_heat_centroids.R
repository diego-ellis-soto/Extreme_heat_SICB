##################################################################
## 1) Load Packages
##################################################################
library(daymetr)
library(dplyr)
library(lubridate)
library(heatwaveR)
library(ggplot2)
library(tidyr)
library(readr)
library(rnaturalearth)  # for world map
library(sf)             # for reading/wrangling shapefiles

valid_centroids_df = data.frame( # Washington ####
                                 centroid_lon = -117.5,
                                 centroid_lat = 48.3,
                                 indivID = 'Washington'
) %>%
  filter(!is.na(centroid_lat), !is.na(centroid_lon))

##################################################################
## 5) Parameters & Output Folder
##################################################################
hist_start   <- 1991
hist_end     <- 2020
compare_year <- 2021

pctile       <- 90  # Important

# Create a main output directory
output_dir <- "outdir/"
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

##################################################################
## 6) Helper Function for Daymet Download (June–July)
##################################################################
get_daymet_jj <- function(year, lat, lon, site_prefix) {
  # Attempt to download from Daymet, handle potential errors/warnings
  tmp <- tryCatch({
    download_daymet(
      site    = paste0(site_prefix, "_", year),
      lat     = lat,
      lon     = lon,
      start   = year,
      end     = year,
      internal = TRUE,
      force   = TRUE   # set FALSE if you want to skip re-downloading
    )
  }, error = function(e) {
    warning(paste("Error downloading DAYMET data for:", site_prefix, year, 
                  "at", lat, "/", lon, "->", e$message))
    return(NULL)  # Return NULL on error
  })
  
  # If download failed, return empty
  if (is.null(tmp)) {
    return(NULL)
  }
  
  # Otherwise subset to June-July
  df <- tmp$data %>%
    mutate(Date = as.Date(paste(year, yday), format = "%Y %j")) %>%
    filter(
      Date >= as.Date(paste0(year, "-06-01")),
      Date <= as.Date(paste0(year, "-07-31"))
    ) %>%
    select(Date, tmax_c = tmax..deg.c.)
  
  return(df)
}

##################################################################
## 7) Prepare a World Map
##    We'll use rnaturalearth for a simple background
##################################################################
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

##################################################################
## 8) Main Loop Over Each Valid Centroid
##    - For each (indivID, species, lat, lon):
##       (A) Download Daymet Data
##       (B) Run heatwaveR
##       (C) Reproduce your script's plots
##       (D) Make a map of June–July 2021 points + centroid
##       (E) Save results
##    - Organize outputs in subfolders by species
##################################################################

all_events_list       <- list()
all_climatology_list  <- list()

# Potential for scalability here and assess multiple centroids ####
for(i in seq_len(nrow(valid_centroids_df))) {
  # i = 1
  # Extract
  this_id  <- valid_centroids_df$indivID[i]
  lat_c    <- valid_centroids_df$centroid_lat[i]
  lon_c    <- valid_centroids_df$centroid_lon[i]
  
  message(sprintf(
    "Assessing extreme heat in  %s: lat=%.5f, lon=%.5f (%d/%d)",
    this_id, lat_c, lon_c, i, nrow(valid_centroids_df)
  ))
  
  # (A) Download historical (1991–2020) data
  df_hist_list <- list()
  
  for (yr in hist_start:hist_end) {
    df_temp <- get_daymet_jj(year = yr, lat = lat_c, lon = lon_c,
                             site_prefix = paste0("hist_", this_id))
    # If we got NULL, skip
    if(!is.null(df_temp)) {
      df_hist_list[[as.character(yr)]] <- df_temp
    }
  }
  # Combine historical
  if(length(df_hist_list) == 0) {
    warning(paste("No historical data for indivID=", this_id, "; skipping."))
    next
  }
  
  df_hist <- bind_rows(df_hist_list, .id = "year") %>%
    mutate(year = as.integer(year)) %>%
    arrange(Date)
  
  # (A2) Download compare_year
  df_2021 <- get_daymet_jj(compare_year, lat_c, lon_c,
                           site_prefix = paste0("compare_", this_id))
  # If missing, skip
  if(is.null(df_2021)) {
    warning(paste("No 2021 data for indivID=", this_id, "; skipping."))
    next
  }
  df_2021 <- df_2021 %>%
    mutate(year = compare_year)
  
  # Combine
  df_combined <- bind_rows(df_hist, df_2021) %>%
    select(Date, year, tmax_c) %>%
    rename(t = Date, temp = tmax_c) %>%
    arrange(t)
  
  # (B) Run heatwaveR
  df_hw_clim <- ts2clm(
    data               = df_combined,
    x                  = t,
    y                  = temp,
    climatologyPeriod  = c("1991-06-01", "2020-07-31"),
    pctile             = pctile
  )
  hw_events <- detect_event(df_hw_clim, x = t, y = temp)
  
  # Mark which individual & species these belong to
  event_df <- hw_events$event %>% 
    mutate(indivID = this_id)
  clim_df  <- hw_events$climatology %>%
    mutate(indivID = this_id)
  
  # Collect in lists for final combination
  all_events_list[[i]]      <- event_df
  all_climatology_list[[i]] <- clim_df
  
  # Build a prefix for file naming
  file_prefix <- this_id
  
  # ---------------------------
  # (C2) Event Line (if events)
  # ---------------------------
  has_events <- nrow(event_df) > 0
  
  if (has_events) {
    p_event_line <- tryCatch({
      event_line(hw_events, spread = 10)
    }, error = function(e) {
      warning(paste("event_line failed for", this_id, "->", e$message))
      return(NULL)
    })
    if(!is.null(p_event_line)) {
      ggsave(
        filename = file.path(output_dir, paste0(file_prefix, "_event_line.png")),
        plot     = p_event_line, width = 7, height = 5
      )
    }
  } else {
    message("No events found for indivID=", this_id, " -> skipping event_line")
  }
  
  # ---------------------------
  # (C4) Consecutive Exceedance Days in 2021
  # ---------------------------
  exceed_days_2021 <- clim_df %>%
    filter(
      format(t, "%Y") == "2021",
      temp > thresh
    ) %>%
    select(t, temp, thresh, indivID)
  
  exceed_days_2021 <- exceed_days_2021 %>%
    arrange(t) %>%
    mutate(
      day_diff  = as.integer(difftime(t, lag(t), units = "days")),
      new_group = ifelse(is.na(day_diff) | day_diff != 1, 1, 0),
      group     = cumsum(new_group)
    )
  
  heatwave_intervals_2021 <- exceed_days_2021 %>%
    group_by(group) %>%
    summarize(
      start_date    = first(t),
      end_date      = last(t),
      duration_days = as.integer(difftime(last(t), first(t), units = "days")) + 1,
      .groups       = 'drop'
    ) %>%
    mutate(indivID = this_id)
  
  
  
  # ---------------------------
  # (D) Make a Map for This Individual's June–July 2021 Points
  #     and the Centroid
  # ---------------------------
  
  # centroid
  
  usa_states <- ne_states(country = "United States of America", returnclass = "sf")
  wa <- usa_states %>% filter(name == "Washington") %>% st_transform(4326)
  
  pts_sf  <- st_as_sf(valid_centroids_df, coords = c("centroid_lon", "centroid_lat"), crs = 4326)
  cent_sf <- st_as_sf(valid_centroids_df, coords = c("centroid_lon", "centroid_lat"), crs = 4326)
  
  # A little padding around WA for nicer framing
  bb <- st_bbox(wa)
  pad_x <- (bb$xmax - bb$xmin) * 0.05
  pad_y <- (bb$ymax - bb$ymin) * 0.05
  xlim <- c(bb$xmin - pad_x, bb$xmax + pad_x)
  ylim <- c(bb$ymin - pad_y, bb$ymax + pad_y)
  
  ggplot() +
    # annotation_map_tile(type = "stamen_terrain", zoomin = 0) +   # terrain-like topography
    geom_sf(data = wa, fill = NA, color = "white", linewidth = 0.8) +
    geom_sf(data = pts_sf, alpha = 0.35, size = 2.0, color = "black") +
    geom_sf(data = cent_sf, size = 4.2, shape = 21, fill = "red", color = "white", stroke = 0.8) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(
      title = "June–July 2021 Locations & Centroid",
      subtitle = "Washington State (terrain basemap)",
      x = NULL, y = NULL
    ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey20")
    )
  
  
  
  if(nrow(valid_centroids_df) > 0) {
    # Build bounding box around the data
    extent_x <- range(c(valid_centroids_df$lon, lon_c), na.rm = TRUE)
    extent_y <- range(c(valid_centroids_df$lat, lat_c), na.rm = TRUE)
    # Add a small margin (in degrees) so we don't cut off points
    margin <- 1
    extent_x[1] <- extent_x[1] - margin
    extent_x[2] <- extent_x[2] + margin
    extent_y[1] <- extent_y[1] - margin
    extent_y[2] <- extent_y[2] + margin
    
    # Plot
    p_map <- ggplot() +
      # World background
      geom_sf(data = world, fill = "lightgray", color = "white") +
      # Individual's points
      geom_point(data = valid_centroids_df, 
                 aes(x = centroid_lon, y = centroid_lat),
                 color = "darkgreen", alpha = 0.5, size = 2) +
      # Centroid
      geom_point(data = data.frame(cx = lon_c, cy = lat_c),
                 aes(x = cx, y = cy),
                 color = "red", size = 3) +
      labs(
        title    = paste0("indivID=", this_id, " (", this_sp, ")"),
        subtitle = "June–July 2021 Locations & Centroid",
        x        = "Longitude",
        y        = "Latitude"
      ) +
      coord_sf(xlim = extent_x, ylim = extent_y) +
      theme_minimal()
    
    # Save map
    ggsave(
      filename = file.path(output_dir, paste0(file_prefix, "_map.png")),
      plot     = p_map, width = 7, height = 5
    )
  }
  
  # (E) (Optional) Save CSVs if desired
  # write.csv(event_df,  file.path(output_dir, paste0(file_prefix, "_events.csv")),        row.names = FALSE)
  # write.csv(clim_df,   file.path(output_dir, paste0(file_prefix, "_climatology.csv")),   row.names = FALSE)
  # write.csv(exceed_days_2021,   file.path(output_dir, paste0(file_prefix, "_exceed_days.csv")), row.names = FALSE)
  # write.csv(heatwave_intervals_2021, file.path(output_dir, paste0(file_prefix, "_intervals.csv")), row.names = FALSE)
}

##################################################################
## 9) Combine All Individuals' Event/Clim Data
##################################################################
all_events_df       <- bind_rows(all_events_list)
all_climatology_df  <- bind_rows(all_climatology_list)

# Inspect
head(all_events_df)
head(all_climatology_df)

# Summaries
all_events_df %>%
  group_by(indivID, species) %>%
  summarize(n_events = n())

##################################################################
## 10) Create a CSV of Heatwave Periods in 2020 & 2021
##     for Subsequent Analysis
##################################################################
# We'll filter events that START/PEAK/END in 2020 or 2021
hw_2020_2021 <- all_events_df %>%
  filter(
    year(date_start) %in% c(2020, 2021) |
      year(date_peak) %in% c(2020, 2021) |
      year(date_end)  %in% c(2020, 2021)
  )

# Save to CSV if you wish
write.csv(hw_2020_2021, file.path(output_dir, "extreme_heat_example_case_figure_1.csv"), row.names = FALSE)

##################################################################
## 11) (Optional) Save the Combined Event/Climatology Data
##################################################################
# write.csv(all_events_df,        file.path(output_dir, "all_events.csv"), row.names = FALSE)
# write.csv(all_climatology_df,   file.path(output_dir, "all_climatology.csv"), row.names = FALSE)

##################################################################
## The folder 'heatwave_analysis_output' now contains subfolders 
## for each species. Inside each species folder, you’ll find:
##  - "ts_plot.png" and "ts_shading.png"
##  - "map.png" showing June–July 2021 points and centroid
##  - (Optional) CSV files if you uncomment them.
##
## A CSV called 'heatwaves_2020_2021.csv' in 'output_dir' 
## summarizes the time periods where a heatwave was detected 
## overlapping 2020 or 2021.
##################################################################  the green and the red locations are too far away


require(mapview)
mapview(cent_sf, col.regions = "red", cex = 6)
