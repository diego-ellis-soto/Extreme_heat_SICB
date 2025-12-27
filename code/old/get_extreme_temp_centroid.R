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

## For world map
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)  # optional, if you want the sf approach

##################################################################
## 2) Read Input Data
##    - Must have columns:
##      indivID, species, timestamp, lat, lon, etc.
##################################################################
animal_df <- read_csv('/Users/diegoellis/Downloads/all_tracking_data_resampled_static_anno_timestamp_for_GEE__mosey_anno.csv')

##################################################################
## 3) Calculate June–July 2021 Centroids per Individual
##    (Also keep species info)
##################################################################
# Initialize an empty list to store results
centroid_list <- list()

# Get unique individual IDs
unique_indiv_ids <- unique(animal_df$indivID)

for (indiv_id in unique_indiv_ids) {
  # Subset data for this individual
  df_sub <- animal_df %>%
    filter(indivID == indiv_id)
  
  # Identify species (assuming the same species for all rows of this individual)
  this_species <- unique(df_sub$species)
  if (length(this_species) > 1) {
    warning(paste("Multiple species found for indivID =", indiv_id, 
                  "- taking first. Values:", paste(this_species, collapse=", ")))
    this_species <- this_species[1]
  }
  
  # Filter to June–July 2021
  df_jj2021 <- df_sub %>%
    filter(
      year(timestamp) == 2021,
      month(timestamp) %in% c(6, 7)
    )
  
  # Calculate mean lat and lon
  centroid_lat <- mean(df_jj2021$lat, na.rm = TRUE)
  centroid_lon <- mean(df_jj2021$lon, na.rm = TRUE)
  
  # Store results in the list
  centroid_list[[indiv_id]] <- data.frame(
    indivID       = indiv_id,
    species       = this_species,
    centroid_lat  = centroid_lat,
    centroid_lon  = centroid_lon
  )
}

# Combine into one data frame
centroids_df <- bind_rows(centroid_list)
print(centroids_df)

##################################################################
## 4) Filter Out NA Lat/Lon Rows
##################################################################
valid_centroids_df <- centroids_df %>%
  filter(!is.na(centroid_lat), !is.na(centroid_lon))

##################################################################
## 5) Parameters & Output Folder
##################################################################
hist_start   <- 1991
hist_end     <- 2020
compare_year <- 2021

pctile       <- 95  # or 90, as you prefer

# Create a main output directory
output_dir <- "heatwave_analysis_output"
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

##################################################################
## 5a) Prepare the World Map (rnaturalearth)
##################################################################
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

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
## 7) Main Loop Over Each Valid Centroid
##    - For each (indivID, species, lat, lon):
##       (A) Download Daymet Data
##       (B) Run heatwaveR
##       (C) Reproduce your script's plots
##       (D) Save results
##    - Organize outputs in subfolders by species
##    - Also create a MAP showing all June–July 2021 points + centroid
##################################################################

all_events_list       <- list()
all_climatology_list  <- list()

for(i in seq_len(nrow(valid_centroids_df))) {
  
  # Extract
  this_id  <- valid_centroids_df$indivID[i]
  this_sp  <- valid_centroids_df$species[i]
  lat      <- valid_centroids_df$centroid_lat[i]
  lon      <- valid_centroids_df$centroid_lon[i]
  
  # Optional progress message
  message(sprintf("Processing: indivID=%s, species=%s, lat=%.5f, lon=%.5f (%d/%d)",
                  this_id, this_sp, lat, lon, i, nrow(valid_centroids_df)))
  
  # Create a species-specific subfolder
  species_dir <- file.path(output_dir, gsub(" ", "_", this_sp))  # replace spaces with underscores
  if(!dir.exists(species_dir)) {
    dir.create(species_dir)
  }
  
  # (i0) For the MAP: 
  #      Subset the individual's data for June–July 2021
  #      We'll overlay these points on the map
  df_sub <- animal_df %>%
    filter(indivID == this_id) %>%
    filter(
      year(timestamp) == 2021,
      month(timestamp) %in% c(6, 7)
    )
  
  # (A) Download historical (1991–2020) data
  df_hist_list <- list()
  for (yr in hist_start:hist_end) {
    df_temp <- get_daymet_jj(year = yr, lat = lat, lon = lon,
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
  df_2021 <- get_daymet_jj(compare_year, lat, lon,
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
    mutate(indivID = this_id, species = this_sp)
  clim_df  <- hw_events$climatology %>%
    mutate(indivID = this_id, species = this_sp)
  
  # Collect in lists for final combination
  all_events_list[[i]]      <- event_df
  all_climatology_list[[i]] <- clim_df
  
  # Build a prefix for file naming
  file_prefix <- paste0("indiv_", this_id, "_sp_", gsub(" ", "_", this_sp))
  
  # ---------------------------
  # (C0) MAP: Points + Centroid
  # ---------------------------
  if(nrow(df_sub) > 0) {
    # We have at least some points
    p_map <- ggplot() +
      geom_sf(data = world, fill = "antiquewhite") +
      # individual's points
      geom_point(
        data = df_sub,
        aes(x = lon, y = lat),
        color = "darkred", alpha = 0.5, size = 2
      ) +
      # centroid
      geom_point(
        data = data.frame(x = lon, y = lat),
        aes(x = x, y = y),
        color = "blue", size = 3
      ) +
      coord_sf(xlim = c(min(df_sub$lon, lon) - 5, max(df_sub$lon, lon) + 5),
               ylim = c(min(df_sub$lat, lat) - 5, max(df_sub$lat, lat) + 5),
               expand = FALSE) +
      labs(
        title = paste0("Map of indivID=", this_id, " (", this_sp, ")"),
        subtitle = "June–July 2021 raw points + centroid",
        x = "Longitude",
        y = "Latitude"
      ) +
      theme_minimal()
    
    ggsave(
      filename = file.path(species_dir, paste0(file_prefix, "_location_map.png")),
      plot     = p_map,
      width    = 7, 
      height   = 5
    )
  } else {
    message("No June–July 2021 points for indivID=", this_id, " -> skipping map")
  }
  
  # ---------------------------
  # (C1) Reproduce the Plots
  # ---------------------------
  
  # Check if there are any events
  has_events <- nrow(event_df) > 0
  
  ## 1. Lolli plot (only if there's at least 1 event)
  if (has_events) {
    p_lolli <- tryCatch({
      lolli_plot(hw_events, metric = "intensity_max")
    }, error = function(e) {
      warning(paste("lolli_plot failed for", this_id, "->", e$message))
      return(NULL)
    })
    if(!is.null(p_lolli)) {
      ggsave(
        filename = file.path(species_dir, paste0(file_prefix, "_lolli_plot.png")),
        plot     = p_lolli, width = 7, height = 5
      )
    }
  } else {
    message("No events found for indivID=", this_id, " -> skipping lolli_plot")
  }
  
  ## 2. Event line (spread=10) (only if events)
  if (has_events) {
    p_event_line <- tryCatch({
      event_line(hw_events, spread = 10)
    }, error = function(e) {
      warning(paste("event_line failed for", this_id, "->", e$message))
      return(NULL)
    })
    if(!is.null(p_event_line)) {
      ggsave(
        filename = file.path(species_dir, paste0(file_prefix, "_event_line.png")),
        plot     = p_event_line, width = 7, height = 5
      )
    }
  } else {
    message("No events found for indivID=", this_id, " -> skipping event_line")
  }
  
  ## 3. Filter data to just 2020–2021 for simpler TMAX vs threshold plot
  df_2021_clim <- hw_events$climatology %>%
    filter(format(t, "%Y") %in% c('2020', '2021'))
  
  p_ts <- ggplot(df_2021_clim, aes(x = t)) +
    geom_line(aes(y = temp, color = "TMAX (°C)")) +
    geom_line(aes(y = thresh, color = "Threshold"), linetype = "dashed") +
    labs(
      title = paste0("indivID=", this_id, " (", this_sp, 
                     "): Jun–Jul 2020-2021 TMAX vs. ", pctile, "th"),
      x = "Date",
      y = "Temperature (°C)",
      color = "Legend"
    ) +
    scale_color_manual(values = c("TMAX (°C)" = "firebrick",
                                  "Threshold" = "darkblue")) +
    theme_minimal()
  
  ggsave(
    filename = file.path(species_dir, paste0(file_prefix, "_ts_plot.png")),
    plot     = p_ts, width = 7, height = 5
  )
  
  ## 4. Days in 2021 Exceeding the Threshold (Group consecutive days)
  #    We can do this even if events are zero, might result in empty table
  exceed_days_2021 <- clim_df %>%
    filter(
      format(t, "%Y") == "2021",
      temp > thresh
    ) %>%
    select(t, temp, thresh, indivID, species)
  
  exceed_days_2021 <- exceed_days_2021 %>%
    arrange(t) %>%
    mutate(
      day_diff = as.integer(difftime(t, lag(t), units = "days")),
      new_group = ifelse(is.na(day_diff) | day_diff != 1, 1, 0),
      group = cumsum(new_group)
    )
  heatwave_intervals_2021 <- exceed_days_2021 %>%
    group_by(group) %>%
    summarize(
      start_date    = first(t),
      end_date      = last(t),
      duration_days = as.integer(difftime(last(t), first(t), units = "days")) + 1,
      .groups       = 'drop'
    ) %>%
    mutate(indivID = this_id, species = this_sp)
  
  ## 5. Plot TMAX vs threshold for 2021 only with shading
  df_2021_clim_strict <- clim_df %>%
    filter(format(t, "%Y") == "2021")
  
  p_shading <- ggplot(df_2021_clim_strict, aes(x = t)) +
    geom_line(aes(y = temp, color = "TMAX (°C)")) +
    geom_line(aes(y = thresh, color = "Threshold"), linetype = "dashed") +
    labs(
      title    = paste0("indivID=", this_id, " (", this_sp, 
                        "): June–July 2021 TMAX vs. ", pctile, "th"),
      x        = "Date",
      y        = "Temperature (°C)",
      color    = "Legend"
    ) +
    scale_color_manual(values = c("TMAX (°C)" = "firebrick",
                                  "Threshold" = "darkblue")) +
    theme_minimal()
  
  if(nrow(heatwave_intervals_2021) > 0) {
    p_shading <- p_shading +
      geom_rect(
        data = heatwave_intervals_2021,
        aes(xmin = start_date, xmax = end_date, ymin = -Inf, ymax = Inf),
        fill = "navyblue",
        alpha = 0.2,
        inherit.aes = FALSE
      ) +
      labs(
        subtitle = "Shaded areas indicate detected heatwave intervals"
      )
  }
  
  ggsave(
    filename = file.path(species_dir, paste0(file_prefix, "_ts_shading.png")),
    plot     = p_shading, width = 7, height = 5
  )
  
  # (D) (Optional) Save CSVs if desired
  # CSV prefix
  # write.csv(event_df,  file.path(species_dir, paste0(file_prefix, "_events.csv")),        row.names = FALSE)
  # write.csv(clim_df,   file.path(species_dir, paste0(file_prefix, "_climatology.csv")),   row.names = FALSE)
  # write.csv(exceed_days_2021,   file.path(species_dir, paste0(file_prefix, "_exceed_days.csv")), row.names = FALSE)
  # write.csv(heatwave_intervals_2021, file.path(species_dir, paste0(file_prefix, "_intervals.csv")), row.names = FALSE)
}

##################################################################
## 8) Combine All Individuals' Event/Clim Data
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
## 9) Save the Combined Data If Desired
##################################################################
# write.csv(all_events_df,        file.path(output_dir, "all_events.csv"), row.names = FALSE)
# write.csv(all_climatology_df,   file.path(output_dir, "all_climatology.csv"), row.names = FALSE)

##################################################################
## 10) Create a CSV for Heatwave Events in 2020 or 2021
##################################################################
# We filter all_events_df so that start/peak/end is in 2020 or 2021
hw_2020_2021 <- all_events_df %>%x
  filter(
    format(date_start, "%Y") %in% c("2020","2021") |
      format(date_peak,  "%Y") %in% c("2020","2021") |
      format(date_end,   "%Y") %in% c("2020","2021")
  )

# Write that CSV
write.csv(hw_2020_2021, file.path(output_dir, "heatwave_events_2020_2021.csv"), row.names = FALSE)

##################################################################
## Done!
## - Each species has a subfolder with PNG plots (including a map).
## - 'heatwave_events_2020_2021.csv' in the main output folder shows 
##   all events that touched 2020 or 2021 for every indivID.
##################################################################
