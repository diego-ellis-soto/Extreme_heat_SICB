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

##################################################################
## 2) Input: Centroid(s)
##################################################################
valid_centroids_df <- data.frame(
  centroid_lon = -117.5,
  centroid_lat = 48.3,
  indivID      = "Washington"
) %>%
  filter(!is.na(centroid_lat), !is.na(centroid_lon))

##################################################################
## 3) Parameters & Output Folder
##################################################################
hist_start   <- 1991
hist_end     <- 2020
compare_year <- 2021

pctiles <- c(85, 90, 95)  # <- 3 thresholds

output_dir <- "outdir/"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

##################################################################
## 4) Helper Function for Daymet Download (June–July)
##################################################################
get_daymet_jj <- function(year, lat, lon, site_prefix) {
  tmp <- tryCatch({
    download_daymet(
      site     = paste0(site_prefix, "_", year),
      lat      = lat,
      lon      = lon,
      start    = year,
      end      = year,
      internal = TRUE,
      force    = TRUE
    )
  }, error = function(e) {
    warning(paste("Error downloading DAYMET data for:", site_prefix, year,
                  "at", lat, "/", lon, "->", e$message))
    return(NULL)
  })
  
  if (is.null(tmp)) return(NULL)
  
  df <- tmp$data %>%
    mutate(Date = as.Date(paste(year, yday), format = "%Y %j")) %>%
    filter(Date >= as.Date(paste0(year, "-06-01")),
           Date <= as.Date(paste0(year, "-07-31"))) %>%
    select(Date, tmax_c = tmax..deg.c.)
  
  df
}

##################################################################
## 5) Background Map Data
##################################################################
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

##################################################################
## 6) Main Loop Over Each Valid Centroid
##################################################################
all_events_list      <- list()
all_climatology_list <- list()

for(i in seq_len(nrow(valid_centroids_df))) {
  
  this_id <- valid_centroids_df$indivID[i]
  lat_c   <- valid_centroids_df$centroid_lat[i]
  lon_c   <- valid_centroids_df$centroid_lon[i]
  
  message(sprintf(
    "Assessing extreme heat in %s: lat=%.5f, lon=%.5f (%d/%d)",
    this_id, lat_c, lon_c, i, nrow(valid_centroids_df)
  ))
  
  ################################################################
  ## (A) Download historical data (1991–2020)
  ################################################################
  df_hist_list <- list()
  for (yr in hist_start:hist_end) {
    df_temp <- get_daymet_jj(
      year        = yr,
      lat         = lat_c,
      lon         = lon_c,
      site_prefix = paste0("hist_", this_id)
    )
    if(!is.null(df_temp)) df_hist_list[[as.character(yr)]] <- df_temp
  }
  
  if(length(df_hist_list) == 0) {
    warning(paste("No historical data for indivID=", this_id, "; skipping."))
    next
  }
  
  df_hist <- bind_rows(df_hist_list, .id = "year") %>%
    mutate(year = as.integer(year)) %>%
    arrange(Date)
  
  ################################################################
  ## (A2) Download compare year (2021)
  ################################################################
  df_2021 <- get_daymet_jj(
    year        = compare_year,
    lat         = lat_c,
    lon         = lon_c,
    site_prefix = paste0("compare_", this_id)
  )
  
  if(is.null(df_2021)) {
    warning(paste("No 2021 data for indivID=", this_id, "; skipping."))
    next
  }
  
  df_2021 <- df_2021 %>% mutate(year = compare_year)
  
  ################################################################
  ## Combine data for heatwaveR
  ################################################################
  df_combined <- bind_rows(df_hist, df_2021) %>%
    select(Date, year, tmax_c) %>%
    rename(t = Date, temp = tmax_c) %>%
    arrange(t)
  
  ################################################################
  ## (B) Run heatwaveR for THREE thresholds (85/90/95)
  ################################################################
  clims_3 <- lapply(pctiles, function(p) {
    clm <- ts2clm(
      data              = df_combined,
      x                 = t,
      y                 = temp,
      climatologyPeriod = c("1991-06-01", "2020-07-31"),
      pctile            = p
    )
    clm$pctile <- p
    clm
  })
  names(clims_3) <- paste0("p", pctiles)
  
  events_3 <- lapply(clims_3, detect_event, x = t, y = temp)
  
  ################################################################
  ## Collect combined event + climatology outputs across pctiles
  ################################################################
  clim_plot_df <- bind_rows(lapply(events_3, function(hw) hw$climatology)) %>%
    mutate(
      pctile = factor(pctile, levels = pctiles),
      year   = as.integer(format(t, "%Y")),
      indivID = this_id
    )
  
  event_all_df <- bind_rows(lapply(events_3, function(hw) hw$event)) %>%
    mutate(
      pctile = factor(pctile, levels = pctiles),
      indivID = this_id
    )
  
  all_events_list[[i]]      <- event_all_df
  all_climatology_list[[i]] <- clim_plot_df
  
  file_prefix <- this_id
  
  ################################################################
  ## (C) Plot thresholds side-by-side (facets)
  ################################################################
  p_thresh_facets <- ggplot(clim_plot_df, aes(x = t)) +
    geom_ribbon(aes(ymin = seas, ymax = thresh), alpha = 0.2) +
    geom_line(aes(y = temp), linewidth = 0.5) +
    geom_line(aes(y = thresh), linetype = "dashed") +
    facet_wrap(~ pctile, nrow = 1) +
    labs(
      title = paste0(this_id, ": June–July temps with percentile thresholds"),
      subtitle = paste0("Climatology: ", hist_start, "–", hist_end,
                        " | Compare year: ", compare_year),
      x = NULL, y = "Tmax (°C)"
    ) +
    theme_minimal()
  
  ggsave(
    filename = file.path(output_dir, paste0(file_prefix, "_thresholds_85_90_95_facets.png")),
    plot     = p_thresh_facets,
    width    = 12, height = 4
  )
  
  ################################################################
  ## (C2) Plot thresholds overlaid (single plot)
  ################################################################
  p_thresh_overlay <- ggplot(clim_plot_df, aes(x = t)) +
    geom_line(aes(y = temp), linewidth = 0.5) +
    geom_line(aes(y = thresh, color = pctile), linetype = "dashed", linewidth = 0.7) +
    labs(
      title = paste0(this_id, ": 85/90/95 thresholds overlaid"),
      subtitle = paste0("Climatology: ", hist_start, "–", hist_end,
                        " | Compare year: ", compare_year),
      x = NULL, y = "Tmax (°C)", color = "Percentile"
    ) +
    theme_minimal()
  
  ggsave(
    filename = file.path(output_dir, paste0(file_prefix, "_thresholds_85_90_95_overlay.png")),
    plot     = p_thresh_overlay,
    width    = 8, height = 5
  )
  
  ################################################################
  ## (C3) OPTIONAL: event_line() for each threshold
  ################################################################
  for(p in pctiles) {
    hw_obj <- events_3[[paste0("p", p)]]
    has_events <- nrow(hw_obj$event) > 0
    
    if(has_events) {
      p_event_line <- tryCatch({
        event_line(hw_obj, spread = 10) +
          ggtitle(paste0(this_id, " event_line (pctile=", p, ")"))
      }, error = function(e) {
        warning(paste("event_line failed for", this_id, "pctile", p, "->", e$message))
        return(NULL)
      })
      
      if(!is.null(p_event_line)) {
        ggsave(
          filename = file.path(output_dir, paste0(file_prefix, "_event_line_p", p, ".png")),
          plot     = p_event_line, width = 7, height = 5
        )
      }
    } else {
      message("No events found for indivID=", this_id, " at pctile=", p, " -> skipping event_line")
    }
  }
  
  ################################################################
  ## (C4) Consecutive exceedance intervals in compare_year for each threshold
  ################################################################
  exceed_intervals_all <- lapply(pctiles, function(p) {
    dfp <- clim_plot_df %>%
      filter(pctile == p,
             format(t, "%Y") == as.character(compare_year),
             temp > thresh) %>%
      arrange(t) %>%
      mutate(
        day_diff  = as.integer(difftime(t, lag(t), units = "days")),
        new_group = ifelse(is.na(day_diff) | day_diff != 1, 1, 0),
        group     = cumsum(new_group)
      )
    
    if(nrow(dfp) == 0) {
      return(tibble(
        pctile = p, start_date = as.Date(NA), end_date = as.Date(NA),
        duration_days = NA_integer_, indivID = this_id
      ))
    }
    
    dfp %>%
      group_by(group) %>%
      summarize(
        start_date    = first(t),
        end_date      = last(t),
        duration_days = as.integer(difftime(last(t), first(t), units = "days")) + 1,
        .groups       = "drop"
      ) %>%
      mutate(pctile = p, indivID = this_id)
  })
  
  heatwave_intervals_compare_year <- bind_rows(exceed_intervals_all)
  
  write.csv(
    heatwave_intervals_compare_year,
    file.path(output_dir, paste0(file_prefix, "_intervals_", compare_year, "_p85_90_95.csv")),
    row.names = FALSE
  )
  
  ################################################################
  ## (D) Map for centroid (fixes your earlier lon/lat column bug)
  ################################################################
  usa_states <- rnaturalearth::ne_states(
    country = "United States of America",
    returnclass = "sf"
  )
  
  wa <- usa_states %>%
    filter(name == "Washington") %>%
    st_transform(4326)
  
  # NOTE: you only have ONE point in valid_centroids_df right now,
  # so "points" and "centroid" are the same. This is just to show it correctly.
  cent_sf <- st_as_sf(
    data.frame(indivID = this_id, centroid_lon = lon_c, centroid_lat = lat_c),
    coords = c("centroid_lon", "centroid_lat"),
    crs = 4326
  )
  
  bb <- st_bbox(wa)
  pad_x <- (bb$xmax - bb$xmin) * 0.05
  pad_y <- (bb$ymax - bb$ymin) * 0.05
  xlim <- c(bb$xmin - pad_x, bb$xmax + pad_x)
  ylim <- c(bb$ymin - pad_y, bb$ymax + pad_y)
  
  p_map_wa <- ggplot() +
    geom_sf(data = wa, fill = NA, color = "white", linewidth = 0.8) +
    geom_sf(data = cent_sf, size = 4.2, shape = 21,
            fill = "red", color = "white", stroke = 0.8) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(
      title = "June–July 2021 Locations & Centroid",
      subtitle = "Washington State outline + centroid",
      x = NULL, y = NULL
    ) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(color = "grey20")
    )
  
  ggsave(
    filename = file.path(output_dir, paste0(file_prefix, "_WA_centroid_map.png")),
    plot     = p_map_wa,
    width    = 7, height = 5
  )
  
  ################################################################
  ## (E) Optional CSV outputs (uncomment if you want)
  ################################################################
  # write.csv(event_all_df, file.path(output_dir, paste0(file_prefix, "_events_p85_90_95.csv")), row.names = FALSE)
  # write.csv(clim_plot_df, file.path(output_dir, paste0(file_prefix, "_climatology_p85_90_95.csv")), row.names = FALSE)
  
}

##################################################################
## 7) Combine All Individuals' Event/Clim Data
##################################################################
all_events_df      <- bind_rows(all_events_list)
all_climatology_df <- bind_rows(all_climatology_list)

head(all_events_df)
head(all_climatology_df)

##################################################################
## 8) Create a CSV of Heatwave Periods in 2020 & 2021
##################################################################
hw_2020_2021 <- all_events_df %>%
  filter(
    year(date_start) %in% c(2020, 2021) |
      year(date_peak) %in% c(2020, 2021) |
      year(date_end)  %in% c(2020, 2021)
  )

write.csv(
  hw_2020_2021,
  file.path(output_dir, "extreme_heat_example_case_figure_1_v2.csv"),
  row.names = FALSE
)

