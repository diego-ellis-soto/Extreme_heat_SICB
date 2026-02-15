##################################################################
## Extreme heat (Daymet + heatwaveR) — stacked ggplot event-line figure
## Includes:
## ✅ Heatwave-duration shading (only during detected events)
## ✅ Color-accessible palette (Okabe–Ito-ish; red/orange kept per request)
## ✅ Legend inside the top panel (automatic: only shown in top facet)
## ✅ Thicker lines + larger text
## ✅ X-axis restricted to real temperature dates (no blank months)
##################################################################

library(daymetr)
library(dplyr)
library(lubridate)
library(heatwaveR)
library(ggplot2)
library(tidyr)
require(sf)
library(grid)
library(gtable)
##################################################################
## Parameters
##################################################################
centroid_lon <- -117.5
centroid_lat <- 48.3
indivID <- "Washington"

hist_start   <- 1991
hist_end     <- 2020
compare_year <- 2021
pctiles <- c(85, 90, 95)

output_dir <- "outdir"
dir.create(output_dir, showWarnings = FALSE)

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

##################################################################
## Helper function — download June–July Daymet
##################################################################
get_daymet_jj <- function(year, lat, lon, site_prefix) {
  
  tmp <- download_daymet(
    site     = paste0(site_prefix, "_", year),
    lat      = lat,
    lon      = lon,
    start    = year,
    end      = year,
    internal = TRUE,
    force    = TRUE
  )
  
  tmp$data %>%
    mutate(Date = as.Date(paste(year, yday), "%Y %j")) %>%
    filter(Date >= as.Date(paste0(year, "-06-01")),
           Date <= as.Date(paste0(year, "-07-31"))) %>%
    select(Date, tmax_c = tmax..deg.c.)
}

##################################################################
## Download historical data
##################################################################
hist_list <- lapply(hist_start:hist_end, function(yr) {
  get_daymet_jj(yr, centroid_lat, centroid_lon, "hist")
})

df_hist <- bind_rows(hist_list, .id = "year") %>%
  mutate(year = as.integer(year))

##################################################################
## Download comparison year
##################################################################
df_comp <- get_daymet_jj(compare_year, centroid_lat, centroid_lon, "comp") %>%
  mutate(year = compare_year)

##################################################################
## Combine for heatwaveR
##################################################################
df_combined <- bind_rows(df_hist, df_comp) %>%
  rename(t = Date, temp = tmax_c) %>%
  arrange(t)

##################################################################
## Compute climatologies + thresholds + events
##################################################################
clims <- lapply(pctiles, function(p) {
  ts2clm(
    data              = df_combined,
    x                 = t,
    y                 = temp,
    climatologyPeriod = c("1991-06-01", "2020-07-31"),
    pctile            = p
  )
})
names(clims) <- paste0("p", pctiles)

events <- lapply(clims, detect_event, x = t, y = temp)

##################################################################
## Build plotting dataframe (climatology) + event intervals
##################################################################
clim_plot_df <- bind_rows(lapply(names(events), function(nm) {
  p <- as.integer(sub("^p","", nm))
  events[[nm]]$climatology %>% mutate(pctile = p)
}))

event_df <- bind_rows(lapply(names(events), function(nm) {
  p <- as.integer(sub("^p","", nm))
  ed <- events[[nm]]$event
  if (nrow(ed) == 0) return(NULL)
  ed %>%
    transmute(
      pctile = p,
      date_start = as.Date(date_start),
      date_end   = as.Date(date_end),
      duration   = duration
    )
}))

##################################################################
## Only comparison year for plotting
##################################################################
plot_df <- clim_plot_df %>%
  mutate(year = year(t)) %>%
  filter(year == compare_year) %>%
  mutate(pctile = factor(pctile, levels = pctiles))

##################################################################
## ✅ X-axis: only where real temperature exists (no blank months)
##################################################################
x_limits <- range(plot_df$t[!is.na(plot_df$temp)], na.rm = TRUE)

# weekly ticks for ~2 months
x_breaks <- seq(
  floor_date(x_limits[1], "week"),
  ceiling_date(x_limits[2], "week"),
  by = "1 week"
)

##################################################################
## ✅ Heatwave duration shading properly:
## Shade only DURING detected heatwave intervals, between thresh and temp.
##################################################################
# Create an interval table (one row per day inside each event interval)
# then join to plot_df to get temp/thresh on those dates.
if (!is.null(event_df) && nrow(event_df) > 0) {
  
  event_days <- bind_rows(lapply(seq_len(nrow(event_df)), function(i) {
    tibble(
      pctile = event_df$pctile[i],
      t = seq(event_df$date_start[i], event_df$date_end[i], by = "day")
    )
  })) %>%
    mutate(pctile = factor(pctile, levels = pctiles))
  
  hw_shade_df <- plot_df %>%
    inner_join(event_days, by = c("pctile", "t")) %>%
    filter(!is.na(temp), !is.na(thresh), temp > thresh) %>%
    transmute(t, pctile, ymin = thresh, ymax = temp)
  
} else {
  hw_shade_df <- plot_df[0, c("t","pctile","thresh","temp")] %>%
    transmute(t, pctile, ymin = thresh, ymax = temp)
}

##################################################################
## Prepare line data
##################################################################
line_df <- plot_df %>%
  select(t, pctile, temp, seas, thresh) %>%
  pivot_longer(
    cols = c(temp, seas, thresh),
    names_to = "series",
    values_to = "value"
  ) %>%
  mutate(series = recode(series,
                         temp   = "Temperature",
                         seas   = "Climatology",
                         thresh = "Threshold"))

##################################################################
## ✅ Legend inside panel (automatic):
## We add points ONLY in the top facet to carry the legend,
## and suppress legends from the main line layers.
##################################################################
top_pctile <- levels(plot_df$pctile)[1]  # "85" given your levels
legend_key_df <- tibble(
  t = x_limits[1],
  value = 0,
  pctile = factor(top_pctile, levels = pctiles),
  series = factor(c("Temperature","Climatology","Threshold"),
                  levels = c("Temperature","Climatology","Threshold"))
)

##################################################################
## ✅ Colors (accessible-ish, while keeping red threshold & orange climatology)
##################################################################
col_temp   <- "#222222"   # near-black
col_clim   <- "#E69F00"   # Okabe–Ito orange
col_thresh <- "#D55E00"   # Okabe–Ito vermillion (still reads as red)

##################################################################
## Plot (stacked panels)
##################################################################
bottom_pctile <- max(pctiles)  # 95 in your case

p <- ggplot() +
  
  # Heatwave shading
  geom_ribbon(
    data = hw_shade_df,
    aes(x = t, ymin = ymin, ymax = ymax),
    fill = col_thresh, alpha = 0.35
  ) +
  
  # ----- Temperature line -----
geom_line(
  data = line_df %>% filter(series == "Temperature",
                            pctile != bottom_pctile),
  aes(x = t, y = value),
  linewidth = 1.35,
  color = col_temp
) +
  
  geom_line(
    data = line_df %>% filter(series == "Temperature",
                              pctile == bottom_pctile),
    aes(x = t, y = value, color = "Temperature"),
    linewidth = 1.35
  ) +
  
  # ----- Climatology -----
geom_line(
  data = line_df %>% filter(series == "Climatology",
                            pctile != bottom_pctile),
  aes(x = t, y = value),
  linewidth = 1.45,
  color = col_clim
) +
  
  geom_line(
    data = line_df %>% filter(series == "Climatology",
                              pctile == bottom_pctile),
    aes(x = t, y = value, color = "Climatology"),
    linewidth = 1.45
  ) +
  
  # ----- Threshold -----
geom_line(
  data = line_df %>% filter(series == "Threshold",
                            pctile != bottom_pctile),
  aes(x = t, y = value),
  linewidth = 1.55,
  color = col_thresh,
  linetype = "dashed"
) +
  
  geom_line(
    data = line_df %>% filter(series == "Threshold",
                              pctile == bottom_pctile),
    aes(x = t, y = value, color = "Threshold"),
    linewidth = 1.55,
    linetype = "dashed"
  ) +
  
  facet_wrap(~ pctile, ncol = 1) +
  
  scale_color_manual(
    values = c(
      Temperature = col_temp,
      Climatology = col_clim,
      Threshold   = col_thresh
    )
  ) +
  
  scale_x_date(
    limits = x_limits,
    breaks = x_breaks,
    date_labels = "%b %d",
    expand = c(0, 0)
  ) +
  
  labs(
    title = paste0(indivID, ": June–July ", compare_year),
    subtitle = paste0("Climatology: ", hist_start, "–", hist_end),
    x = NULL,
    y = "Tmax (°C)",
    color = NULL
  ) +
  
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 15),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    
    # Legend appears naturally in bottom facet
    legend.position = "bottom"
  )


##################################################################
## (MAP) State outline + centroid point (auto-detect state)
## - No USA map
## - No hard-coded state name
## - Centroid can be mean of multiple locations (if provided)
##################################################################

library(sf)
library(rnaturalearth)
library(dplyr)

# ---------------------------------------------------------------
# 1) Provide locations (one or many). If you only have one point,
#    just keep one row.
# ---------------------------------------------------------------
locations_df <- data.frame(
  lon = centroid_lon,
  lat = centroid_lat
)

# If you *do* have multiple locations, use them instead, e.g.:
# locations_df <- data.frame(lon = c(...), lat = c(...))

# ---------------------------------------------------------------
# 2) Compute centroid as MEAN of locations (flexible)
# ---------------------------------------------------------------
centroid_lon_mean <- mean(locations_df$lon, na.rm = TRUE)
centroid_lat_mean <- mean(locations_df$lat, na.rm = TRUE)

cent_sf <- st_as_sf(
  data.frame(indivID = indivID, lon = centroid_lon_mean, lat = centroid_lat_mean),
  coords = c("lon", "lat"),
  crs = 4326
)

# Optional: also keep the original points if you want to show them
pts_sf <- st_as_sf(locations_df, coords = c("lon", "lat"), crs = 4326)

# ---------------------------------------------------------------
# 3) Load US states and find which state contains the centroid
# ---------------------------------------------------------------
usa_states <- rnaturalearth::ne_states(
  country = "United States of America",
  returnclass = "sf"
) %>%
  st_transform(4326)

# Find containing state (spatial join)
state_hit <- st_join(cent_sf, usa_states, join = st_within, left = FALSE)

if (nrow(state_hit) == 0) {
  stop("Centroid point is not within any US state polygon (check lon/lat or CRS).")
}

# Extract that state's polygon
state_name <- state_hit$name[1]
state_poly <- usa_states %>% filter(name == state_name)

# ---------------------------------------------------------------
# 4) Plot state + centroid (and optionally all points)
# ---------------------------------------------------------------
bb <- st_bbox(state_poly)
pad_x <- (bb$xmax - bb$xmin) * 0.05
pad_y <- (bb$ymax - bb$ymin) * 0.05

p_state_centroid <- ggplot() +
  
  # State polygon with black outline
  geom_sf(
    data = state_poly,
    fill = "grey90",
    color = "black",     # ← black outline
    linewidth = 1.0      # ← thicker border for visibility
  ) +
  
  # Optional: original locations
  geom_sf(
    data = pts_sf,
    color = "grey30",
    size = 2,
    alpha = 0.7
  ) +
  
  # BIG centroid point
  geom_sf(
    data = cent_sf,
    shape = 21,
    size = 14,            # ← much bigger centroid
    fill = "red3",
    color = "white",
    stroke = 2         # thicker edge around centroid
  ) +
  
  coord_sf(
    xlim = c(bb$xmin - pad_x, bb$xmax + pad_x),
    ylim = c(bb$ymin - pad_y, bb$ymax + pad_y),
    expand = FALSE
  ) +
  labs(
    title = paste0(state_name, " (state detected from centroid)"),
    subtitle = paste0(
      "Centroid (mean): ",
      round(centroid_lat_mean, 3), ", ",
      round(centroid_lon_mean, 3)
    ),
    x = NULL, y = NULL
  ) +
  theme_void(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 12)
  )

ggsave(
  filename = file.path(output_dir, paste0(indivID, "_state_centroid_", run_stamp, ".png")),
  plot = p_state_centroid,
  width = 6.5, height = 5.2, dpi = 300
)


##################################################################
## Save plot + save posthoc object (no overwriting)
##################################################################
ggsave(
  filename = file.path(output_dir,
                       paste0(indivID, "_eventline_stacked_accessible_", run_stamp, ".png")),
  plot = p,
  width = 10, height = 10, dpi = 300
)

saveRDS(
  list(plot_df = plot_df, event_df = event_df, ggplot_object = p),
  file.path(output_dir,
            paste0(indivID, "_plot_object_", run_stamp, ".rds"))
)

print(p)
# 
# # print(p_state_centroid)
# # 
# # map_grob <- ggplotGrob(p_state_centroid)
# 
# p <- p +
#   theme(
#     plot.subtitle = element_text(margin = margin(r = 120))  # increase if needed
#   )
# 
# # Convert plots to grobs
# g_main <- ggplotGrob(p)
# 
# # Clean map (no margins) so it sits nicely as an inset
# g_map <- ggplotGrob(
#   p_state_centroid +
#     theme_void() +
#     theme(plot.margin = margin(0, 0, 0, 0))
# )
# 
# # Find the "subtitle" row (falls back to "title" if no subtitle)
# use_row <- if ("subtitle" %in% g_main$layout$name) "subtitle" else "title"
# row_info <- g_main$layout[g_main$layout$name == use_row, ][1, ]
# 
# t <- row_info$t
# b <- row_info$b
# l <- row_info$l
# r <- row_info$r
# 
# # Place the map in the top-right portion of that row
# # (tweak 'w_cols' to make it wider/narrower)
# w_cols <- 4
# g_main <- gtable_add_grob(
#   g_main,
#   grobs = g_map,
#   t = t, b = b,
#   l = r - w_cols, r = r,
#   z = Inf
# )
# 
# # Draw (and save)
# outfile <- file.path(output_dir, paste0(indivID, "_eventline_with_header_map_", run_stamp, ".png"))
# png(outfile, width = 900, height = 1050, res = 150)
# grid.newpage()
# grid.draw(g_main)
# dev.off()
# 
# # Optional: show in viewer
# grid.newpage()
# grid.draw(g_main)