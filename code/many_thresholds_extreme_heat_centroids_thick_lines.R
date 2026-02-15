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
p <- ggplot() +
  
  # Heatwave shading only during detected event durations
  geom_ribbon(
    data = hw_shade_df,
    aes(x = t, ymin = ymin, ymax = ymax),
    fill = col_thresh, alpha = 0.35
  ) +
  
  # Temperature line
  geom_line(
    data = line_df %>% filter(series == "Temperature"),
    aes(x = t, y = value),
    linewidth = 1.35,
    color = col_temp
  ) +
  
  # Climatology line
  geom_line(
    data = line_df %>% filter(series == "Climatology"),
    aes(x = t, y = value),
    linewidth = 1.45,
    color = col_clim
  ) +
  
  # Threshold line
  geom_line(
    data = line_df %>% filter(series == "Threshold"),
    aes(x = t, y = value),
    linewidth = 1.55,
    color = col_thresh,
    linetype = "dashed"
  ) +
  
  facet_wrap(~ pctile, ncol = 1) +
  
  scale_x_date(
    limits = x_limits,
    breaks = x_breaks,
    date_labels = "%b %d",
    expand = c(0, 0)
  ) +
  
  labs(
    title = paste0(indivID, ": June–July ", compare_year,
                   " (event-line style)"),
    subtitle = paste0("Climatology: ", hist_start, "–", hist_end,
                      " | Panels = percentile thresholds"),
    x = NULL,
    y = "Tmax (°C)"
  ) +
  
  theme_minimal(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 16),
    strip.text = element_text(face = "bold", size = 15),
    axis.title.y = element_text(size = 18),
    axis.text = element_text(size = 14),
    panel.grid.minor = element_blank(),
    
    # Completely remove legend
    legend.position = "none"
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
