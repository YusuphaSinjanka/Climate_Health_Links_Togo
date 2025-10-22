# 02_maps.R — rectangle maps with neat Togo overlay; explicit dates for GeoTIFFs

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr)
  library(ggplot2); library(scales); library(RColorBrewer)
})

dir.create("out/maps", recursive = TRUE, showWarnings = FALSE)

# ---- expects from 00_config.R ----
# T2M (°C), RH (%), Pmm (mm/month), togo_country, togo_regions
stopifnot(exists("T2M"), exists("RH"), exists("Pmm"),
          exists("togo_country"), exists("togo_regions"))

# Periods
baseline <- c(1981, 2010)
recent   <- c(2011, 2025)

# DATE HANDLING (for multiband GeoTIFFs without time) 
# We know the ERA5-Land stack is monthly from 1980-01 to 2025-09 → 549 bands
make_band_dates <- function(r, start = as.Date("1980-01-01")) {
  seq(from = start, by = "1 month", length.out = nlyr(r))
}
r_dates <- function(r) make_band_dates(r)
in_span <- function(d, span) dplyr::between(as.integer(format(d, "%Y")), span[1], span[2])

# helpers
annual_mean <- function(r, span) {
  d <- r_dates(r)
  i <- which(in_span(d, span))
  stopifnot(length(i) > 0)
  r_sel <- r[[i]]
  yfac  <- factor(format(d[i], "%Y"))
  terra::tapp(r_sel, yfac, mean, na.rm = TRUE) |> terra::mean(na.rm = TRUE)
}

# per-decade slope using simple time index; robust to monthly/annual
slope_decade <- function(r) {
  d <- r_dates(r)
  stopifnot(length(d) > 1)
  # time in (fractional) years, slope converted to per-decade
  x <- as.numeric(as.integer(format(d, "%Y")) + (as.integer(format(d, "%m")) - 0.5)/12)
  fun <- function(v) {
    if (all(is.na(v))) return(NA_real_)
    b <- tryCatch(stats::coef(stats::lm(v ~ x))[2] * 10, error = function(e) NA_real_)
    b
  }
  terra::app(r, fun)
}

# nice gg theme
theme_map <- theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_line(colour = "#e3e3e3", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.key.width = unit(2.5, "lines"))

# convert raster to df & plot over full rectangle (pad by half-pixel so edges aren’t clipped)
plot_r <- function(r, title, label, palette, limits, diverging = FALSE, fname) {
  df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
  names(df)[3] <- "val"
  
  ex  <- terra::ext(r)
  res <- terra::res(r)
  xlim <- c(ex$xmin - res[1]/2, ex$xmax + res[1]/2)
  ylim <- c(ex$ymin - res[2]/2, ex$ymax + res[2]/2)
  
  p <- ggplot(df, aes(x, y, fill = val)) +
    geom_raster(interpolate = FALSE) +
    geom_sf(data = togo_country, inherit.aes = FALSE,
            fill = NA, colour = "black", linewidth = 0.5) +
    geom_sf(data = togo_regions, inherit.aes = FALSE,
            fill = NA, colour = "black", linewidth = 0.25) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = title, x = NULL, y = NULL, fill = label) +
    theme_map
  
  if (diverging) {
    p <- p + scale_fill_gradient2(
      low = palette[1], mid = palette[2], high = palette[3],
      midpoint = 0, limits = limits, oob = scales::squish,
      labels = label_number(accuracy = NULL),
      guide = guide_colorbar(ticks.colour = "black"))
  } else {
    p <- p + scale_fill_gradientn(
      colours = palette, limits = limits, oob = scales::squish,
      labels = label_number(accuracy = NULL),
      guide = guide_colorbar(ticks.colour = "black"))
  }
  
  ggsave(fname, p, width = 6.6, height = 9.6, dpi = 300)
  message("Wrote: ", fname)
}

# palettes (your request)
pal_temp_seq <- brewer.pal(9, "YlOrRd")            # Temperature levels
pal_prec_seq <- brewer.pal(9, "YlGnBu")            # Precipitation levels
pal_rh_seq   <- brewer.pal(9, "Blues")             # RH levels (clear, readable)

# Diverging (for trends/deltas) centered on 0
pal_temp_div <- c("#4575b4", "#f7f7f7", "#d73027") # cool–white–warm
pal_prec_div <- c("#2166ac", "#f7f7f7", "#1a9850") # dry–white–wet
pal_rh_div   <- c("#2166ac", "#f7f7f7", "#b2182b") # down–white–up

# build layers 
# Baseline & recent means
t2m_base   <- annual_mean(T2M, baseline)    # °C
rh_base    <- annual_mean(RH,  baseline)    # %
p_base     <- annual_mean(Pmm, baseline)    # mm/month (avg monthly across years)

t2m_recent <- annual_mean(T2M, recent)
rh_recent  <- annual_mean(RH,  recent)
p_recent   <- annual_mean(Pmm, recent)

# Trends per decade (1980–2025)
t2m_trend <- slope_decade(T2M)              # °C/decade
rh_trend  <- slope_decade(RH)               # pp/decade

# Precip trend on annual totals (reduces seasonality)
d_all <- r_dates(Pmm)
yrfac <- factor(format(d_all, "%Y"))
P_ann <- terra::tapp(Pmm, yrfac, sum, na.rm = TRUE)      # mm/year
# Rebuild a date vector for annual layers (Jan 1 of each year present)
P_ann_dates <- as.Date(paste0(levels(yrfac), "-01-01"))
# Reuse slope_decade logic on annual data:
slope_decade_annual <- function(r_annual, dates) {
  x <- as.numeric(as.integer(format(dates, "%Y")))
  fun <- function(v) {
    if (all(is.na(v))) return(NA_real_)
    b <- tryCatch(stats::coef(stats::lm(v ~ x))[2] * 10, error = function(e) NA_real_)
    b
  }
  terra::app(r_annual, fun)
}
p_trend <- slope_decade_annual(P_ann, P_ann_dates)       # mm/decade

# Deltas (recent - baseline)
t2m_delta   <- t2m_recent - t2m_base                     # °C
rh_delta    <- rh_recent - rh_base                       # pp
p_delta_pct <- 100 * (p_recent - p_base) / p_base        # %

# pretty limits (data-driven but tidy)
rng_t  <- quantile(values(t2m_base), c(0.02, 0.98), na.rm = TRUE)
rng_rh <- quantile(values(rh_base),  c(0.02, 0.98), na.rm = TRUE)
rng_p  <- quantile(values(p_base),   c(0.02, 0.98), na.rm = TRUE)
t2m_abs_rng <- range(pretty(rng_t))
rh_abs_rng  <- range(pretty(rng_rh))
p_abs_rng   <- range(pretty(rng_p))

t2m_trend_rng <- c(-0.5, 0.5)    # °C/decade
rh_trend_rng  <- c(-5, 5)        # pp/decade
p_trend_rng   <- c(-120, 120)    # mm/decade (annual totals)

t2m_delta_rng <- c(-1.5, 1.5)    # °C
rh_delta_rng  <- c(-8, 8)        # pp
p_delta_prng  <- c(-50, 50)      # %

# plot
# Means (baseline & recent)
plot_r(t2m_base,
       sprintf("Togo rectangle · 2 m Temperature (°C) · Mean %d–%d", baseline[1], baseline[2]),
       "°C", pal_temp_seq, limits = t2m_abs_rng, diverging = FALSE,
       fname = "out/maps/t2m_mean_annual_baseline.png")

plot_r(t2m_recent,
       sprintf("Togo rectangle · 2 m Temperature (°C) · Mean %d–%d", recent[1], recent[2]),
       "°C", pal_temp_seq, limits = t2m_abs_rng, diverging = FALSE,
       fname = "out/maps/t2m_mean_annual_recent.png")

plot_r(rh_base,
       sprintf("Togo rectangle · Relative Humidity (%%) · Mean %d–%d", baseline[1], baseline[2]),
       "%", pal_rh_seq, limits = rh_abs_rng, diverging = FALSE,
       fname = "out/maps/rh_mean_annual_baseline.png")

plot_r(rh_recent,
       sprintf("Togo rectangle · Relative Humidity (%%) · Mean %d–%d", recent[1], recent[2]),
       "%", pal_rh_seq, limits = rh_abs_rng, diverging = FALSE,
       fname = "out/maps/rh_mean_annual_recent.png")

plot_r(p_base,
       sprintf("Togo rectangle · Precipitation (mm/month) · Mean %d–%d", baseline[1], baseline[2]),
       "mm/month", pal_prec_seq, limits = p_abs_rng, diverging = FALSE,
       fname = "out/maps/precip_mean_annual_baseline.png")

plot_r(p_recent,
       sprintf("Togo rectangle · Precipitation (mm/month) · Mean %d–%d", recent[1], recent[2]),
       "mm/month", pal_prec_seq, limits = p_abs_rng, diverging = FALSE,
       fname = "out/maps/precip_mean_annual_recent.png")

# Trends (per decade)
plot_r(t2m_trend,
       "Togo rectangle · Temperature trend (°C/decade) · 1980–2025",
       "°C/decade", pal_temp_div, limits = t2m_trend_rng, diverging = TRUE,
       fname = "out/maps/t2m_trend_decade.png")

plot_r(rh_trend,
       "Togo rectangle · RH trend (pp/decade) · 1980–2025",
       "pp/decade", pal_rh_div, limits = rh_trend_rng, diverging = TRUE,
       fname = "out/maps/rh_trend_decade.png")

plot_r(p_trend,
       "Togo rectangle · Precip trend (mm/decade, annual totals) · 1980–2025",
       "mm/decade", pal_prec_div, limits = p_trend_rng, diverging = TRUE,
       fname = "out/maps/precip_trend_decade.png")

# Deltas recent vs baseline
plot_r(t2m_delta,
       sprintf("Togo rectangle · ΔTemperature (°C) · (%d–%d) − (%d–%d)", recent[1], recent[2], baseline[1], baseline[2]),
       "°C", pal_temp_div, limits = t2m_delta_rng, diverging = TRUE,
       fname = "out/maps/t2m_delta_recent_vs_base.png")

plot_r(rh_delta,
       sprintf("Togo rectangle · ΔRH (pp) · (%d–%d) − (%d–%d)", recent[1], recent[2], baseline[1], baseline[2]),
       "pp", pal_rh_div, limits = rh_delta_rng, diverging = TRUE,
       fname = "out/maps/rh_delta_recent_vs_base.png")

plot_r(p_delta_pct,
       sprintf("Togo rectangle · ΔPrecip (%%) · (%d–%d) vs (%d–%d)", recent[1], recent[2], baseline[1], baseline[2]),
       "%", pal_prec_div, limits = p_delta_prng, diverging = TRUE,
       fname = "out/maps/precip_pct_anom_annual.png")



