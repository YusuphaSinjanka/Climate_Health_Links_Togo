# 03_cycles_baseline_vs_recent.R
suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tidyr)
  library(ggplot2); library(stringr)
})

# setup
dir.create("out/ts", recursive = TRUE, showWarnings = FALSE)
source("00_config.R")  # provides togo_country, togo_regions (WGS84)

# Input rasters exported from WSL step 6
f_T  <- file.path("data","T2M_C_RECT_alpha.tif")
f_RH <- file.path("data","RH_pct_RECT_alpha.tif")
f_P  <- file.path("data","P_mm_month_RECT_alpha.tif")
stopifnot(file.exists(f_T), file.exists(f_RH), file.exists(f_P))

T2M <- rast(f_T); RH <- rast(f_RH); Pmm <- rast(f_P)

# Time axis
nly <- nlyr(T2M)
dates <- seq(as.Date("1980-01-01"), by = "1 month", length.out = nly)
stopifnot(nlyr(RH) == nly, nlyr(Pmm) == nly)

baseline <- c(1981, 2010)
recent   <- c(2011, 2025)
recent_cap <- c(2011, 2025) # used in plot labels

# helperS
layer_dates <- function(r) dates
idx_period  <- function(d, yrs) which(format(d, "%Y") %in% as.character(seq(yrs[1], yrs[2])))

# national monthly mean series from raster
nat_monthly <- function(r) {
  m <- terra::extract(r, terra::vect(togo_country), fun = mean, na.rm = TRUE, exact = TRUE)[, -1]
  tibble(date = layer_dates(r), mean = as.numeric(m[1, ])) |>
    mutate(month = as.integer(format(date, "%m")), year = as.integer(format(date, "%Y")))
}

# regional monthly mean series from raster
reg_monthly <- function(r) {
  ex <- terra::extract(r, terra::vect(togo_regions), fun = mean, na.rm = TRUE, exact = TRUE)
  d  <- as_tibble(ex)
  # glue region names
  regname_col <- intersect(c("region","NAME_1","ADM1_EN","ADM1_FR","ADM1_NAME","NAME"), names(togo_regions))[1]
  regnames <- togo_regions[[regname_col]]
  d$region <- regnames[d$ID]
  d |> select(-ID) |>
    pivot_longer(-region, names_to = "lyr", values_to = "mean") |>
    mutate(idx = as.integer(gsub("[^0-9]", "", lyr)),
           date = dates[idx],
           month = as.integer(format(date, "%m")),
           year  = as.integer(format(date, "%Y"))) |>
    select(region, date, year, month, mean)
}

# monthly climatology for a period (national)
clim_nat <- function(df, yrs) {
  i <- idx_period(df$date, yrs)
  df[i, ] |>
    group_by(month) |>
    summarise(mean = mean(mean, na.rm = TRUE), .groups = "drop")
}

# monthly climatology for a period (regional)
clim_reg <- function(df, yrs) {
  i <- idx_period(df$date, yrs)
  df[i, ] |>
    group_by(region, month) |>
    summarise(mean = mean(mean, na.rm = TRUE), .groups = "drop")
}

# compute monthly series 
T_nat_m  <- nat_monthly(T2M)
RH_nat_m <- nat_monthly(RH)
P_nat_m  <- nat_monthly(Pmm)

T_reg_m  <- reg_monthly(T2M)
RH_reg_m <- reg_monthly(RH)
P_reg_m  <- reg_monthly(Pmm)

# build climatologies (baseline vs recent)
T_nat_base  <- clim_nat(T_nat_m,  baseline)
T_nat_rec   <- clim_nat(T_nat_m,  recent)
RH_nat_base <- clim_nat(RH_nat_m, baseline)
RH_nat_rec  <- clim_nat(RH_nat_m, recent)
P_nat_base  <- clim_nat(P_nat_m,  baseline)
P_nat_rec   <- clim_nat(P_nat_m,  recent)

T_reg_base  <- clim_reg(T_reg_m,  baseline)
T_reg_rec   <- clim_reg(T_reg_m,  recent)
RH_reg_base <- clim_reg(RH_reg_m, baseline)
RH_reg_rec  <- clim_reg(RH_reg_m, recent)
P_reg_base  <- clim_reg(P_reg_m,  baseline)
P_reg_rec   <- clim_reg(P_reg_m,  recent)

# plotting helpers (robust factor handling)
plot_nat_cycle <- function(base_df, rec_df, title, ylab, outpng) {
  recent_label <- sprintf("%d–%d", recent_cap[1], recent_cap[2])
  
  df <- bind_rows(
    mutate(base_df, period = "base"),
    mutate(rec_df,  period = "recent")
  ) |>
    mutate(month_lab = factor(month.abb[month], levels = month.abb),
           period = factor(period, levels = c("base","recent"),
                           labels = c("1981–2010", recent_label)))
  
  p <- ggplot(df, aes(month_lab, mean, group = period, linetype = period)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.6) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(title = title, x = NULL, y = ylab, linetype = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major = element_line(colour = "#eaeaea", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  ggsave(outpng, p, width = 8, height = 4.8, dpi = 300)
  message("Saved: ", outpng)
}

plot_reg_cycles <- function(df_base, df_rec, title, ylab, outpng) {
  recent_label <- sprintf("%d–%d", recent_cap[1], recent_cap[2])
  
  df <- bind_rows(
    mutate(df_base, period = "base"),
    mutate(df_rec,  period = "recent")
  ) |>
    mutate(month_lab = factor(month.abb[month], levels = month.abb),
           period = factor(period, levels = c("base","recent"),
                           labels = c("1981–2010", recent_label)))
  
  p <- ggplot(df, aes(month_lab, mean,
                      colour = region,
                      linetype = period,
                      group = interaction(region, period))) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.5) +
    scale_linetype_manual(values = c("solid", "dashed")) +
    labs(title = title, x = NULL, y = ylab, colour = "Region", linetype = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major = element_line(colour = "#eaeaea", linewidth = 0.3),
          panel.grid.minor = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(face = "bold"))
  ggsave(outpng, p, width = 10, height = 6, dpi = 300)
  message("Saved: ", outpng)
}

# make plots
plot_nat_cycle(T_nat_base, T_nat_rec,
               "Togo · 2 m Temperature — annual cycle (baseline vs recent)",
               "°C", "out/ts/NAT_cycle_T_baseline_vs_recent.png")

plot_nat_cycle(RH_nat_base, RH_nat_rec,
               "Togo · Relative Humidity — annual cycle (baseline vs recent)",
               "%",  "out/ts/NAT_cycle_RH_baseline_vs_recent.png")

plot_nat_cycle(P_nat_base, P_nat_rec,
               "Togo · Precipitation — annual cycle (baseline vs recent)",
               "mm/month", "out/ts/NAT_cycle_Pmm_baseline_vs_recent.png")

plot_reg_cycles(T_reg_base, T_reg_rec,
                "Regions · 2 m Temperature — annual cycle (baseline vs recent)",
                "°C", "out/ts/REG_cycle_T_baseline_vs_recent.png")

plot_reg_cycles(RH_reg_base, RH_reg_rec,
                "Regions · Relative Humidity — annual cycle (baseline vs recent)",
                "%",  "out/ts/REG_cycle_RH_baseline_vs_recent.png")

plot_reg_cycles(P_reg_base, P_reg_rec,
                "Regions · Precipitation — annual cycle (baseline vs recent)",
                "mm/month", "out/ts/REG_cycle_Pmm_baseline_vs_recent.png")

message("All cycle plots saved in out/ts")
