# 01_tables.R  — Descriptive tables (national + regional + climatology)
suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tidyr); library(readr)
  library(lubridate); library(purrr)
})

# Must load config first (shapes, file paths, helpers, periods, palettes)
# source("00_config.R")  # <- already sourced in your console; keep commented here

dir.create("out/tables", recursive = TRUE, showWarnings = FALSE)

# 1) Load stacks 
stopifnot(file.exists(f_T), file.exists(f_RH), file.exists(f_P))

T2M <- read_stack_monthly(f_T)         # °C
RH  <- read_stack_monthly(f_RH)        # %
Pmm <- read_stack_monthly(f_P)         # mm/month

# Make sure we have dates on layers
dates_T  <- as.Date(terra::time(T2M))
dates_RH <- as.Date(terra::time(RH))
dates_P  <- as.Date(terra::time(Pmm))
stopifnot(length(dates_T) > 0, length(dates_RH) > 0, length(dates_P) > 0)

# For annual aggregation
yr_T  <- year(dates_T)
yr_RH <- year(dates_RH)
yr_P  <- year(dates_P)

# terra vectors from config
togo_sv    <- togo_state_sv
regions_sv <- togo_regions_sv

# helpers 
# Country mean time series (monthly)
country_series <- function(r, dates) {
  v <- terra::extract(r, togo_sv, fun = mean, na.rm = TRUE, exact = TRUE)[, -1]
  tibble(date = dates, value = as.numeric(v[1, ]))
}

# Regional mean time series (monthly)
region_series <- function(r, dates) {
  ex <- terra::extract(r, regions_sv, fun = mean, na.rm = TRUE, exact = TRUE)
  tibble(region = togo_regions$region) |>
    bind_cols(as_tibble(ex[, -1])) |>
    pivot_longer(-region, names_to = "layer", values_to = "value") |>
    mutate(idx = as.integer(gsub("[^0-9]", "", layer)),
           date = dates[idx]) |>
    select(region, date, value)
}

# Annual aggregate from monthly
annualize <- function(df, var, agg = c("mean","sum")) {
  agg <- match.arg(agg)
  df |>
    mutate(year = year(date)) |>
    group_by(across(any_of(c("region","year")))) |>
    summarise(value = if (agg == "mean") mean(value, na.rm = TRUE)
              else sum(value, na.rm = TRUE),
              .groups = "drop") |>
    mutate(var = var)
}

# Baseline/recent split helpers
is_base   <- function(y) dplyr::between(y, baseline[1], baseline[2])
is_recent <- function(y) dplyr::between(y, recent[1],   recent[2])

# 2) NATIONAL ANNUAL SUMMARIES
T_nat_m  <- country_series(T2M, dates_T)  %>% mutate(var="T (°C)")
RH_nat_m <- country_series(RH,  dates_RH) %>% mutate(var="RH (%)")
P_nat_m  <- country_series(Pmm, dates_P)  %>% mutate(var="P (mm/month)")

T_nat_y  <- annualize(T_nat_m,  "T (°C)",     "mean") %>% rename(ann_value = value)
RH_nat_y <- annualize(RH_nat_m, "RH (%)",     "mean") %>% rename(ann_value = value)
P_nat_y  <- annualize(P_nat_m,  "P (mm/year)","sum")  %>% rename(ann_value = value)

nat_all  <- bind_rows(T_nat_y, RH_nat_y, P_nat_y)

nat_summary <- nat_all %>%
  group_by(var) %>%
  summarise(
    years_min = min(year), years_max = max(year),
    ann_min   = min(ann_value, na.rm = TRUE),
    ann_q25   = quantile(ann_value, 0.25, na.rm = TRUE),
    ann_mean  = mean(ann_value, na.rm = TRUE),
    ann_q75   = quantile(ann_value, 0.75, na.rm = TRUE),
    ann_max   = max(ann_value, na.rm = TRUE),
    base_mean = mean(ann_value[is_base(year)], na.rm = TRUE),
    recent_mean = mean(ann_value[is_recent(year)], na.rm = TRUE),
    delta_recent_minus_base = recent_mean - base_mean,
    .groups = "drop"
  )

write_csv(nat_summary, "out/tables/national_annual_summary.csv")

# 3) REGIONAL ANNUAL SUMMARIES 
T_reg_m  <- region_series(T2M, dates_T)  %>% mutate(var="T (°C)")
RH_reg_m <- region_series(RH,  dates_RH) %>% mutate(var="RH (%)")
P_reg_m  <- region_series(Pmm, dates_P)  %>% mutate(var="P (mm/month)")

T_reg_y  <- T_reg_m  %>% annualize("T (°C)",      "mean") %>% rename(ann_value = value)
RH_reg_y <- RH_reg_m %>% annualize("RH (%)",      "mean") %>% rename(ann_value = value)
P_reg_y  <- P_reg_m  %>% annualize("P (mm/year)", "sum")  %>% rename(ann_value = value)

reg_all <- bind_rows(T_reg_y, RH_reg_y, P_reg_y)

reg_summary <- reg_all %>%
  group_by(region, var) %>%
  summarise(
    years_min = min(year), years_max = max(year),
    ann_min   = min(ann_value, na.rm = TRUE),
    ann_q25   = quantile(ann_value, 0.25, na.rm = TRUE),
    ann_mean  = mean(ann_value, na.rm = TRUE),
    ann_q75   = quantile(ann_value, 0.75, na.rm = TRUE),
    ann_max   = max(ann_value, na.rm = TRUE),
    base_mean = mean(ann_value[is_base(year)], na.rm = TRUE),
    recent_mean = mean(ann_value[is_recent(year)], na.rm = TRUE),
    delta_recent_minus_base = recent_mean - base_mean,
    .groups = "drop"
  ) %>%
  arrange(region, var)

write_csv(reg_summary, "out/tables/regions_annual_summary.csv")

# 4) BASELINE MONTHLY CLIMATOLOGY (NATIONAL)
clim_monthly <- function(df_month, label) {
  df_month %>%
    filter(is_base(year(date))) %>%
    mutate(month = month(date)) %>%
    group_by(month) %>%
    summarise(
      mean  = mean(value, na.rm = TRUE),
      q10   = quantile(value, 0.10, na.rm = TRUE),
      q90   = quantile(value, 0.90, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(var = label)
}

clim_nat <- bind_rows(
  clim_monthly(T_nat_m,  "T (°C)"),
  clim_monthly(RH_nat_m, "RH (%)"),
  clim_monthly(P_nat_m,  "P (mm/month)")
) %>%
  mutate(month_name = factor(month.abb[month], levels = month.abb)) %>%
  select(var, month, month_name, mean, q10, q90)

write_csv(clim_nat, "out/tables/national_baseline_monthly_climatology.csv")

message("Tables written:",
        "\n  - out/tables/national_annual_summary.csv",
        "\n  - out/tables/regions_annual_summary.csv",
        "\n  - out/tables/national_baseline_monthly_climatology.csv")
