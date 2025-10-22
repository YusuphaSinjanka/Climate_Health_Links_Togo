# 04_extremes_monthly.R
suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tidyr)
  library(lubridate); library(readr); library(stringr); library(purrr)
})

dir.create("out/csv", recursive = TRUE, showWarnings = FALSE)

# 0) Inputs & shapes (expects you've run 00_config.R previously)
#    - If not, we set safe defaults here.

if (!exists("baseline")) baseline <- c(1981, 2010)
if (!exists("recent"))   recent   <- c(2011, 2025)

# Files produced in WSL and copied to ./data
f_T  <- file.path("data","T2M_C_RECT_alpha.tif")
f_RH <- file.path("data","RH_pct_RECT_alpha.tif")
f_P  <- file.path("data","P_mm_month_RECT_alpha.tif")

stopifnot(file.exists(f_T), file.exists(f_RH), file.exists(f_P))

# Load rasters
T2M <- rast(f_T)           # °C, monthly 1980-01 .. 2025-09 (549 bands)
RH  <- rast(f_RH)          # %, same time axis
Pmm <- rast(f_P)           # mm/month, same time axis

# Dates from band count (ERA5-Land monthly starts 1980-01)
n_bands <- nlyr(T2M)
dates   <- seq(as.Date("1980-01-01"), by = "1 month", length.out = n_bands)

# 1) Shapes — read your working GPKG (from 00_config step you fixed)

gpkg_path <- "data/togo2.gpkg"
stopifnot(file.exists(gpkg_path))

# The layer names you told me to use after the fix
lyr0 <- "ADM_ADM_0"  # country
lyr1 <- "ADM_ADM_1"  # regions

togo_country <- sf::st_read(gpkg_path, layer = lyr0, quiet = TRUE) |> sf::st_make_valid()
togo_regions <- sf::st_read(gpkg_path, layer = lyr1, quiet = TRUE) |> sf::st_make_valid()

# Normalize CRS
togo_country <- sf::st_transform(togo_country, 4326)
togo_regions <- sf::st_transform(togo_regions, 4326)

# Ensure a clean region name column
name_candidates <- c("region","NAME_1","ADM1_EN","ADM1_FR","NL_NAME_1","NAME","ADM1_NAME","ADMIN1")
has <- intersect(name_candidates, names(togo_regions))
stopifnot(length(has) > 0)
togo_regions <- togo_regions |>
  dplyr::mutate(region = stringr::str_to_title(.data[[has[1]]])) |>
  dplyr::mutate(region = dplyr::case_when(
    stringr::str_detect(region, "Savan") ~ "Savanes",
    stringr::str_detect(region, "Kara")  ~ "Kara",
    stringr::str_detect(region, "Centr") ~ "Centrale",
    stringr::str_detect(region, "Plate") ~ "Plateaux",
    stringr::str_detect(region, "Marit") ~ "Maritime",
    TRUE ~ region
  ))

# Terra vectors
togo_sv    <- terra::vect(togo_country)
regions_sv <- terra::vect(togo_regions)

# 2) Helpers
layer_dates <- function(r) {
  # Reuse the global dates since band counts match
  # (If you ever change the time range per raster, replace with terra::time)
  seq(as.Date("1980-01-01"), by = "1 month", length.out = nlyr(r))
}

idx_period <- function(d, yrs) which(format(d, "%Y") %in% as.character(seq(yrs[1], yrs[2])))

# Extract NATIONAL monthly means
nat_monthly <- function(r) {
  d <- layer_dates(r)
  m <- terra::extract(r, togo_sv, fun = mean, na.rm = TRUE, exact = TRUE)[,-1]
  tibble(date = d, value = as.numeric(m[1, ]))
}

# Extract REGIONAL monthly means
reg_monthly <- function(r) {
  d  <- layer_dates(r)
  ex <- terra::extract(r, regions_sv, fun = mean, na.rm = TRUE, exact = TRUE)
  ex <- ex %>% mutate(region = togo_regions$region) %>% relocate(region)
  ex %>%
    pivot_longer(-region, names_to = "layer", values_to = "value") %>%
    mutate(idx = as.integer(gsub("[^0-9]", "", layer)),
           date = d[idx]) %>%
    select(region, date, value)
}

# Count months per YEAR satisfying a predicate
count_months_by_year_nat <- function(df_month, predicate) {
  df_month %>%
    mutate(hit = predicate(value), year = year(date)) %>%
    group_by(year) %>%
    summarise(count = sum(hit, na.rm = TRUE), .groups = "drop") %>%
    mutate(var = unique(df_month$var))
}

count_months_by_year_reg <- function(df_month, predicate) {
  df_month %>%
    mutate(hit = predicate(value), year = year(date)) %>%
    group_by(region, year) %>%
    summarise(count = sum(hit, na.rm = TRUE), .groups = "drop") %>%
    mutate(var = unique(df_month$var))
}

# Robust period summary: averages of yearly counts for baseline vs recent
period_summary <- function(df_counts, base = baseline, rec = recent) {
  out <- df_counts %>%
    mutate(period = case_when(
      between(year, base[1], base[2]) ~ "base",
      between(year, rec[1],  rec[2])  ~ "recent",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(period)) %>%
    group_by(across(any_of(c("region","var","period")))) %>%
    summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = period, values_from = mean_count)
  
  if (!"base"   %in% names(out)) out$base   <- NA_real_
  if (!"recent" %in% names(out)) out$recent <- NA_real_
  
  out %>%
    mutate(
      `1981-2010` = base,
      `2011-2025` = recent,
      delta_recent_minus_base = recent - base
    ) %>%
    select(any_of(c("region","var")), `1981-2010`, `2011-2025`, delta_recent_minus_base)
}


# 3) Build monthly series (national & regional)

T_nat_m  <- nat_monthly(T2M) %>% mutate(var = "Hot months (T ≥ 30°C)")
RH_nat_m <- nat_monthly(RH)  %>% mutate(var = "Humid months (RH ≥ 80%)")
P_nat_m  <- nat_monthly(Pmm) %>% mutate(var = "Very wet (P ≥ 200 mm)")  # we will also do very-dry separately

T_reg_m  <- reg_monthly(T2M) %>% mutate(var = "Hot months (T ≥ 30°C)")
RH_reg_m <- reg_monthly(RH)  %>% mutate(var = "Humid months (RH ≥ 80%)")
P_reg_m  <- reg_monthly(Pmm)%>% mutate(var = "Very wet (P ≥ 200 mm)")


# 4) Define “extreme month” predicates (tune thresholds here)

is_hot     <- function(x) x >= 30            # °C
is_humid   <- function(x) x >= 80            # %
is_verywet <- function(x) x >= 200           # mm/month
is_verydry <- function(x) x <= 20            # mm/month (and ≥ 0 by construction)


# 5) Count months per year (national + regional)

T_nat_y   <- count_months_by_year_nat(T_nat_m,  is_hot)
RH_nat_y  <- count_months_by_year_nat(RH_nat_m, is_humid)
P_nat_wy  <- count_months_by_year_nat(P_nat_m,  is_verywet)

# for very-dry precip, rebuild with appropriate label
P_nat_dm  <- nat_monthly(Pmm) %>% mutate(var = "Very dry (P ≤ 20 mm)")
P_nat_dy  <- count_months_by_year_nat(P_nat_dm, is_verydry)

T_reg_y   <- count_months_by_year_reg(T_reg_m,  is_hot)
RH_reg_y  <- count_months_by_year_reg(RH_reg_m, is_humid)

P_reg_wy  <- count_months_by_year_reg(P_reg_m,  is_verywet)
P_reg_dm  <- reg_monthly(Pmm) %>% mutate(var = "Very dry (P ≤ 20 mm)")
P_reg_dy  <- count_months_by_year_reg(P_reg_dm, is_verydry)

# 6) Period summaries (baseline vs recent) — robust to missing years

T_nat_sum  <- period_summary(T_nat_y)
RH_nat_sum <- period_summary(RH_nat_y)
P_nat_wsum <- period_summary(P_nat_wy)
P_nat_dsum <- period_summary(P_nat_dy)

T_reg_sum  <- period_summary(T_reg_y)
RH_reg_sum <- period_summary(RH_reg_y)
P_reg_wsum <- period_summary(P_reg_wy)
P_reg_dsum <- period_summary(P_reg_dy)


