# 00_config.R  
suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(stringr)
})

# SHAPES
gpkg_path <- "data/togo2.gpkg" 
lyr0 <- "ADM_ADM_0"             # country
lyr1 <- "ADM_ADM_1"             # regions

stopifnot(file.exists(gpkg_path))
stopifnot(!is.na(lyr0), !is.na(lyr1))

togo_country <- sf::st_read(gpkg_path, layer = lyr0, quiet = TRUE) |>
  sf::st_make_valid() |>
  sf::st_transform(4326)

togo_regions <- sf::st_read(gpkg_path, layer = lyr1, quiet = TRUE) |>
  sf::st_make_valid() |>
  sf::st_transform(4326)

# Pick a reasonable region name column and normalize
name_candidates <- c("region","NAME_1","ADM1_EN","ADM1_FR","NL_NAME_1","NAME","ADM1_NAME","ADMIN1")
reg_col <- intersect(name_candidates, names(togo_regions))[1]
stopifnot(!is.na(reg_col))

togo_regions <- togo_regions |>
  mutate(region = str_to_title(.data[[reg_col]])) |>
  mutate(region = case_when(
    str_detect(region, "Savan") ~ "Savanes",
    str_detect(region, "Kara")  ~ "Kara",
    str_detect(region, "Centr") ~ "Centrale",
    str_detect(region, "Plate") ~ "Plateaux",
    str_detect(region, "Marit") ~ "Maritime",
    TRUE ~ region
  ))

togo_state <- sf::st_union(togo_country) # single (multi)polygon

# terra vectors for masking/cropping if needed elsewhere
togo_country_sv <- terra::vect(togo_country)
togo_regions_sv <- terra::vect(togo_regions)
togo_state_sv   <- terra::vect(togo_state)

# DATA FILES (masked GeoTIFFs from WSL pipeline)
exp_dir <- "export"
f_T  <- file.path(exp_dir, "T2M_C_TOGO_masked.tif")       # °C
f_RH <- file.path(exp_dir, "RH_pct_TOGO_masked.tif")      # %
f_P  <- file.path(exp_dir, "P_mm_month_TOGO_masked.tif")  # mm/month

# Quick presence check (non-fatal; other scripts can stopifnot)
invisible({
  if (!all(file.exists(f_T, f_RH, f_P))) {
    message("Note: one or more GeoTIFFs not found in ./export. ",
            "Run the WSL preprocessing if needed.")
  }
})

# Helper to read a multiband GeoTIFF and attach monthly time (1980-01 to …)
read_stack_monthly <- function(path, start = as.Date("1980-01-01")) {
  r <- terra::rast(path)
  n <- terra::nlyr(r)
  # assign monthly timestamps (GeoTIFFs don't carry time)
  try({
    terra::time(r) <- seq.Date(start, by = "1 month", length.out = n)
  }, silent = TRUE)
  r
}

# PLOTTING THEME & PALETTES 
theme_clean <- ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(colour = "#eaeaea", linewidth = 0.3),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position  = "bottom",
    plot.title       = ggplot2::element_text(face = "bold")
  )

# Suggested palettes (consistent, readable)
pal_temp_seq <- c("#4575b4","#74add1","#abd9e9","#fee090","#f46d43","#d73027")   # °C
pal_temp_div <- c("#313695","#4575b4","#74add1","#f7f7f7","#f46d43","#d73027")   # deltas/trends
pal_rh_seq   <- c("#2166ac","#67a9cf","#d1e5f0","#fddbc7","#ef8a62","#b2182b")  # %
pal_rh_div   <- c("#2166ac","#67a9cf","#d1e5f0","#f7f7f7","#fddbc7","#ef8a62","#b2182b")
pal_prec_seq <- c("#f7fcf0","#ccebc5","#7bccc4","#2b8cbe","#084081")            # mm/month
pal_prec_div <- c("#7b3294","#c2a5cf","#e7d4e8","#f7f7f7","#d9f0d3","#a6dba0","#008837") # % / trends

# PERIODS
baseline <- c(1981, 2010)
recent   <- c(2011, 2025)

# Make output folders if missing (used by other scripts)
dir.create("out", showWarnings = FALSE)
dir.create("out/maps", recursive = TRUE, showWarnings = FALSE)
dir.create("out/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("out/ts", recursive = TRUE, showWarnings = FALSE)

message("00_config.R loaded. Shapes ready, file paths set, palettes & theme available.")

