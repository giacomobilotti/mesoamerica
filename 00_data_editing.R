#### Load and convert data for further use ----
# the output data are provided as supplement

# load libraries
library(readxl)
library(sf)
library(dplyr)

# set up directories
sourcedir <- file.path('data', 'raw_data')
targetdir <- file.path('data', 'derived_data')

## Yautepec survey data ----

# data is stored with metadata from the source reference (Smith et al 2020) in the 
# Yautepec-Survey folder

# For the analyses, we will only need the data stored in the 2-PhaseRecode sheet 
# of the YV-DataTables-2020.xlsx

YV_Data <- read_excel(
  file.path(sourcedir, 'Yautepec-Survey', 'YV-DataTables-2020.xlsx'),
  sheet = "2-PhaseRecode"
  )

# # take a look
# View(YV_Data)

## convert to sf
# CRS is UTM system. Zone 14, WGS 1984, which is EPSG:32614 

YV_sf <- st_as_sf(YV_Data, coords = c('Este', 'Norte'), crs = st_crs(32614) )

## Assign absolute chronology (start/end)
# Chronology table after Smith et al 2020

# Helper function to convert BC/AD to BP
# It is an adaptation of the rcarbon::BCADtoBP function 

dates_to_bp <- function (x) {
  index <- !is.na(x)
  if (any(x[index] == 0)) {
    stop("0 BC/AD is not a valid year.")
  }
  res <- matrix(c(x, rep(NA, length(x))), ncol = 2)
  res[index & x > 0, 2] <- abs(res[index & x > 0, 1] - 1950)
  res[index & x < 0, 2] <- abs(res[index & x < 0, 1] - 1949) # e.g.,: 100 BC is 100 years away from 1 AD not 0 (which does not exist)
  return(res[, 2])
}

chrono_tbl <- data.frame(
  Fase = c("Col", "S", "Lp2", "Lp1", "MP", "EP", "Epi",
           "LC", "MC", "EC", "TF", "LF", "MF", "EF"),
  start = dates_to_bp(c(1650, 1520, 1440, 1300, 1150, 850, 600,
               450, 300, 200, -100, -500, -1100, -1500)),
  end   = dates_to_bp(c(1820, 1650, 1520, 1440, 1300, 1150, 850,
               600, 450, 300, 200, -100, -500, -1100))
)

# join tables
YV_sf <- YV_sf |>
  left_join(chrono_tbl[, c("Fase", "start", "end")], by = "Fase")

# save 
st_write(
  obj = YV_sf, 
  file = file.path(targetdir, 'sites.gpkg'),
  layer = 'yautepec',
  append = FALSE,
)

## Postclassic cities
# this dataset includes several cities from Mexico and Central America. See 
# Raw data is stored in the PostclassicCitySize folder alongside the original metadata

# there are different files and tables

## load spatial data 
pc_sites <- read.csv(
  file = file.path(sourcedir, 'PostclassicCitySize', 'PostclassicCitySize-RC_Georef.csv')
) 

## Load epicentre data
# Samp 1
epi_sites <- read_excel(
  file.path(sourcedir, 'PostclassicCitySize', 'PostclassicCitySize-Data.xlsx'),
  sheet = "Samp1-SizeOnly", skip = 1)
# remove empty row and subset cols
epi_sites <- epi_sites[-1,c(1:3, 5, 7:10)]
colnames(epi_sites)[3:4] <- c('ha', 'ha_epi')

## Additional data
# Samp 3

samp3 <- read_excel(
  file.path(sourcedir, 'PostclassicCitySize', 'PostclassicCitySize-Data.xlsx'),
  sheet = "Samp3") |>
  filter(!is.na(area)) 
samp3 <- samp3[,-8]

## combine non spatial data
# The samp3 has 7 more entries, add the epicentre area to 
samp_epi <- samp3 |>
  left_join(
    epi_sites %>% select(Code, ha_epi),
    by = "Code"
  )

# ## To actually check if the data is the same you can run this
# compare_df <- full_join(
#   samp3,
#   epi_sites,
#   by = "Code",
#   suffix = c(".samp", ".epi")
# )
# # check values  
# compare_df <- compare_df |>
#   mutate(
#     area_match = dplyr::near(area, ha, tol = 1e-6),
#     pop1_match = dplyr::near(pop1.samp, pop1.epi, tol = 1e-6),
#     pop2_match = dplyr::near(pop2.samp, pop2.epi, tol = 1e-6),
#     pop3_match = dplyr::near(pop3.samp, pop3.epi, tol = 1e-6)
#     )
# 
# mismatches <- compare_df |>
#   filter(!(area_match & pop1_match & pop2_match & pop3_match & pop4_match)) |>
#   select(Code, starts_with("area"), starts_with("pop"))
# # In two cases the entries do not match: Calixtlahuaca (c02) and Mayapan (y10)
# # Codes only in samp3
# only_samp3 <- setdiff(samp3$Code, epi_sites$Code)
# # 7 sites
# # Codes only in epi_sites
# only_epi <- setdiff(epi_sites$Code, samp3$Code)
# # no entries

## Join spatial data
samp_epi_sf <- samp_epi|>
  right_join(
    pc_sites |> select(id, site, region, epicenter_ha, site_area_ha, population, longitude, latitude),
    by = c("Code" = "id")
  )
# fill NA in epicentre data if in either of the two columns the information is recorded
samp_epi_sf <- samp_epi_sf |>
  mutate(
    Site = coalesce(
      as.character(Site),
      as.character(site)
    ),
    ha_epi = coalesce(
      as.numeric(ha_epi),
      as.numeric(epicenter_ha)
    ),
    area = coalesce(
      as.numeric(area),
      as.numeric(site_area_ha)
    ),
    population = coalesce(
      as.numeric(population),
      as.numeric(pop1)
    )
  ) |>
  select(-site, -epicenter_ha, -site_area_ha, -pop1)
## Save 
# all sites
write.csv(samp_epi_sf, file.path(targetdir, 'postclassic_cities.csv'))
# only 62 sites have coordinates
samp_epi_sf[!is.na(samp_epi_sf$latitude),] |>
  st_as_sf(
    coords = c('longitude', 'latitude'),
    crs = st_crs(4326)
  ) |>
  st_write(
    dsn = file.path(targetdir, 'sites.gpkg'),
    layer = 'postclassic',
    append = FALSE
    )