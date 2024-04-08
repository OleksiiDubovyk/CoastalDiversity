# PREAMBLE ----

# adds timemark to RStudio console statements
h <- taskCallbackManager()
h$add(function(expr, value, ok, visible) {
  options("prompt"=format(Sys.time(), "%H:%M:%S> "));
  return(TRUE) },
  name = "simpleHandler")

# clean the environment
rm(list = ls())

# LIBRARIES ---- 
library(tidyverse)
library(lubridate)

# MAIN ----

## Data import ----

dets <- read_csv("./detections.csv")
depls <- read_csv("./deployments.csv")

# Sunrise/sunset calculations ----

source("suntime.R")

# transform local time into a degree: say midnight is 0=2pi, sunrise should be pi/2, noon should be pi, and sunset should be 3/2 pi

time2deg <- function(time, sunrise, sunset){
  # all args in dttm format
  # time <- ymd_hms(time)
  # sunrise <- ymd_hms(sunrise)
  # sunset <- ymd_hms(sunset)
  daylength <- difftime(sunset, sunrise, units = "s") %>% as.numeric()
  noon <- round(sunrise + daylength/2)
  midnight <- noon - 60*60*12
  d_sunrise <- difftime(sunrise, midnight, units = "s") %>% as.numeric()
  d_sunset <- difftime(sunset, midnight, units = "s") %>% as.numeric()
  d_noon <- 60*60*12
  d_time <- difftime(time, midnight, units = "s") %>% as.numeric()
  if (time >= midnight & time < sunrise){
    return(d_time/d_sunrise * 0.5*pi)
  }else if (time >= sunrise & time < noon){
    return((d_time - d_sunrise)/(d_noon - d_sunrise) * 0.5*pi + 0.5*pi)
  }else if (time >= noon & time < sunset){
    return((d_time - d_noon)/(d_sunset - d_noon) * 0.5*pi + pi)
  }else{
    return((d_time - d_sunset)/(24*60*60 - d_sunset) * 0.5*pi + 1.5*pi)
  }
}

# Time shift dates - when EST UTC-05:00 switched to EDT UTC-04:00 or vice versa

timedate2offset <- function(dttm){
  # function to figure out what the time offset is based on the date
  # will NOT work correctly for any date before 2022-11-05 and after 2024-11-02
  edt2est2021 <- ymd_hm("2021-11-06 02:00")
  est2edt2022 <- ymd_hm("2022-03-12 02:00")
  edt2est2022 <- ymd_hm("2022-11-05 02:00")
  est2edt2023 <- ymd_hm("2023-03-11 02:00")
  edt2est2023 <- ymd_hm("2023-11-04 02:00")
  est2edt2024 <- ymd_hm("2024-03-09 02:00")
  edt2est2024 <- ymd_hm("2024-11-02 02:00")
  if (dttm < edt2est2021) {
    return(NA)
  }else if (dttm > edt2est2021 & dttm <= est2edt2022){
    return(-5)
  }else if (dttm > est2edt2022 & dttm <= edt2est2022){
    return(-4)
  }else if (dttm > edt2est2022 & dttm <= est2edt2023){
    return(-5)
  }else if (dttm > est2edt2023 & dttm <= edt2est2023){
    return(-4)
  }else if (dttm > edt2est2023 & dttm <= est2edt2024){
    return(-5)
  }else if (dttm > est2edt2024 & dttm <= edt2est2024){
    return(-4)
  }else{
    return(NA)
  }
}

# dets_suntimes <- sapply(dets$DateTime, function(x){
#   suntime(date = date(x), 
#           lat = 36.8794, lon = -76.2892, 
#           utc_offset = timedate2offset(x))
# }) %>% t() %>% as_tibble() # takes time, about 10 minutes
# colnames(dets_suntimes) <- c("sunrise", "sunset")
# dets_suntimes <- dets_suntimes %>% mutate(
#   sunrise = ymd_hms(sunrise),
#   sunset = ymd_hms(sunset)
# ) %>%
#   select(sunrise, sunset)
# 
# dets <- bind_cols(dets, dets_suntimes)

# det_sintime <- sapply(1:nrow(dets), function(i){
#   sin(time2deg(time = dets$DateTime[i], sunrise = dets$sunrise[i], sunset = dets$sunset[i]))
# })
# 
# det_costime <- sapply(1:nrow(dets), function(i){
#   cos(time2deg(time = dets$DateTime[i], sunrise = dets$sunrise[i], sunset = dets$sunset[i]))
# })
# 
# dets <- dets %>%
#   mutate(sintime = det_sintime,
#          costime = det_costime)

# 
# write_csv(dets, "./dets_w_suntime.csv")

dets <- read_csv("./dets_w_suntime.csv")

## Getting traits ----

# dets %>%
#   select(Species, Guild) %>%
#   distinct() %>% 
#   write_csv("taxa.csv", col_names = F)

# ^ Irosh added scientific names to taxa.csv, do not run because line 135 will overwrite the file - OD

taxa <- read_csv("taxa.csv", col_names = F) %>% 
  mutate(eng = X1, sci = X3) %>% 
  select(eng, sci) %>% 
  distinct() %>%
  filter(!is.na(sci))

traits_birds <- read_tsv("BirdFuncDat.txt")
traits_mammals <- read_tsv("MamFuncDat.txt")

traits_birds <-  taxa %>% inner_join(traits_birds, by = c("sci" = "Scientific")) # yes, it's not good to overwrite the tibbles, but we're not gonna need them - OD
traits_mammals <-  taxa %>% inner_join(traits_mammals, by = c("sci" = "Scientific"))

traits_all <- bind_rows(
  traits_birds %>%
    select(intersect(colnames(traits_birds), colnames(traits_mammals))) %>%
    mutate(class = "bird"),
  traits_mammals %>%
    select(intersect(colnames(traits_birds), colnames(traits_mammals))) %>%
    mutate(class = "mammal")
) %>% # join birds and mammals with only traits for which data available for both
  mutate(d_inv = `Diet-Inv`, 
         d_vend = `Diet-Vend`,
         d_vect = `Diet-Vect`,
         d_vfish = `Diet-Vfish`,
         d_vunk = `Diet-Vunk`,
         d_scav = `Diet-Scav`,
         d_fru = `Diet-Fruit`,
         d_nect = `Diet-Nect`,
         d_seed = `Diet-Seed`,
         d_plant = `Diet-PlantO`,
         bodymass = `BodyMass-Value`) %>% # rename variables with easier abbrs
  select(eng, sci, d_inv, d_vend, d_vect, d_vfish, d_vunk, d_scav, d_fru, d_nect, d_seed, d_plant, bodymass) # keep the important vars
