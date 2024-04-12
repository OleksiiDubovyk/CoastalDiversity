# PREAMBLE ----

# adds timemark to RStudio console statements
h <- taskCallbackManager()
h$add(function(expr, value, ok, visible) {
  options("prompt"=format(Sys.time(), "%H:%M:%S> "));
  return(TRUE) },
  name = "simpleHandler")

# clean the environment
rm(list = ls())

# handy function, opposite of %in%

'%!in%' <- function(x,y)!('%in%'(x,y))

# LIBRARIES ---- 
library(data.table)
library(tidyverse)
library(lubridate)
library(caret)
library(progress)

# MAIN ----

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Data import ---------------------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dets <- read_csv("./detections.csv") %>%
  filter(Species %!in% c("NEED NEW LABEL", "Unknown"))
depls <- read_csv("./deployments.csv")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sunrise/sunset calculations ------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

source("suntime.R")

# transform local time into a degree: say midnight is 0=2pi, sunrise should be pi/2, noon should be pi, and sunset should be 3/2 pi

time2deg <- function(time, sunrise, sunset){
  # all args in dttm format
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

# Code used to generate the final version of dets_w_suntime.csv ~~~~~~~~~~~~~~~~~~~~~~~~~~
#
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
# remove(dets_suntimes)
# 
# det_degtime <- sapply(1:nrow(dets), function(i){
#   time2deg(time = dets$DateTime[i], 
#            sunrise = dets$sunrise[i], 
#            sunset = dets$sunset[i]) %% (2*pi) # %% to ensure that deg is between 0 and 2pi
# })
# # 
# dets <- dets %>%
#   mutate(degtime = det_degtime)
# 
# det_sintime <- sapply(1:nrow(dets), function(i){
#   sin(dets$degtime[i])
# })
# 
# det_costime <- sapply(1:nrow(dets), function(i){
#   cos(dets$degtime[i])
# })
# 
# dets <- dets %>%
#   mutate(sintime = det_sintime,
#          costime = det_costime)
# 
# remove(det_costime, det_sintime, det_degtime)
# 
# write_csv(dets, "./dets_w_suntime.csv")
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dets <- read_csv("./dets_w_suntime.csv")

# timing convention is the following:
# 0 = astronomic midnight
# pi/2 = sunrise
# pi = astronomic noon
# 3*pi/2 = sunset
# 2*pi = astronomic midnight

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Getting traits ------------------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# dets %>%
#   select(Species, Guild) %>%
#   distinct() %>% 
#   write_csv("taxa.csv", col_names = F)

# ^ Irosh added scientific names to taxa.csv, do not run because line 142 will overwrite the file - OD

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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Timing of observations -----------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Bird observations as a function of time
birds <- dets %>%
  filter(Guild %in% c("Perching Bird", "Wading Bird", "Marsh Bird", "Raptor", "Seabird", "Shorebird", "Waterfowl"))

ggplot(birds, aes(x = degtime)) +
  geom_histogram(color = "black", bins = 24, boundary = 0) +
  coord_polar("x", start = pi, direction = 1) +
  scale_x_continuous(limits = c(0, 2*pi), breaks = seq(0, 1.5*pi, 0.5*pi), 
                     labels = c("Midnight", "Sunrise", "Noon", "Sunset")) +
  labs(title = "Bird Observations vs. Time") +
  xlab("") +
  ylab("Count")

#Mammal observations as a function of time
mammals <- dets %>%
  filter(Guild %in% c("Mesomammal", "Small Mammal", "Large Mammal"))

ggplot(mammals, aes(x = degtime)) +
  geom_histogram(color = "black", bins = 24, boundary = 0) +
  coord_polar("x", start = pi, direction = 1) +
  scale_x_continuous(limits = c(0, 2*pi), breaks = seq(0, 1.5*pi, 0.5*pi), 
                     labels = c("Midnight", "Sunrise", "Noon", "Sunset")) +
  labs(title = "Mammal Observations vs. Time") +
  xlab("") +
  ylab("Count")

# All observations as a function of time
dets %>%
  ggplot(aes(x = degtime)) +
  geom_histogram(color = "black", bins = 24, boundary = 0) +
  coord_polar("x", start = pi, direction = 1) +
  scale_x_continuous(limits = c(0, 2*pi), breaks = seq(0, 1.5*pi, 0.5*pi), 
                     labels = c("Midnight", "Sunrise", "Noon", "Sunset")) +
  labs(title = "All observations vs. Time") +
  xlab("") +
  ylab("Count")

# Observations of species as a function of time

for (s in dets$Species %>% unique()){
  g <- dets %>%
    filter(Species == s) %>%
    ggplot(aes(x = degtime)) +
    geom_histogram(color = "black", bins = 24, boundary = 0) +
    coord_polar("x", start = pi, direction = 1) +
    scale_x_continuous(limits = c(0, 2*pi), breaks = seq(0, 1.5*pi, 0.5*pi), 
                       labels = c("Midnight", "Sunrise", "Noon", "Sunset")) +
    labs(title = paste(s, "vs. Time")) +
    xlab("") +
    ylab("Count")
  print(g)
}

# Group observations by time bins, species, and locations

dets %>%
  mutate(timebin = degtime %/% (2*pi/24)) %>%
  group_by(Site, Species, Tide, timebin) %>%
  summarise(n = sum(Count), .groups = "drop")

# Get species richness per time bin

dets %>%
  mutate(timebin = degtime %/% (2*pi/24)) %>%
  group_by(Site, Species, Tide, timebin) %>%
  summarise(n = sum(Count), .groups = "drop") %>%
  group_by(Site, Tide, timebin) %>%
  summarise(sprich = n(), .groups = "drop") %>%
  # ggplot part below
  ggplot(aes(x = timebin, y = sprich)) +
  geom_jitter() +
  geom_smooth() + 
  geom_vline(xintercept = c(0, 6, 12, 18, 24)) +
  xlab("Time bin") +
  ylab("Species richness")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Tides ---------------------------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tides <- read_csv("tides.csv")

# next step takes time, ~ 5 min
# tides <- tides %>%
#   mutate(DateTime = paste(Date, `Time (GMT)`) %>% ymd_hms(),
#          level = `Verified (ft)`) %>%
#   select(DateTime, level) %>%
#   mutate(off = sapply(.$DateTime, function(x) timedate2offset(x))) %>% # getting offset from UTC
#   mutate(DateTime = DateTime + off*60*60) %>%# adjust to local time
#   select(DateTime, level)
# write_csv(tides, "tides.csv")

# obs_tides <- dets %>%
#   select(DateTime, Tide) %>%
#   setDT()
# 
# setDT(tides)
# 
# obs_tides[, level := tides[.SD, on = "DateTime", roll = "nearest", level]]
# 
# obs_tides <- as_tibble(obs_tides)
# 
# # K-nearest neighbors classifier
# validationIndex <- createDataPartition(obs_tides$Tide, p = 0.70, list = FALSE)
# train <- obs_tides[validationIndex,]
# test <- obs_tides[-validationIndex,]
# trainControl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
# metric <- "Accuracy"
# fit.knn <- train(Tide ~ level, data = train, method = "knn",
#                  metric = metric, trControl = trainControl)
# knn.k1 <- fit.knn$bestTune
# plot(fit.knn) # K = 5 seems the best
# prediction <- predict(fit.knn, newdata = test)
# 
# tides$Tide <- predict(fit.knn, newdata = tides)
# tides <- tides %>% as_tibble()
# write_csv(tides, "tides.csv")

depls_dated <- depls %>%
  mutate(ymdstr = sapply(Deployment, function(x){
    strsplit(x, split = "_")[[1]][1]
  })) %>%
  mutate(y = substring(ymdstr, 1, 4),
         m = substring(ymdstr, 5, 6),
         d = substring(ymdstr, 7, 8)) %>%
  mutate(startdate = ymd_hms(paste(paste(y, m, d, sep = "-"), "00:00:00"))) %>%
  select(Site, Type, Deployment, total_minutes_final, startdate) %>%
  mutate(enddate = startdate + total_minutes_final*60)

find_tide_distr <- function(start, stop){
  if (start == stop){
    tibble(High = 0, Low = 0, Mid = 0)
  }else{
    tides %>%
      filter(DateTime >= start & DateTime < stop) %>%
      group_by(Tide) %>%
      summarise(n = n()) %>% 
      pivot_wider(names_from = Tide, values_from = n)
  }
}

depls_tide <- lapply(1:nrow(depls_dated), function(i) {
  find_tide_distr(depls_dated$startdate[i], depls_dated$enddate[i])
}) %>% bind_rows()

# this tibble has start, end date of deployments, total amount pics taken, and ratio of hrs of different stages of tide
depls_dated <- bind_cols(depls_dated, depls_tide)

### on the second thought, I could approach it differently. too lazy to clean it up tho

dets_empty <- tibble()

for (i in 1:nrow(depls_dated)){
  var_times <- depls_dated$startdate[i] + seq(0, depls_dated$total_minutes_final[i]*60, 60)
  
  var_empty_dets <- tibble(
    Site = depls_dated$Site[i],
    Deployment = depls_dated$Deployment[i],
    Type = depls_dated$Type[i],
    DateTime = var_times,
    Species = "none"
  )
  
  dets_empty <- bind_rows(dets_empty, var_empty_dets)
}

tides$Tide <- as_factor(tides$Tide)

dets_empty <- dets_empty %>% 
  mutate(DateHour = round_date(DateTime, unit = "h")) %>%
  left_join(tides, c("DateHour" = "DateTime")) %>%
  select(Site, Deployment, Type, DateTime, Species, level, Tide)

# now we have an index of all pictures ever taken and their tides (plus/minus 30 min)

# do the same as with detections - get solar standardized time (will take loooonger - like 2 days)

pb <- progress_bar$new(
  format = "[:bar] :percent eta: :eta elapsed: :elapsed :current step",
  clear = FALSE, total = nrow(dets_empty), width = 100)

dets_empty_suntimes <- sapply(dets_empty$DateTime, function(x){
  pb$tick()
  suntime(date = date(format(x, format="%Y-%m%-%d %T")),
          lat = 36.8794, lon = -76.2892,
          utc_offset = timedate2offset(format(x, format="%Y-%m%-%d %T")))
}) %>% t() %>% as_tibble() # takes time
colnames(dets_empty_suntimes) <- c("sunrise", "sunset")
dets_empty_suntimes <- dets_empty_suntimes %>% mutate(
  sunrise = ymd_hms(sunrise),
  sunset = ymd_hms(sunset)
) %>%
  select(sunrise, sunset)

dets_empty <- bind_cols(dets_empty, dets_empty_suntimes)
remove(dets_empty_suntimes)

det_empty_degtime <- sapply(1:nrow(dets_empty), function(i){
  time2deg(time = dets_empty$DateTime[i],
           sunrise = dets_empty$sunrise[i],
           sunset = dets_empty$sunset[i]) %% (2*pi) # %% to ensure that deg is between 0 and 2pi
})

dets_empty <- dets_empty %>%
  mutate(degtime = det_empty_degtime)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Rarefaction per location/tide/time bin ------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



