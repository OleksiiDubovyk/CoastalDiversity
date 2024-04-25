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
library(ggplot2)

# MAIN ----

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Data import ---------------------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dets <- read_csv("./detections.csv") %>%
  filter(Species %!in% c("NEED NEW LABEL", "Unknown"))
depls <- read_csv("./deployments.csv") %>% filter(total_minutes_final > 0)

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
#                                                                                       \/
#                                                                                       \/
dates_dets <- tibble(
  date = dets %>% .$DateTime %>% date() %>% unique()
)

dets_suntimes <- sapply(dates_dets$date, function(x){
  suntime(date = x,
          lat = 36.8794, lon = -76.2892,
          utc_offset = timedate2offset(x))
}) %>% t() %>% as_tibble()
colnames(dets_suntimes) <- c("sunrise", "sunset")
dets_suntimes <- dets_suntimes %>% mutate(
  sunrise = ymd_hms(sunrise),
  sunset = ymd_hms(sunset)
) %>%
  select(sunrise, sunset)

dates_dets <- bind_cols(dates_dets, dets_suntimes)
remove(dets_suntimes)

dets <- dets %>%
  mutate(date = date(DateTime)) %>%
  left_join(dates_dets, by = "date") %>%
  select(-date)

det_degtime <- sapply(1:nrow(dets), function(i){
  time2deg(time = dets$DateTime[i],
           sunrise = dets$sunrise[i],
           sunset = dets$sunset[i]) %% (2*pi) # %% to ensure that deg is between 0 and 2pi
})

dets <- dets %>%
  mutate(degtime = det_degtime)

det_sintime <- sapply(1:nrow(dets), function(i){
  sin(dets$degtime[i])
})

det_costime <- sapply(1:nrow(dets), function(i){
  cos(dets$degtime[i])
})

dets <- dets %>%
  mutate(sintime = det_sintime,
         costime = det_costime)

remove(det_costime, det_sintime, det_degtime)
# 
# write_csv(dets, "./dets_w_suntime.csv")
# dets <- read_csv("./dets_w_suntime.csv")
#                                                                                       /\
#                                                                                       /\
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# timing convention is the following:
# 0 = astronomic midnight
# pi/2 = sunrise
# pi = astronomic noon
# 3*pi/2 = sunset
# 2*pi = astronomic midnight

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Getting traits ------------------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

remove(traits_birds, traits_mammals)

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
  mutate(enddate = ymd_hms(paste(paste(y, m, d, sep = "-"), "12:00:00"))) %>%
  select(Site, Type, Deployment, total_minutes_final, enddate) %>%
  mutate(startdate = enddate - total_minutes_final*60 + 60)

dets_empty <- tibble()

for (i in 1:nrow(depls_dated)){
  var_times <- depls_dated$startdate[i] + seq(from = 0, to = (-1 + depls_dated$total_minutes_final[i])*60, by = 60)
  
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
  left_join((tides %>%
               group_by(DateTime) %>%
               summarise(level = min(level), Tide = Tide[1])), c("DateHour" = "DateTime")) %>%
  select(Site, Deployment, Type, DateTime, Species, level, Tide)

# now we have an index of all pictures ever taken and their tides (plus/minus 30 min)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Degtimes of nondetections -------------------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dates_nondets <- tibble(
  date = dets_empty %>% .$DateTime %>% date() %>% unique()
)

nondets_suntimes <- sapply(dates_nondets$date, function(x){
  suntime(date = x,
          lat = 36.8794, lon = -76.2892,
          utc_offset = timedate2offset(x))
}) %>% t() %>% as_tibble()
colnames(nondets_suntimes) <- c("sunrise", "sunset")
nondets_suntimes <- nondets_suntimes %>% mutate(
  sunrise = ymd_hms(sunrise),
  sunset = ymd_hms(sunset)
) %>%
  select(sunrise, sunset)

dates_nondets <- bind_cols(dates_nondets, nondets_suntimes)
remove(nondets_suntimes)

dets_empty <- dets_empty %>%
  mutate(date = date(DateTime)) %>%
  left_join(dates_nondets, by = "date") %>%
  select(-date)

nondet_degtime <- sapply(1:nrow(dets_empty), function(i){
  time2deg(time = dets_empty$DateTime[i],
           sunrise = dets_empty$sunrise[i],
           sunset = dets_empty$sunset[i]) %% (2*pi) # %% to ensure that deg is between 0 and 2pi
}) # takes 30ish mins

dets_empty <- dets_empty %>%
  mutate(degtime = nondet_degtime)

# dets_empty %>% write_csv("dets_empty.csv")
# 
# dets_empty <- read_csv("dets_empty.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# End of Degtimes of nondetection ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### ED trying something here, but doesn't account for frequency of all tide levels

dets_tides <- dets
dets_tides$DateTime <- as.POSIXct(dets$DateTime)

# Function to find the closest tide level for each datetime
find_closest_tide <- function(datetime) {
  closest_index <- which.min(abs(as.numeric(datetime - tides$DateTime)))
  return(tides$level[closest_index])
}

# Add a new column level2 to det_tides with the closest tide level for each DateTime
dets_tides$level2 <- sapply(dets_tides$DateTime, find_closest_tide)

# All Dets vs Tides
ggplot(dets_tides, aes(x = level2)) +
  geom_histogram(binwidth = .25, color = "black", fill = "skyblue") +
  labs(title = "Number of Detections by Tide Level",
       x = "Tide Level",
       y = "Count")

#Bird dets vs Tides
birds <- dets_tides %>%
  filter(Guild %in% c("Perching Bird", "Wading Bird", "Marsh Bird", "Raptor", "Seabird", "Shorebird", "Waterfowl"))

ggplot() +
  geom_histogram(data = birds, aes(x = level2), binwidth = 0.25, color = "black", fill = "skyblue") +
  labs(title = "Number of Bird Detections by Tide Level",
       x = "Tide Level",
       y = "Count") +
  ggtitle("Number of Bird Detections by Tide Level")

#Mammal dets vs Tides
mammals <- dets_tides %>%
  filter(Guild %in% c("Mesomammal", "Small Mammal", "Large Mammal"))

ggplot() +
  geom_histogram(data = mammals, aes(x = level2), binwidth = 0.25, color = "black", fill = "skyblue") +
  labs(title = "Number of Mammal Detections by Tide Level",
       x = "Tide Level",
       y = "Count") +
  ggtitle("Number of Mammal Detections by Tide Level")

#Using categorical tide levels

SpeciesCount <- dets %>%
  filter(Species != "Unknown") %>%  #NOTE ignoring unknowns
  group_by(Type, Species) %>%
  summarize(Total_Observations = n()) %>%
  pivot_wider(names_from = Type, values_from = Total_Observations, names_prefix = "Count_") %>%
  mutate(Totals = rowSums(select(., starts_with("Count")), na.rm = TRUE))


# Filter only species with n >= 5
SpeciesCount_filtered <- SpeciesCount %>%
  filter(Totals >= 5)

# Reorder levels of Tide variable
dets$Tide <- factor(dets$Tide, levels = c("High", "Mid", "Low"))

#Plot
ggplot(dets, aes(x = Guild, fill = Tide)) +
  geom_bar(position = "fill", width = 0.95) +  # Adjust width of bars and space between bars
  labs(title = "Detections Across Tide Levels",
       x = "Guild",
       y = "Proportion") +
  scale_fill_manual(values = c("High" = "darkblue", "Mid" = "blue", "Low" = "lightblue")) +  # Gradient color scale
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),  # Angle x-axis labels and adjust text size
        plot.margin = margin(l = 2, r = 2, unit = "cm")) +  # Adjust plot margin
  coord_flip() +  # Flip coordinates to make the plot wider
  theme(plot.title = element_text(hjust = 0.5))  # Center plot title

# OD: Ella's code doesn't account for the baseline tide level

# Chi-squared test for observations vs null tide levels

null_tides <- depls_dated %>%
  summarise(High = sum(High), Mid = sum(Mid), Low = sum(Low)) %>% unlist()

birds_tides <- dets %>%
  filter(Guild %in% c("Perching Bird", "Wading Bird", "Marsh Bird", "Raptor", "Seabird", "Shorebird", "Waterfowl")) %>%
  .$Tide %>% summary()

mammals_tides <- dets %>%
  filter(Guild %in% c("Mesomammal", "Small Mammal", "Large Mammal")) %>%
  .$Tide %>% summary()

tibble(null_tides, birds_tides) %>% chisq.test()
tibble(null_tides, mammals_tides) %>% chisq.test()



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Rarefaction per location/tide/time bin ------------------------------------------------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Master DataSet

mdets <- dets %>%
  select(colnames(dets)[colnames(dets) %in% colnames(dets_empty)]) %>%
  mutate(DateTime = floor_date(DateTime, "minute")) %>%
  mutate(uniq = paste(Site, Type, DateTime))

mndets <- dets_empty %>%
  mutate(DateTime = floor_date(DateTime, "minute")) %>%
  mutate(uniq = paste(Site, Type, DateTime))

mndets <- mndets %>%
  filter(uniq %!in% mdets$uniq)

mds <- bind_rows(
  mdets, mndets
) %>%
  select(- uniq) %>%
  mutate(timebin = degtime %/% (2*pi/24))

remove(mdets, mndets)

# Rarefaction

permestSR <- function(loc, tide, time, data = mds, n = 0.05, nperm = 100, species = T){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Permutational estimate of observed species richness
  #
  # Args:
  # loc - location code(s), vector if multiple
  # tide - tide level, either "Low", "Mid", or "High", vector if multiple
  # time - timebin, any integer from 0 to 24, vector if multiple
  # data - dataset with both detections and non-detections
  # n - level at which to estimate diversity, where 0 = no observations, 1 = all observations
  # nperm - number of permutations
  #
  # Output:
  # Vector, estimated species richness and 95% CI
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  locs <- unique(data$Site)
  tides <- unique(data$Tide)
  times <- unique(data$timebin)
  
  if (!missing(loc)){
    locs <- loc
  }
  if (!missing(tide)){
    tides <- tide
  }
  if (!missing(time)){
    times <- time
  }
  
  sr <- numeric(nperm)
  
  data <- data %>%
    filter(Site %in% locs,
           Tide %in% tides,
           timebin %in% times)
  
  n <- round(nrow(data) * n)
  
  # pb <- progress_bar$new(
  #   format = "[:bar] :percent eta: :eta elapsed: :elapsed :current step",
  #   clear = FALSE, total = nperm, width = 100)
  for (i in 1:nperm){
    index <- sample(1:nrow(data), size = n, replace = F)
    iter_data <- data[index,]
    if (species){
      sr[i] <- iter_data %>%
        filter(Species != "none") %>%
        .$Species %>%
        unique() %>%
        length()
    }else{
      sr[i] <- iter_data %>%
        filter(Species != "none") %>%
        nrow()
    }
    # pb$tick()
  }
  
  return(c("n" = n, 
           "Est" = mean(sr), 
           "CI_0.025" = unname(quantile(sr, 0.025)), 
           "CI_0.975" = unname(quantile(sr, 0.975))))
  # return(sr)
  
}

permrarSR <- function(nperm = 100, step = 0.05, ...){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Permutational estimate of observed species richness
  #
  # Args:
  # loc - location code
  # tide - tide level, either "Low", "Mid", or "High"
  # time - timebin, any integer from 0 to 24
  # data - dataset with both detections and non-detections
  # nperm - number of permutations 
  # step - step size of sample size, where 0 = no observations, 1 = all observations
  #
  # Output:
  # Tibble with the columns representing n, mean species richness, and 95% CI
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  steps <- seq(0, 1, step)
  out <- sapply(steps, function(n){
    permestSR(n = n, nperm = nperm, ...)
  }) %>% t() %>% as_tibble()
  return(out)
}

testrar <- permrarSR(step = 0.01)

testrar %>%
  ggplot(aes(x = n)) +
  geom_line(aes(y = CI_0.025), color = "gray") +
  geom_line(aes(y = CI_0.975), color = "gray") +
  geom_line(aes(y = Est), color = "black") +
  labs(x = "Sample size", y = "Species richness")

# permutational approach seems to be way to slow for the large datasets (1.5 hrs for 2.6 mln rows)
# there might be a way to develop the analytical approach

# what is the relationship of # individuals sampled vs # photos?

testrarind <- tibble(n = runif(1000, 0, nrow(mds)))
testrarind$ind <- sapply(testrarind$n, function(x){
  mds[sample(1:nrow(mds), x, F),] %>%
    filter(Species != "none") %>%
    nrow()
})

testrarind %>%
  ggplot(aes(x = n, y = ind)) +
  geom_point()

# almost perfectly linear!

lm(ind ~ n, testrarind) %>% summary()

# the line eqn is ind = 1.189e+00 + 6.978e-03 * n

n2ind <- function(n) ifelse(n > 0, round(1.189e+00 + 6.978e-03 * n), 0)
ind2n <- function(ind) ifelse(ind > 0, round((ind - 1.189e+00) / 6.978e-03), 0)

