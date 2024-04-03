# (c) Oleksii Dubovyk, oadubovyk@gmail.com
# Dept of Biological Sciences, Old Dominion University, Norfolk, VA 23529, USA
# Apr 2, 2024

suntime <- function(date, lat, lon, utc_offset){
  # based on https://gml.noaa.gov/grad/solcalc/calcdetails.html
  
  # date as 'yyyy-mm-dd'
  
  #functions for radians/degrees transformation
  rad2deg <- function(rad) {(rad * 180) / (pi)}
  deg2rad <- function(deg) {(deg * pi) / (180)}
  
  #function for clock time
  s2hms <- function(secs){
    td <- seconds_to_period(secs)
    sprintf('%02d:%02d:%02d', td@hour, minute(td), second(td))
  }
  
  date <- ymd(date)
  
  tsm <- seq(0.1, 24, 0.1)/24 # time since midnight
  
  jd <- (difftime(ymd(date), ymd("1899-12-30"), units = "d") %>% as.numeric()) + 2415018.5 + tsm - utc_offset/24 # julian date
  jc <- (jd - 2451545)/36525 # julian century
  
  gmlsd <- (280.46646 + jc * (36000.76983 + jc*0.0003032)) %% 360 # geom mean long sun (degrees)
  gmasd <- 357.52911 + jc*(35999.05029 - 0.0001537*jc) # geom mean anom sun (degrees)
  
  eeo <- 0.016708634 - jc*(0.000042037+0.0000001267*jc)# eccent Earth orbit
  
  sec <- sin(deg2rad(gmasd))*(1.914602 - jc*(0.004817 + 0.000014*jc)) + sin(deg2rad(2*gmasd))*(0.019993 - 0.000101*jc) + sin(deg2rad(3*gmasd))*0.000289 # Sun Eqn of center
  
  stl <- gmlsd + sec # Sun true lon
  sta <- gmasd + sec # Sun true anom
  
  srv <- (1.000001018*(1 - eeo^2))/(1+eeo*cos(deg2rad(sta)))# Sun rad vector (AUs)
  sal <- stl-0.00569-0.00478*sin(deg2rad(125.04 - 1934.136*jc)) # Sun approx lon (deg)
  
  moe <- 23+(26+((21.448-jc*(46.815+jc*(0.00059-jc*0.001813))))/60)/60# mean obliq ecliptic (deg)
  oc <- moe+0.00256*cos(deg2rad(125.04-1934.136*jc)) # obliq corr (deg)
  
  sra <- rad2deg(atan2(cos(deg2rad(oc))*sin(deg2rad(sal)), cos(deg2rad(sal)))) # Sun Rt ascen (deg)
  sdecl <- rad2deg(asin(sin(deg2rad(oc))*sin(deg2rad(sal)))) # Sun declin
  
  y <- tan(deg2rad(oc/2))^2
  
  eot <- 4*rad2deg(y*sin(2*deg2rad(gmlsd))-2*eeo*sin(deg2rad(gmasd))+4*eeo*y*sin(deg2rad(gmasd))*cos(2*deg2rad(gmlsd))-0.5*y*y*sin(4*deg2rad(gmlsd))-1.25*eeo*eeo*sin(2*deg2rad(gmasd))) # Eqn of time
  has <- rad2deg(acos(cos(deg2rad(90.833))/(cos(deg2rad(lat))*cos(deg2rad(sdecl)))-tan(deg2rad(lat))*tan(deg2rad(sdecl))))# HA sunrise (deg)
  
  solar_noon <- (720-4*lon-eot+utc_offset*60)/1440
  # solar_noon_time <- sapply(solar_noon, function(x) s2hms(round(x*24*60*60)))
  
  solar_dif <- has*4/1440
  sunrise <- solar_noon - solar_dif
  sunset <- solar_noon + solar_dif
  
  sunrise_time <- s2hms(mean(sunrise*60*60*24) %>% round())
  sunset_time <- s2hms(mean(sunset*60*60*24) %>% round())
  
  return(c(paste(date, sunrise_time), paste(date, sunset_time)))
  
}
