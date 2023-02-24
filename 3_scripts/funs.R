library(here)
library(sp)
library(rworldmap)
library(spData)
library(tidyverse)
library(data.table, include.only = c("fread", "fwrite"))

#quick read function for ease of reading
library(data.table, include.only = c("fread", "fwrite"))
qread = function(x) as.data.frame(fread(here(x))) 

## ggplot theme for easier visualization

theme.larger =   theme(axis.title = element_text(size = rel(1.8)),
                       axis.text = element_text(size = rel(1.8)),
                       strip.text = element_text(size = rel(1.8)),
                       plot.title = element_text(size = rel(1.8)),
                       legend.text = element_text(size = rel(1.8)),
                       legend.title = element_text(size = rel(1.8)),
)

coords_2country = function(long,#longitude (vector)
                           lat #latitude (vector)
)
  # Function to calculate the country for each pair of longitudes and latitudes.
{  
  points = data.frame(long, lat)
  countriesSP <- getMap(resolution='low')
  #countriesSP <- getMap(resolution='high') #you could use high res map from rworldxtra if you were concerned about detail
  
  # converting points to a SpatialPoints object
  # setting CRS directly to that from rworldmap
  pointsSP = SpatialPoints(points, proj4string=CRS(proj4string(countriesSP)))  
  
  
  # use 'over' to get indices of the Polygons object containing each point 
  indices = over(pointsSP, countriesSP)
  
  #indices$continent   # returns the continent (6 continent model)
  # indices$REGION   # returns the continent (7 continent model)
  # indices$ADMIN  #returns country name
  indices$ISO3 # returns the ISO3 code
}


coords_2state <- function(long,
                          lat,
                          states = spData::us_states,
                          name_col = "NAME") {
  pointsDF = data.frame(long = long, lat = lat)
  ## Convert points data.frame to an sf POINTS object
  pts <- st_as_sf(pointsDF, coords = 1:2, crs = 4326)
  
  ## Transform spatial data to some planar coordinate system
  ## (e.g. Web Mercator) as required for geometric operations
  states <- st_transform(states, crs = 3857)
  pts <- st_transform(pts, crs = 3857)
  
  ## Find names of state (if any) intersected by each point
  state_names <- states[[name_col]]
  ind <- as.integer(st_intersects(pts, states))
  return(state_names[ind])
}


### checking trips
## NOTE: this is ONLY written to work with a single species at a time, and based on 
## using trip-wrangling data generated in final-data-integration.R and stored in 2_data_wrangling/trip-wrangling/

trip_abs = function(dat, #full data set
                    code.cur, #code to work with
                    infer.messy.levels = c("GENUS", "SUBFAMILY", "FAMILY", "COMPLEX")){ #what levels of unidentified to not include in infering zeroes
  dat$inferred = FALSE
  dat.com = dat %>%  
    filter(gather == "community") %>% 
    filter(!is.na(code))
  trip.temp = as.data.frame(fread(here("2_data_wrangling/trip-wrangling/trip-templates.csv")))
  dict.zeroes = as.data.frame(fread(here("2_data_wrangling/trip-wrangling/dictionary-zeroes.csv")))
  dict.zeroes = dict.zeroes[dict.zeroes$code == code.cur,]
  dict.zeroes = dict.zeroes[dict.zeroes$tax.messy.level %in% infer.messy.levels,]
  dat.usecodes = c(code.cur, dict.zeroes$code.messy)
  dat.obs = dat.com[dat.com$code %in% dat.usecodes,]
  events.obs = unique(dat.obs$event.id)
  events.no.obs = unique(dat.com$event.id)
  events.no.obs = events.no.obs[!events.no.obs %in% events.obs]
  cat(paste0("Total of ", length(events.obs), " observations of taxa (or related unknowns), and\n", length(events.no.obs), " inferred zero events.\n"))
  #How many inferred zeroes are we skipping because of our infer.messy.levels?
  dat.obs.exact = dat.com[dat.com$code %in% code.cur,]
  events.obs.exact = unique(dat.obs.exact$event.id)
  events.no.obs.exact = unique(dat.com$event.id)
  events.no.obs.exact = events.no.obs.exact[!events.no.obs.exact %in% events.obs.exact]
  cat(paste0("Because of our current infer.messy.levels, we skipped inferring zeroes for ", 
             length(events.no.obs.exact)-length(events.no.obs),
             " events\n"))
  #identify all the trips that were NOT observed in the data
  trip.abs = trip.temp[trip.temp$event.id %in% events.no.obs,]
  trip.abs$code = code.cur
  trip.abs$name = paste0("absence for ", code.cur)
  trip.abs$inferred = TRUE
  trip.abs$common = dat[dat$code == code.cur,"common"][1]
  trip.abs = trip.abs[,names(dat)]
  dat.use = dat[dat$code == code.cur, ]
  dat.use = dat.use[!is.na(dat.use$code),]
  dat.use = rbind(dat.use, trip.abs)
  return(list(dat = dat.use,
              events.missed.messy = length(events.no.obs.exact)-length(events.no.obs))
  )
}


#creates a dataset from a GU code, adding absences from community gathering
#events that didn't report that species. Use infer.messy.levels to handle zero inference 
#when there is a report of an uknown at genus level, family 
make_dataset = function(code.cur, 
                        use.range,
                        infer.messy.levels,
                        name.pretty = NULL){
  if(is.null(name.pretty)){name.pretty = code.cur}
  dat = as.data.frame(fread(here("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv")))
  dat.trips = trip_abs(dat, code.cur, infer.messy.levels)
  dat.cur = dat.trips$dat
  events.missed.messy = dat.trips$events.missed.messy
  ## add code to restrict to range
  ## just use polygon of range map
  if(use.range){
    print("Constraining data to range")
    dat.cur = use_range(dat = dat.cur, code.cur = code.cur)
  }
  fwrite(dat.cur, 
         here("2_data_wrangling/cleaned by code", paste0(name.pretty,".csv")))
  return(events.missed.messy)
}

use_range = function(dat, #data frame of any data of interest, with columns `lon` and `lat` for longitude and latitude.
                     code.cur){
  dat.sp = dat # separate object for spatial coordinates, because I'm paranoid about interactions with other code #use polygon range map
  range.map = rgdal::readOGR(here(paste0("2_data_wrangling/range-maps/relabeled/", code.cur, ".kml")))
  sp::coordinates(dat.sp) <- ~ lon + lat
  sp::proj4string(dat.sp) <- sp::proj4string(range.map)
  overlaid <- sp::over(dat.sp, as(range.map, "SpatialPolygons"))
  ind.within <- which(!is.na(overlaid))
  dat = dat[ind.within,]
  return(dat)
}

## quickly filter to avoid plotting multiple points with identical lat/lon (makes mapping much faster)
viz_filter = function(dat, reso = 4){
  stopifnot(c("lat","lon") %in% names(dat))
  ind.bad = duplicated(round(dat[,c("lat","lon")], reso))
  return(dat[!ind.bad,])
}

loc_viz = function(dat, reso = 4){
  ## quick viz of map to check for outliers etc
  state.map.data <- map('state', fill = TRUE, plot = FALSE)%>% st_as_sf()
  ggplot() +
    geom_sf(data= state.map.data)+
    geom_point(data = viz_filter(dat, reso = reso),
               aes(x = lon, y = lat), size=1, col = 'blue')
}

## data filtering to summarize by lat/lon rounding, originally for shiny app
dat_filt = function(dat, code.cur, dat.source, rnd.digits = 0){
  ## filter data to the nearest lat/lon (or finer scale for rnd.digits>0)
  ## and summarizer number of years of coverage (nyear), number of observations (nobs)
  ##  [includes reported 0s but not inferred 0s], and number of butterflies reported total (ncount)
  ##  Note that for nyear, we're not counting years in which only 0s were reported.
  if(dat.source %in% c("pollard", "Pollard")){
    dat.cur = dat %>% 
      filter(code == code.cur) %>% 
      filter(source != "NFJ") %>% 
      filter(source != "MASSBfly") %>% 
      mutate(lon.rnd = round(lon, digits= rnd.digits),
             lat.rnd = round(lat, digits= rnd.digits))
  }else{
    dat.cur = dat %>% 
      filter(code == code.cur) %>% 
      filter(source == dat.source) %>% 
      mutate(lon.rnd = round(lon, digits= rnd.digits),
             lat.rnd = round(lat, digits= rnd.digits))
  }
  dat.year = dat.cur %>%
    filter(count > 0) %>% 
    select(year, lon.rnd, lat.rnd) %>% 
    unique() %>% 
    group_by(lon.rnd, lat.rnd) %>% 
    summarize(nyear = n()) %>% 
    ungroup()
  dat.obs = dat.cur %>%
    select(lon.rnd, lat.rnd, date) %>% 
    unique() %>% 
    group_by(lon.rnd, lat.rnd) %>% 
    summarize(nobs = n()) %>% 
    ungroup()
  dat.count = dat.cur %>% 
    group_by(lon.rnd, lat.rnd) %>% 
    summarize(ncount = sum(count)) %>% 
    ungroup()
  dat.sum = inner_join(dat.year, dat.obs, by = c("lon.rnd", "lat.rnd"))
  dat.sum = inner_join(dat.sum, dat.count, by = c("lon.rnd", "lat.rnd"))
  return(dat.sum)
}

## data filtering to summarize by lat/lon rounding for scatterpie, originally for shiny app
dat_filt_pie = function(dat, code.cur, rnd.digits = 0){
  dat %>% 
    filter(code == code.cur) %>%
    # filter(source == dat.source) %>%
    mutate(lon.rnd = round(lon, digits= rnd.digits),
           lat.rnd = round(lat, digits= rnd.digits)) %>% 
    group_by(name, lon.rnd, lat.rnd) %>% 
    summarise(value = sum(count)) %>% #scatterpie is pretty jank, values have to be in column called "value"
    ungroup()
}













