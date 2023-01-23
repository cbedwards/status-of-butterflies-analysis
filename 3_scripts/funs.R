library(here)
library(sp)
library(rworldmap)
library(spData)

#quick read function for ease of reading
library(data.table, include.only = c("fread", "fwrite"))
qread = function(x) as.data.frame(fread(here(x))) 


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
library(data.table, include.only = c("fread", "fwrite"))
trip_abs = function(dat){
  trip.temp = as.data.frame(fread(here("2_data_wrangling/trip-wrangling/trip-templates.csv")))
  usecodes = as.data.frame(fread(here("2_data_wrangling/trip-wrangling/trip-codes.csv")))
  dat$inferred = FALSE
  dat.usecode = unique(dat$code)[unique(dat$code) %in% usecodes$x]
  if(length(dat.usecode>0)){
    print("integrating trip absences")
    for(i in 1:length(dat.usecode)){
      #identify all the trips that were NOT observed in the data
      code.cur = dat.usecode[i]
      dat.cur = dat[dat$code == code.cur,]
      trip.abs = trip.temp[!(trip.temp$event.id %in% dat.cur$event.id),]
      trip.abs$code = code.cur
      trip.abs$name = paste0("absence for ", code.cur)
      trip.abs$inferred = TRUE
      trip.abs$common = dat[dat$code == code.cur,"common"][1]
      trip.abs = trip.abs[,names(dat)]
      dat = rbind(dat,trip.abs)
    }
  }else{
    print("taxa was never present for community sampling, skipping adding absences") 
  }
  return(dat)
}


#creates a dataset from 1 or more GU codes, adding absences from community gathering
#events that didn't report that species. Can give multiple GU codes
make_dataset = function(code, name.pretty = NULL){
  if(is.null(name.pretty)){name.pretty = code}
  dat = as.data.frame(fread(here("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv")))
  dat.cur = trip_abs(dat[dat$code %in% code,])
  fwrite(dat.cur, 
         here("2_data_wrangling/cleaned by code", paste0(name.pretty,".csv")))
}


## quickly filter to avoid plotting multiple points with identical lat/lon (makes mapping much faster)
viz_filter = function(dat, reso = 4){
  stopifnot(c("lat","lon") %in% names(dat))
  ind.bad = duplicated(round(dat[,c("lat","lon")], reso))
  return(dat[ind.bad,])
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













