## For each species code, find the square kilometer of species range in each region.
## Note that this is approximate, as it is based on the range blobs, but serves our general purpose
## of relative weighting based on area.

library(here)
library(tidyverse)
library(sf)
library(terra)

## load in FWS regional maps, also the list of region names used in the other data analysis from the csv.
usfws = st_read(here("2_data_wrangling/FWS_National_Legacy_Regional_Boundaries/FWS_National_Legacy_Regional_Boundaries.shp"))
usfws = terra::vect(usfws)

## identify all codes with range maps present:
maps.list = list.files(here("2_data_wrangling/range-maps/relabeled/"))
maps.list = grep("*.kml", maps.list, value = TRUE)
maps.list = gsub(".kml", "", maps.list)


## make storage object
res.area = NULL
##loop over each code

for(cur.code in maps.list){
  range.map = st_read(here(paste0("2_data_wrangling/range-maps/relabeled/", cur.code, ".kml")))
  range.map = terra::vect(range.map)
  
  
  cur.area = data.frame(region = unique(usfws$REGNAME))
  cur.area$km2 = -99
  cur.area$code = cur.code
  for(i in 1:nrow(cur.area)){
    overlap = terra::intersect(range.map, usfws[usfws$REGNAME==cur.area$region[i],])
    if(length(overlap)>0){
      cur.area$km2[i] = terra::expanse(overlap, unit = "km")
    }else{
      cur.area$km2[i] = 0
    }
  }
  res.area = rbind(res.area, cur.area)
}

## clear out entries with no abundance, the alaska region that we're not using
res.area = res.area %>% 
  filter(region != "Alaska Region") %>% 
  filter(km2>0)


## fix naming conventions - make it match other analysis
## For identifying mismatches and fixes, I read in the region-by-state file. 
## This is not needed now that I have the code to fix itcode. 
# regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))
# dat.regionnames = unique(regions.dict$region)

res.area$region = gsub(" Region", "", res.area$region)
res.area$region = gsub("Mountain Prairie", "Mountain-Prarie", res.area$region)

write_csv(res.area,
          here("2_data_wrangling/range-area-by-regions.csv"))

