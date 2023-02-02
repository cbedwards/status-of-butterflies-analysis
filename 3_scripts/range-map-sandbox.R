# https://eriqande.github.io/rep-res-eeb-2017/plotting-spatial-data-with-ggplot.html

library(here)
library(tidyverse)
library(mgcv)
library(sf)
library(sp)
library(geometry)
library(viridis)
library(scales)
library(ggpubr)
library(geometry)
library(formula.tools)
library(rgdal)
source(here("3_scripts/funs.R"))
source(here("3_scripts/funs-fitting.R"))
source(here("3_scripts/model-fitting-function.R"))
library(ggspatial)
state.map.data <- maps::map('state', fill = TRUE, plot = FALSE) %>%
  sf::st_as_sf()
temp = rgdal::readOGR(here("2_data_wrangling/range-maps/EUPHPHA.kml"))
# temp = st_read(here("2_data_wrangling/range-maps/EUPHPHA.kml"))

make_dataset(code = "EUPHPHA")
dat = qread(paste0("2_data_wrangling/cleaned by code/EUPHPHA.csv"))
dat.plot = dat %>% 
  mutate(present = count>0,
         lat = round(lat, 2),
         lon = round(lon, 2)) %>% 
  select(lat, lon, present) %>% 
  unique()

ggplot()+
  geom_sf(data = state.map.data)+
  layer_spatial(data = temp, fill = 'cornflowerblue', alpha = 0.2)+
  geom_point(data = dat.plot, aes(x = lon, y = lat, col = present), shape=15, alpha = .5)+
  theme_minimal()+
  ggtitle("Baltimore checkerspot butterfly")



dat.sp = SpatialPoints(as.matrix(dat.plot[,c('lat','lon')]),
                       proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
out = over(dat.sp, temp)


dat.sf = st_as_sf(dat.plot, coords = c("lon", "lat"), crs = st_crs(temp))
# out = st_contains(dat.sf$geometry, temp)
# plot(temp)
# 
# library(rgeos)
# gIntersects(dat.sf, temp)
