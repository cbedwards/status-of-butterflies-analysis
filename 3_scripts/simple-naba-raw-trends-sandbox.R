library(tidyverse)
library(here)
library(mgcv)
library(sf)
library(sp)
library(geometry)
library(viridis)
library(scales)
library(ggpubr)
library(rgdal)
library(geometry)
library(quantreg)
source(here("3_scripts/funs.R"))
source(here("3_scripts/funs-fitting.R"))
source(here("3_scripts/model-fitting-function.R"))


# gr = .005
# (exp(gr*20)-1)*100
# 
# dat.sum = dat %>% 
#   group_by(source, year) %>% 
#   summarize(nrecord = n()) %>% 
#   ungroup()
# 
# dat.sites = 
#   
#   ggplot(dat.sum %>% filter(source == "Illinois Butterfly Monitoring Network"), aes(x = year, y = nrecord))+
#   geom_path()

cur.code = "EVECOM"

make_dataset(cur.code, use.range = T,
             infer.messy.levels = c("GENUS", "SUBFAMILY", "FAMILY", "COMPLEX"),
             name.pretty = NULL)
dat.spec = read_csv(paste0("2_data_wrangling/cleaned by code/", cur.code, ".csv"))

## Grab FWS regions, designated by state based on this map: https://www.fws.gov/about/regions
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))

## add regions and a factor version of regions to the data  
dat.spec = left_join(dat.spec, regions.dict, by = "state")
### Misc ------------

dat.use = dat.spec %>% 
  filter(source == "ATACHAM") %>% 
  group_by(region, year) %>% 
  summarise(ave = mean(count)) %>% 
  ungroup()

ggplot(dat.use, aes(x = year, y = ave))+
  geom_path()+
  geom_smooth(method = "lm")+
  facet_wrap(.~region)


rq(doy ~ year, tau = c(.1, .5, .9), data = dat.spec %>%  filter(source == "NFJ"))

res = NULL
print(cur.code)
for(cur.region in na.omit(unique(dat.spec$region))){
  dat.fit = dat.spec %>% 
    filter(source == "NFJ") %>% 
    filter(region == cur.region) %>% 
    filter(year>2000)
    # filter(doy > 100) %>% 
    # filter(doy <200)
  out.lm = glm.nb(count ~ year, data = dat.fit)
  res.cur = data.frame(region = cur.region, 
                       trend = coef(out.lm)[2],
                       n = nrow(dat.fit))
  res = rbind(res, res.cur)
}
res

