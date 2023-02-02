## prototype script for identifying species that are too rare for our general analysis
##   approach, but have sufficiently rich timeseries at individual sites that we 
##   can potentially fit trends at a site-by-site basis. 

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
source(here("3_scripts/funs.R"))
source(here("3_scripts/funs-fitting.R"))
source(here("3_scripts/model-fitting-function.R"))
dat = qread(paste0("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv"))


## using GU code

## identify GU codes of interest:
threshold = 200 #min number of observations for "normal" approach
code.tab = table(dat$code)
gu.code.interest = names(code.tab[code.tab < threshold])

## let's look at the first one
dat.cur = dat %>% 
  filter(code == gu.code.interest[1])

cat(paste0("total observations: ", nrow(dat.cur), "\n"))
sites.best = rev(sort(table(dat.cur$site)))[1:min(5,length(unique(dat$site)))]
cat("observations per site for the most data-rich sites: ")
print(sites.best)


## look at data coverage across years for the best sites
dat.sitelevel = dat %>% 
  filter(site %in% names(sites.best)) %>% 
  group_by(site) %>% 
  summarize(number.of.years = length(unique(year)))
cat("number of years covered per site for the most data-rich sites: ")
print(dat.sitelevel)


## plot the site with the most data
ggplot(dat.cur %>% filter(site == names(sites.best)[1]),
       aes(x = doy, y = count))+
  geom_point(size = 1.5) +
  facet_wrap(. ~ year) +
  ggtitle(names(sites.best)[[1]])

## using scientific name (e.g. for subspecies)

# example: red-spotted purple
sciname = "Limenitis arthemis astyanax"

dat.cur = dat %>% 
  filter(name == sciname)

cat(paste0("total observations: ", nrow(dat.cur), "\n"))
sites.best = rev(sort(table(dat.cur$site)))[1:min(5,length(unique(dat$site)))]
cat("observations per site for the most common sites:\n ")
print(sites.best)

dat.sitelevel = dat %>% 
  filter(site %in% names(sites.best)) %>% 
  group_by(site) %>% 
  summarize(number.of.years = length(unique(year)))
print(dat.sitelevel)

## look at data coverage across years for the best sites
dat.sitelevel = dat %>% 
  filter(site %in% names(sites.best)) %>% 
  group_by(site) %>% 
  summarize(number.of.years = length(unique(year)))
cat("number of years covered per site for the most data-rich sites: ")
print(dat.sitelevel)


## plot the site with the most data
ggplot(dat.cur %>% filter(site == names(sites.best)[1]),
       aes(x = doy, y = count))+
  geom_point(size = 1.5) +
  facet_wrap(. ~ year) +
  ggtitle(names(sites.best)[[1]])

