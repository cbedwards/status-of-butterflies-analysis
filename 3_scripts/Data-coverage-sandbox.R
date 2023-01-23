## script for exploring issues with zeroes, effort through time, data coverage through time
## for analysis team meeting on 1/23/23


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


## Checking "inferred" zeroes ------------

## read in all data
dat = qread(paste0("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv"))

## Each data collection event is -- hypothetically -- represented with a unique event.id
##  so each NABA great circle per year, each pollard data collection site-day, etc.
length(unique(dat$event.id))   ##~112,308 unique ids. So we should infer at most 112308 
## zeros for a given species.

## make + read inferred zeroes for a single species. Let's start with cabbage whites
code.cur = "PIERAP"
make_dataset(code = code.cur)
dat.cur = qread(paste0("2_data_wrangling/cleaned by code/",code.cur, ".csv"))

## How many inferred zeroes?
sum(dat.cur$inferred) ##We HAD been getting !250k. Now we're getting 27347. This is good.

## Check again for rare species: frosted elfin
code.cur = "CALLIRU"
make_dataset(code = code.cur)
dat.cur = qread(paste0("2_data_wrangling/cleaned by code/",code.cur, ".csv"))

## How many inferred zeroes?
sum(dat.cur$inferred) ## close to the number of unique IDs. So that's good.
length(unique(dat.cur$event.id))

## Plotting theme

theme.larger =   theme(axis.title = element_text(size = rel(1.8)),
                       axis.text = element_text(size = rel(1.8)),
                       strip.text = element_text(size = rel(1.8)),
                       plot.title = element_text(size = rel(1.8)),
                       legend.text = element_text(size = rel(1.8)),
                       legend.title = element_text(size = rel(1.8)),
)

## Effort visualization
## plot effort x time faceted by program
dat.plot = dat %>% 
  group_by(year, source) %>%
  summarize(party.minutes = mean(party.minutes, na.rm=T),
            duration = mean(duration, na.rm=T))
  
ggplot(dat.plot %>% filter(!is.na(duration)), aes(x = year))+
  geom_point(aes(y = duration))+
  geom_path(aes(y = duration))+
  facet_wrap(. ~ source, scales = "free")+
  ylim(c(0,NA))+
  ggtitle("effort as \'duration\'")+
  theme.larger+
  theme(strip.text = element_text(size=rel(1.3)))

ggplot(dat.plot %>% filter(!is.na(party.minutes)), aes(x = year))+
  geom_point(aes(y = party.minutes))+
  geom_path(aes(y = party.minutes))+
  facet_wrap(. ~ source, scales = "free")+
  ylim(c(0,NA))+
  ggtitle("effort as \'party minutes\'")+
  theme_bw()+
  theme.larger

dat.plot = dat %>% 
  group_by(year, source) %>% 
  summarize(data.count = length(unique(event.id)))

ggplot(dat.plot, aes(x = year, y = data.count, group = source, color = source))+
  geom_path()+
  theme(legend.position = "none")+
  theme.larger

  facet_wrap(. ~ source)+
