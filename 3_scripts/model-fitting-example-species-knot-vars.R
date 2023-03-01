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
source(here("3_scripts/funs.R"))
source(here("3_scripts/funs-fitting.R"))
source(here("3_scripts/model-fitting-function.R"))

### Example: cabbage white butterfly -------
##Variables to adjust species and knots 
species.id<-"ATACAM" #species code for analysis
k.lat<-40 #static # knots across latitude
k.lon<-40 #static # knots accross longitude
k.doy<-4  #static # knots across doy if pheno smoothing
pheno<-TRUE #True or False for  smoothing across doy (check smooth structure on line 59)
savefigs<-FALSE #True or False for saving ggplot and summary output files; see line 98
## First, the model_runner() function can do most of the steps, and spit out a bunch of
## useful results + figures.

## Let's define the formula, the random effect generating function, and 
## the knot-making function.

sitere_maker = function(dat){
  ##identify the sites NOT to lump
  temp = dat %>% 
    group_by(site) %>% 
    summarize(n.obs = n()) %>% 
    filter(n.obs > 2)
  cat(paste0("combining sites with low data into one random effect, total of ", nrow(temp)," observations affected \n"))
  site.re = rep("light on data", nrow(dat))
  site.re[dat$site %in% temp$site] = dat$site[dat$site %in% temp$site]
  ## adding the following line to reduce model complexity. With this line, random effects
  ## are basically just there to capture effort
  # site.re[dat$effort.universal.type!="site-based"] = "effort reported"
  site.re = as.factor(site.re)
  if(sum(site.re == "light on data")>0){
    site.re = relevel(site.re, ref = "light on data")
  }
  return(as.factor(site.re))
}

## default knot placement:
knot_maker = function(dat){
  return(list())
}

## formula for the model fit
if(pheno) {
  form.use = formula(count ~ te(lat, lon, by = year, k = c(k.lat, k.lon), bs = c("cr", "cr")) +
                       effort.universal:effort.universal.type + 
                       sourcefac +
                       s(doy,  k=k.doy) + #by=as.factor(region)
                       s(site.refac, bs = 're')) ## can specify form listed above or use formula() to write it directly here.
  ## we also need to give it a region dictionary
  
} else {
form.use = formula(count ~ te(lat, lon, by = year, k = c(k.lat, k.lon), bs = c("cr", "cr")) +
                     effort.universal:effort.universal.type + 
                     sourcefac +
                     #s(doy, by=as.factor(region), k=k.doy) +
                     s(site.refac, bs = 're')) ## can specify form listed above or use formula() to write it directly here.
## we also need to give it a region dictionary

}

## Grab FWS regions, designated by state based on this map: https://www.fws.gov/about/regions
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))


out = model_runner(code.cur = species.id,
                   form.use = form.use, #formula for the model
                   knot_maker = knot_maker, #knot-generating function. We can leave it be
                   sitere_maker = sitere_maker, #function for defining custom "sites"
                   regions.dict = regions.dict, #dictionary for identifying regions. 
                   fit.family = "nb", #error distribution for model fit. nb = negative binomial
                   min.year=2020,
                   use.range = TRUE, #constrain analysis to the range maps provided by Eliza. TRUE is good here.
                   use.inferred = TRUE, #infer zeroes from community-gathering trips? TRUE is good here
                   infer.messy.levels = c("GENUS", "SUBFAMILY", "FAMILY", "COMPLEX"), #When inferring zeroes, how do we deal with "Unknown in genus of interest"? 
                   ## specified levels of unknowns are ignored when inferring zeroes.
                   geography.constrain = FALSE, #superceeded by the range maps. If TRUE, restricts analysis to convex hull of non-zero observations
                   use.only.source = NULL, #Can specify using data only from individual sources. NULL means use all sources
                   n.threads.use = 2, # how many cores to us?
                   do.pheno = pheno) #if including phenology in your model, set to TRUE to calculate abundance separately for each day of year.


out$fig.abund
out$fig.trend
summary(out$loc.plot)

if(savefigs) {
if(pheno) {
ggsave(out$fig.abund, file=paste0(species.id,"-",k.lat,"x",k.lon,"-p",k.doy,"ab.png"))
ggsave(out$fig.trend, file=paste0(species.id,"-",k.lat,"x",k.lon,"-p",k.doy,"tr.png"))
write_csv(as.data.frame(summary(out$loc.plot)), file=paste0(species.id,"-p",k.lat,"x",k.lon,"-",k.doy,".csv"))
} else {
  ggsave(out$fig.abund, file=paste0(species.id,"-",k.lat,"x",k.lon,"-","ab.png"))
  ggsave(out$fig.trend, file=paste0(species.id,"-",k.lat,"x",k.lon,"-","tr.png"))
  write_csv(as.data.frame(summary(out$loc.plot)), file=paste0(species.id,"-",k.lat,"x",k.lon,"-",".csv"))
}
}