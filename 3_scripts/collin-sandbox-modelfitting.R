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
form.use = formula(count ~ te(lat, lon, by = year, k = c(5, 5), bs = c("cr", "cr")) +
                     effort.universal+#:effort.universal.type +
                     #effort.universal +
                     s(doy, k = 6) +
                     # sourcefac +
                     s(site.refac, bs = 're')) ## can specify form listed above or use formula() to write it directly here.
## we also need to give it a region dictionary

## Grab FWS regions, designated by state based on this map: https://www.fws.gov/about/regions
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))

##Note: regions for FWS

out = model_runner(code.cur = "ATACAM",
                   form.use = form.use, #formula for the model
                   knot_maker = knot_maker, #knot-generating function. We can leave it be
                   sitere_maker = sitere_maker, #function for defining custom "sites"
                   regions.dict = regions.dict, #dictionary for identifying regions. 
                   fit.family = "nb", #error distribution for model fit. nb = negative binomial
                   use.range = TRUE, #constrain analysis to the range maps provided by Eliza. TRUE is good here.
                   use.inferred = TRUE, #infer zeroes from community-gathering trips? TRUE is good here
                   infer.messy.levels = c("GENUS", "SUBFAMILY", "FAMILY", "COMPLEX"), #When inferring zeroes, how do we deal with "Unknown in genus of interest"? 
                   ## specified levels of unknowns are ignored when inferring zeroes.
                   geography.constrain = FALSE, #superceeded by the range maps. If TRUE, restricts analysis to convex hull of non-zero observations
                   use.only.source = NULL, #Can specify using data only from individual sources. NULL means use all sources
                   n.threads.use = 2, # how many cores to us?
                   # pheno.window = c(0.001, 0.999),
                   regions.use = "Northeast",
                   min.year = -9999,
                   do.pheno = TRUE) #if including phenology in your model, set to TRUE to calculate abundance separately for each day of year.

out$fig.abund
out$fig.trend


form.use = formula(count ~ te(lat, lon, by = year, k = c(1, 1), bs = c("cr", "cr")) +
                     effort.universal+#:effort.universal.type +
                     #effort.universal +
                     s(doy, k = 6) +
                     # sourcefac +
                     s(site.refac, bs = 're')) ## can specify form listed above or use formula() to write it directly here.
## we also need to give it a region dictionary

## Grab FWS regions, designated by state based on this map: https://www.fws.gov/about/regions
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))

##Note: regions for FWS

out3 = model_runner(code.cur = "ATACAM",
                   form.use = form.use, #formula for the model
                   knot_maker = knot_maker, #knot-generating function. We can leave it be
                   sitere_maker = sitere_maker, #function for defining custom "sites"
                   regions.dict = regions.dict, #dictionary for identifying regions. 
                   fit.family = "nb", #error distribution for model fit. nb = negative binomial
                   use.range = TRUE, #constrain analysis to the range maps provided by Eliza. TRUE is good here.
                   use.inferred = TRUE, #infer zeroes from community-gathering trips? TRUE is good here
                   infer.messy.levels = c("GENUS", "SUBFAMILY", "FAMILY", "COMPLEX"), #When inferring zeroes, how do we deal with "Unknown in genus of interest"? 
                   ## specified levels of unknowns are ignored when inferring zeroes.
                   geography.constrain = FALSE, #superceeded by the range maps. If TRUE, restricts analysis to convex hull of non-zero observations
                   use.only.source = NULL, #Can specify using data only from individual sources. NULL means use all sources
                   n.threads.use = 2, # how many cores to us?
                   # pheno.window = c(0.001, 0.999),
                   regions.use = "Northeast",
                   min.year = -9999,
                   do.pheno = TRUE) #if including phenology in your model, set to TRUE to calculate abundance separately for each day of year.

out3$fig.abund
out3$fig.trend

plot_grid(out$fig.trend, out3$fig.trend)

# temp = predict(out$fitted.model, newdata = out$data, type = 'response')
# summary(temp)


sitere_maker = function(dat){
  ##identify the sites NOT to lump
  temp = dat %>% 
    group_by(site, year) %>%
    summarize(mean.cont = mean(count)) %>% 
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

# sitere_maker(dat.spec)


## formula for the model fit
form.use = formula(count ~ te(lat, lon, by = year, k = c(10, 10), bs = c("cr", "cr")) +
                     effort.universal:effort.universal.type +
                     effort.universal +
                     s(doy, k = 6) +
                     # sourcefac +
                     s(site.refac, bs = 're')) ## can specify form listed above or use formula() to write it directly here.
## we also need to give it a region dictionary

## Grab FWS regions, designated by state based on this map: https://www.fws.gov/about/regions
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))

##Note: regions for FWS

out2 = model_runner(code.cur = "ATACAM",
                    form.use = form.use, #formula for the model
                    knot_maker = knot_maker, #knot-generating function. We can leave it be
                    sitere_maker = sitere_maker, #function for defining custom "sites"
                    regions.dict = regions.dict, #dictionary for identifying regions. 
                    fit.family = "nb", #error distribution for model fit. nb = negative binomial
                    use.range = TRUE, #constrain analysis to the range maps provided by Eliza. TRUE is good here.
                    use.inferred = TRUE, #infer zeroes from community-gathering trips? TRUE is good here
                    infer.messy.levels = c("GENUS", "SUBFAMILY", "FAMILY", "COMPLEX"), #When inferring zeroes, how do we deal with "Unknown in genus of interest"? 
                    ## specified levels of unknowns are ignored when inferring zeroes.
                    geography.constrain = FALSE, #superceeded by the range maps. If TRUE, restricts analysis to convex hull of non-zero observations
                    use.only.source = NULL, #Can specify using data only from individual sources. NULL means use all sources
                    n.threads.use = 2, # how many cores to us?
                    # pheno.window = c(0.001, 0.999),
                    regions.use = NULL,
                    min.year = -9999,
                    do.pheno = TRUE) #if including phenology in your model, set to TRUE to calculate abundance separately for each day of year.
out2$fig.abund
out2$fig.trend

library(cowplot)
plot_grid(out$fig.abund, out2$fig.abund)
plot_grid(out$fig.trend, out2$fig.trend)
