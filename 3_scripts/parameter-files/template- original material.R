## parameter file to be sourced in 3_scripts/model-runner-parameterized.R
## 
## 
## ## function to make site id for random effects.
## treat any site with <5 years of non-zero data as a generic site
sitere_maker = function(dat){
  temp = dat %>% 
    group_by(site, year) %>% 
    summarize(has.nonzero = any(count>0)) %>% 
    filter(has.nonzero == TRUE) %>% 
    group_by(site) %>% 
    summarize(nyear.nonzero = n()) %>% 
    filter(nyear.nonzero < 5)
  site.re = rep("light on data", nrow(dat))
  site.re[dat$site %in% temp$site] = dat$site[dat$site %in% temp$site]
  ## adding the following line to reduce model complexity. With this line, random effects
  ## are basically just there to capture effort.
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

## default knot placement, but forcing doy to wrap from day 365.5 to day 0.5 (good baseline when including doy as a cyclic spline)
# knot_maker = function(dat){
#   knot.list = list(doy = c(0.5, 365.5))
# } 

# adaptive knot placement for doy - in addition to wrapping from day 365.5 to day 0.5, remaining knots are placed
# based on quantiles of NON-INFERRED data. Length of the `probs` sequence should be k-2
# knot_maker = function(dat){
#   doy.knots = c(.5,
#                 as.numeric(quantile(dat$doy[dat$inferred==FALSE],
#                                     probs = seq(.2,.8, length = 4))),
#                 365.5)
#   return(list(doy = doy.knots))
# }

# adaptive knot placement for doy (as above) AND lat/lon, again based on the quantiles of the
# non-inferred data. note that the length for the probs sequence in lat and lon should equal k

# knot_maker = function(dat){
#   doy.knots = c(.5,
#                 as.numeric(quantile(dat$doy[dat$inferred==FALSE],
#                                     probs = seq(.2,.8, length = 4))),
#                 365.5)
#   lat.knots = as.numeric(quantile(dat$lat[dat$inferred==FALSE],
#                                   probs = seq(0,1, length = 10)))
#   lon.knots = as.numeric(quantile(dat$lon[dat$inferred==FALSE],
#                                   probs = seq(0,1, length = 10)))
#   return(list(doy = doy.knots,
#               lat = lat.knots,
#               lon = lon.knots))
# }


## Parameters for looping ----------

do.summary = TRUE #If true, create summary text file for the whole combination of runs. 
## expect to use this once we settle on a fit.

## Grab semi-USGS regions, designated by state based on this map: https://www.fws.gov/about/regions
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))

## specifying run name (to help identify/distinguish results files for different parameterizations)
run.suffix = "full-run-allsources-redo" ## Change this for whatever you're trying out
## specify whether or not to use inferred 0s.
use.inferred = TRUE

## should we constrain data to within the range maps? Seems like a good idea when we have them
use.range = TRUE

## Specifying taxa:
specs.do.all = TRUE #If TRUE, try to apply model to all species with 400+data points
# If false, use the default four species defined here:

## should we constrain model to only region (convex hull) of non-inferred points?
geography.constrain = FALSE

## formula:
## "full" model
form.use = formula(count ~ te(lat, lon, by = year, k = c(5, 5), bs = c("cr", "cr")) +
                     effort.universal:effort.universal.type +
                     sourcefac +
                     s(site.refac, bs = 're')) ## can specify form listed above or use formula() to write it directly here.
do.pheno = FALSE
## NFJ-only model
# form.use = formula(count ~ te(lat, lon, by = year, k = c(5, 5), bs = c("cr", "cr")) +
#                     effort.universal +
#                    s(site.refac, bs = 're')) ## can specify form listed above or use formula() to write it directly here.
do.pheno = FALSE

## NOTE: if NOT fitting a doy term in the model, you can set `do.pheno` to FALSE, 
## greatly speading up estimation of abundance + trends from the fitted model.
## `do.pheno` MUST BE SET TO TRUE if there is a "doy" term in the model, so that 
## abundance and trend are estimated from the full activity curve, not from 2 day per year.

## select a subset of the data sources to use
##   if left as NULL, will use all sources
use.only.source = NULL

## How many threads to use when fitting? reduce if your computer is struggling.
n.threads.use = 4


## choosing species

#Set of four example species:
#note that the only important column is $code. I include specname for ease of reading.
specs.do = data.frame(code = c("PIERAP", 
                               "NYMVAU",
                               "VANCAR",
                               "EUPHPHA"),
                      specname = c("cabbage white",
                                   "Compton tortoiseshell",
                                   "painted lady",
                                   "Baltimore checkerspot")
)
# specs.do = data.frame(code = c("EUPHPHA"),
#                       specname = c("Baltimore checkerspot")
# )

# alternately, apply approach to ALL species
# Currently: requiring we have a map.
if(specs.do.all){
  dat = qread(paste0("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv"))
  #count the number of years, sites, and observation events with > 0 butterflies reported (per species)
  dat.sum = dat %>% 
    filter(presence == 1) %>% 
    group_by(code) %>% 
    summarise(nyear = length(unique(year)),
              nsite = length(unique(site)),
              nsource = length(unique(source)),
              nsource.rich = sum(table(source)>10),
              neffort.type = length(unique(effort.universal.type)),
              neffort.type.rich =  sum(table(effort.universal.type)>10),
              nobs = n()
    ) %>% 
    ungroup()
  dat.use = dat.sum %>% ## identifying thresholds for minimum years of data, minimum sites, minimum obs 
    filter(nobs >= 100) %>% 
    filter(nyear >= 10 ) %>% 
    filter(nsite >= 10 ) %>% 
    filter(nsource > 1 ) %>% 
    filter(neffort.type > 1)
  dat.use = dat.use[!(grepl("-", dat.use$code)),]
  dat.use = dat.use %>% 
    filter(code != "TBD") %>% 
    filter(code != "BFLY") %>% 
    filter(code != "NONE")
  codes.use = dat.use$code
  maps.list = list.files(here("2_data_wrangling/range-maps/relabeled/"))
  maps.list = grep("*.kml", maps.list, value = TRUE)
  maps.list = gsub(".kml", "", maps.list)
  codes.use = codes.use[codes.use %in% maps.list]
  specs.do = data.frame(code = codes.use)
  specs.do$specname = codes.use
}
problem.codes = c("CALLDUM", "CALLPOL", "CHLOPAL", "ERYZAR", "EUPYBIM")
## quick and dirty solution: skip those five species.
specs.do = specs.do[-((which(specs.do$code %in% problem.codes))),]