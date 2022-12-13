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
source(here("3_scripts/funs.R"))
source(here("3_scripts/funs-fitting.R"))

### Example: cabbage white butterfly -------
## Inferring zeroes from community counts, not constraining geography

## Make data with inferred 0s, read it in
make_dataset(code = "PIERAP", name.pretty = "PIERAP")
dat = qread(paste0("2_data_wrangling/cleaned by code/","PIERAP", ".csv"))
## remove any observations with NA counts
dat = dat %>% 
  filter(!is.na(count))
dim(dat)
head(dat)
## we need to turn source into a factor for gams
dat$sourcefac = as.factor(dat$source)


## Grab semi-USGS regions, designated by state based on this map: https://www.fws.gov/about/regions
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))
dat = left_join(dat, regions.dict, by = "state")
## Currently we have some observations that are at a coast or beyond it, and 
## our state identification is giving those a state of ""
# loc_viz(dat %>% filter(state == ""))
# We need to cut out observations without a state (and thus without a region)
dat = dat %>% 
  filter(!is.na(dat$region))

dat$regionfac = as.factor(dat$region)

## Activity by region and day of year, colored by whether the observation is inferred
## or actively reported
ggplot(dat) + 
  geom_point(aes(x = doy, y = count, col = inferred), shape = 1)+
  facet_wrap(~regionfac)+ggtitle("Counts across regions")+
  xlab("day of year")+
  theme_minimal()
## How many of our observations are actively observed vs inferred?
ggplot(dat) +
  geom_histogram(aes(x=doy, fill = count>0)) +
  facet_wrap(~ inferred)+
  ggtitle("Histogram of # records across doy\nseparated by whether inferred or not (TRUE = inferred)")+
  xlab("day of year")+
  theme_minimal()

## Model formula:

# Feel free to play around. Some features that I think/have found to be important:
#   0) For interactions between different terms, we want to default to te(), which fit tensor
#      product smooths. This handles interactions between variables that are not on the same scale
#      (the default interactions in s() assume wiggliness is the same for a 1-unit change
#      in any dimension). 
#   1) We have observations across the year, even when a given species isn't active.
#        Including a day of year ("doy") smooth allows the model to capture across-year variation. 
#        BUT the doy smooth can sometimes misbehave, and create "U" shapes that produce
#        wildly innacurate estimates of abundance. One solution that has worked in a number
#        of test cases is a doy x region interaction, with region as a random effect.
#        In general, I've found that using a cyclic spline (bs = "cc") for day of year
#        also helps minimize model misbehavior, since it forces the beginning and end of the 
#        year to match (in both value and derivative).
#   2) Using the `by = year` argument allows us to have an interaction between a smooth
#       and year as a continuous variable. Notably, we are also using this to semi-directly estimate
#       population growth rate, since by=year fits a log-linear relationship of count across years
#      (assuming we use a distribution with a log link like negative binomial, which we are doing)

form.use = formula(count ~ te(doy, regionfac, k = c(6),  bs = c("cc", "re")) +
                     # te(doy, by = year, k = c(6),  bs = c("cc")) +
                     te(lat, lon, by = year, k = c(10, 10), bs = c("cr", "cr")) 
                   + sourcefac + s(regionfac, bs = "re"))

##   I'm going to specify knot locations(below), but if I didn't we would want to 
##   specify k in the formulas, as here:
# form.use = formula(count ~ te(doy, regionfac, k = c(6),  bs = c("cc", "re")) +
#                      te(doy, by = year, k = c(6),  bs = c("cc")) +
#                      te(lat, lon, by = year, k = c(10, 10), bs = c("cr", "cr")) + sourcefac)

## Specifying knots
## IF we're including inferred zeroes and not constraining geography, we may want
## to modify the default knot placement. mgcv:gam by default uses quantiles of the
## data to place knots (ie if data is more dense we want higher density of knots there).
## The inclusion of inferred zeroes can screw that up (tons of apparent data-density where
## nothing is actually happening).
## We can specify knots based on the quantiles of the ACTUAL data 
## For cyclic splines (probably the best choice for day of year), we need to specify
## the endpoints that loop around to each other. I use day = 0.5 and day = 364.5,
## such that Dec 31 and Jan 1 are effectively 1 day apart (except for leap years)
##
## If we're doing this, the gam will use the list to determine the NUMBER of knots 
## as well 
## 
## Another option is to NOT specific knots. I've found specifying them works a bit
## better in my test cases, but feel free to experiment.
doy.knots = c(.5, 
              as.numeric(quantile(dat$doy[dat$inferred==FALSE],
                                  probs = seq(.05,.95, length = 4))),
              365.5)
lat.knots = as.numeric(quantile(dat$lat[dat$inferred==FALSE],
                                probs = seq(0,1, length = 10)))
lon.knots = as.numeric(quantile(dat$lon[dat$inferred==FALSE],
                                probs = seq(0,1, length = 10)))
## actual list of knots that will be used:
# knots.list = list(doy = doy.knots,
#                   lat = lat.knots,
#                   lon = lon.knots)
knots.list = list(doy = doy.knots)
# knots.list = list(doy = doy.knots,
#                   lat = lat.knots,
#                   lon = lon.knots)

## On my computer, takes ~30 sec for cabbage white. YMMV.
fit = bam(form.use, #bam is a variant of gam for very large data sets: lower memory 
          # use, sometimes faster, built-in multi-threading
          data = dat,
          method="fREML", #Fast REML calculations
          family="nb", #negative binomial
          knots=knots.list,# using our knot placements from above
          discrete = TRUE, #speeds up operations
          nthreads = 4) #speeds up operations by using multiple threads
# YOU MAY NEED TO REDUCE THIS for laptops
# gam.check(fit) #this does some model diagnostics for gams, but depending on the model
# can take a LONG time. 



## abund_mapper() calculates average abundance across all years on a spatial grid.
##  Code in funs-fitting.R. outputs list, with $fig being a fiture, and $data
##  giving calculated data.
out.abund = abund_mapper(dat, fit, regions.dict, dat.constrain = FALSE)
out.abund$fig

## trend_plotter() calculates population across a spatial grid. Note that it presumes
## we're fitting a loglinear relationship to year (so comparing two years gives you
## the exact estimate of growth rate for all years). Again, $fig is a figure, $data is
## the estimated trend through space
out.trend = trend_plotter(dat, fit, regions.dict, color.zoom = FALSE)
out.trend$fig #+ scale_fill_viridis()
## The figure uses a color gradient based on distance above/below 0. This is useful 
## when thinking about actual growth rates, less useful when looking for nuanced spatial 
## differences in trends. + scale_fill_viridis() changes the color scheme

## When including a smooth across DOY, it's important to check that it's not being
## fit as U-shaped. 

## plot activity curves for the lat/lon of max data density
demo_activity_plots(dat, fit, regions.dict)
## plot activity curves for the lat/lon of maximum estimated density
## (if activity is estimated to be U-shaped, estimated density tends to explode)
pt.maxabund = out.abund$data[which.max(out.abund$data$abund.index),]
activity_plotter(dat, fit, regions.dict, 
                 lat.plot = pt.maxabund$lat, 
                 lon.plot = pt.maxabund$lon,
                 allyears = TRUE,
                 source.adaptive = FALSE)

## Model evaluation: compare the average NFJ abundance across the last X years
## vs the average predicted abundance at those locations across the same years
NFJ_compare(dat, fit, regions.dict, nyears = 10)

## To compare with abundance map, plot NFJ trends by region:
NFJ_regional_trends(dat, regions.dict)

