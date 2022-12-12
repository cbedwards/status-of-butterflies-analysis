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


##################################
##  Looping over multiple taxa  ##
##################################
## The following code systematizes fitting the same model to a range of taxa and 
## saving maps, diagnostic plots, and model summary info into a pdf report in 4_res/fit-summaries 
## with hopefully-identifiable names. Because some of the figures are pretty time-consuming to load
## as PDFs, this code also saves jpg versions in 4_res/fit-summaries/jpgs-fastviews.
## 
## In theory, only material in "Parameters for looping" and possibly "Models I'm exploring" needs to be changed to specify
## new models, data, etc
## The one exception is "Specifying Knots", which is in the middle of the loop, as 
## adaptive knot placement needs to update for each new data set. If you want to tweak
## knot placement, you'll want to modify the code leading up to and including `knots.list = ...`

{ # adding bracket to simplify running everything below.
  
  ### Models I'm exploring: --------
  # ## baseline formula:
  form.4 = as.formula(count ~ te(doy, lat, lon, k = c(7,  5, 5),  bs = c("cc", "cr", "cr")) +
                        te(lat, lon, by = year, k = c(5, 5), bs = c("cr", "cr")) + sourcefac)
  ## with shrinkage for doy x space
  form.5 = as.formula(count ~ te(doy, lat, lon, k = c(7,  5, 5),  bs = c("cc", "cs", "cs")) +
                        te(lat, lon, by = year, k = c(5, 5), bs = c("cr", "cr")) + sourcefac)
  ## activity curve static in each region
  form.6 = as.formula(count ~ s(doy,by = regionfac, k = c(7),  bs = c("cc")) +
                        te(lat, lon, by = year, k = c(5, 5), bs = c("cr", "cr")) + sourcefac)
  form.6.2 = as.formula(count ~ s(doy, by = regionfac, k = c(7),  bs = c("cc")) +
                          te(lat, lon, by = year, k = c(10, 10), bs = c("cr", "cr")) + sourcefac)
  ## Trying to allow doy to vary across years and regions
  form.6.2.1 = as.formula(count ~ te(doy, year, by = regionfac, k = c(5, 3),  bs = c("cc", "cr")) +
                            te(lat, lon, by = year, k = c(10, 10), bs = c("cr", "cr")) + sourcefac)
  form.6.2.2 = as.formula(count ~ te(doy, by = regionfac, k = c(7),  bs = c("cc")) +
                            te(doy, by = year, k = c(7),  bs = c("cc")) +
                            te(lat, lon, by = year, k = c(10, 10), bs = c("cr", "cr")) + sourcefac)
  
  form.6.3 = formula(count ~ te(doy, regionfac, k = c(6),  bs = c("cc", "re")) +
                       te(doy, by = year, k = c(6),  bs = c("cc")) +
                       te(lat, lon, by = year, k = c(10, 10), bs = c("cr", "cr")) + 
                       sourcefac + s(regionfac, bs = "re"))
  
  ## Gaussian activity curve static in each region
  form.7 = formula(count ~ regionfac + doy:regionfac + I(doy^2):regionfac +
                     te(lat, lon, by = year, k = c(10, 10), bs = c("cr", "cr")) + sourcefac)
  
  
  ## Parameters for looping ----------
  
  ## specifying run name (to help identify/distinguish results files for different parameterizations)
  run.suffix = "quadraticdoy" ## Change this for whatever you're trying out
  ## specify whether or not to use inferred 0s.
  use.inferred = TRUE
  
  ## Specifying taxa:
  specs.do.all = FALSE #If TRUE, try to apply model to all species with 400+data points
  # If false, use the default four species defined here:
  
  ## should we constrain model to only region (convex hull) of non-inferred points?
  geography.constrain = FALSE 
  
  ## formula:
  form.use = form.7 ## can specify form listed above or use formula() to write it directly here.
  
  
  
  ## How many threads to use when fitting? reduce if your computer is struggling.
  n.threads.use = 4
  
  ## actual looping --------
  
  #Set of four example species:
  specs.do = data.frame(code = c("PIERAP", 
                                 "NYMVAU",
                                 "VANCAR",
                                 "EUPHPHA"),
                        specname = c("cabbage white",
                                     "Compton tortoiseshell",
                                     "painted lady",
                                     "Baltimore checkerspot")
  )
  
  # alternately, apply approach to ALL species
  if(specs.do.all){
    dat = qread(paste0("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv"))
    code.tab = table(dat$code)
    code.tab = code.tab[code.tab >=400]
    code.tab = code.tab[names(code.tab != "?-SP")]
    code.tab = code.tab[names(code.tab != "TBD")]
    code.tab = code.tab[!grepl("-", names(code.tab))|
                          grepl("^S-", names(code.tab))]
    sum(dat$code %in% names(code.tab))
    specs.do = data.frame(code = names(code.tab))
    specs.do$specname = specs.do$code
  }
  
  
  for(i.spec in 1:nrow(specs.do)){ 
    code.cur = specs.do$code[i.spec]
    # code.cur = "DANPLE"
    specname.cur = specs.do$specname[i.spec]
    # specname.cur = "monarch butterfly"
    
    if(use.inferred){
      make_dataset(code.cur, specname.cur)
      dat = qread(paste0("2_data_wrangling/cleaned by code/",specname.cur, ".csv"))
    }else{
      ## curious if the inferred 0s are messing things up
      dat = qread(paste0("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv"))
      dat = dat %>%
        filter(code == code.cur)
    }
    dat = dat %>% 
      filter(!is.na(count))
    dim(dat)
    head(dat)
    
    dat$sourcefac = as.factor(dat$source)
    
   
    
    ## Grab semi-USGS regions, designated by state based on this map: https://www.fws.gov/about/regions
    regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))
    
    dat = left_join(dat, regions.dict, by = "state")
    ## FOR TESTING PURPOSES let's remove region = NA
    dat = dat %>% 
      filter(!is.na(dat$region))
    
    dat$regionfac = as.factor(dat$region)
    
    
    gp.counts = ggplot(dat) + 
      geom_point(aes(x = doy, y = count, col = inferred), shape = 1)+
      facet_wrap(~regionfac)+ggtitle("Counts across regions")+
      xlab("day of year")+
      theme_minimal()
    gp.hist = ggplot(dat) +
      geom_histogram(aes(x=doy, fill = count>0)) +
      facet_wrap(~ inferred)+
      ggtitle("Histogram of # records across doy\nseparated by whether inferred or not (TRUE = inferred)")+
      xlab("day of year")+
      theme_minimal()
    
    ## Get species names
    names.vec = unique(dat$name)
    ##remove "absence for" entries in list of species names -- these are to denote inferences
    names.vec = names.vec[!grepl("absence for ", names.vec)]
    plot.title = paste0(code.cur, ": ", paste0(names.vec, collapse = ", "))
    
    ## constraining geography if called for
    ## Thinking here is that splines can misbehave if we have many observations (of 0)
    ## in biologically uninteresting areas (ie lots of 0s on the east coast for a 
    ## west-coast-only species). This can happen because of our conferred zeroes
    ## 
    ## Constraining data might be best with range maps, but a quick-and-dirty solution
    ## is to only use observations within the convex hull of our non-inferred data.
    ## 
    ## Potential elaboration: add a buffer of ~X lat/lon around each non-inferred data point.
    ## 
    if(geography.constrain == TRUE){
      shape.hull = convhulln(dat[dat$inferred==FALSE,c("lon", "lat")])
      dat = dat[inhulln(shape.hull, as.matrix(dat[, c("lon","lat")])),]
    }
    
    ## Specifying knots -- THIS MIGHT NEED TO BE UPDATED
    ## Knot placement:
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
    # doy.knots = c(.5, 
    #               as.numeric(quantile(dat$doy[dat$inferred==FALSE],
    #                                   probs = seq(.05,.95, length = 4))),
    #               365.5)
    # lat.knots = as.numeric(quantile(dat$lat[dat$inferred==FALSE],
    #                                 probs = seq(0,1, length = 10)))
    # lon.knots = as.numeric(quantile(dat$lon[dat$inferred==FALSE],
    #                                 probs = seq(0,1, length = 10)))
    ## actual list of knots that will be used:
    # knots.list = list(doy = doy.knots,
    #                   lat = lat.knots,
    #                   lon = lon.knots)
    knots.list = list(doy = doy.knots) #if allowing lat/lon knots to be placed automatically
    # knots.list = list(doy = c(0.5, 365.5)) #default for not specifying knots - forces smooths to connect
    #             Dec 31 and Jan 1
    #knots.list = list() #if allowing ALL knots to be bplaced automatically (not recommended with cylic spline)    #
    
    print(paste0(specname.cur, " (", code.cur, "), Inferring zeros = ", use.inferred))
    time.start = proc.time()
    fit = bam(form.use,
              data = dat,
              method="fREML", 
              family="nb",
              knots=knots.list,
              discrete = TRUE,
              nthreads = n.threads.use)
    time.end = proc.time()
    print(time.end-time.start)
    print(fit$formula)
    
    if(FALSE){
      gam.check(fit)
    }
    
    ## generate plots --------------
    ## plot abundance map
    out.abund = abund_mapper(dat, fit, regions.dict, dat.constrain = geography.constrain)
    out.abund$fig
    
    ## plot trends
    out.trend = trend_plotter(dat, fit, regions.dict, 
                              dat.constrain = geography.constrain)
    out.trend$fig #+ scale_fill_viridis()
    
    ## plot activity curves for point of max data
    gp.activity.maxdata = demo_activity_plots(dat, fit, regions.dict)
    ## plot activity curves for point of max estimated density (diangostic for unreasonable activity curves)
    pt.maxabund = out.abund$data[which.max(out.abund$data$abund.index),]
    gp.activity.maxabund = activity_plotter(dat, fit, regions.dict, 
                                            lat.plot = pt.maxabund$lat, 
                                            lon.plot = pt.maxabund$lon,
                                            allyears = TRUE,
                                            source.adaptive = FALSE)
    
    ## Plot NFJ abundance by region
    gp.nfj.abund = NFJ_compare(dat, fit, regions.dict)
    gp.nfj.abund
    
    ## Plot NFJ trends by region
    gp.nfj.trend = NFJ_regional_trends(dat, regions.dict)
    gp.nfj.trend
    
    ### saving ---------
    ## identify next available version number (to avoid overwriting)
    cur.files = list.files(here(paste0("4_res/fit-summaries/")))
    cur.file.code = cur.files[grepl(code.cur, cur.files)]
    cur.file.code = cur.file.code[grepl("[.]pdf", cur.file.code)]
    if(length(cur.file.code)>0){
      cur.file.code = gsub("[.]pdf", "", cur.file.code)
      cur.nums = as.numeric(gsub(paste0(code.cur, "-", run.suffix,"-V"), "", cur.file.code))
      use.num = max(cur.nums)+1
    }else{
      use.num = 1
    }
    
    {pdf(here(paste0("4_res/fit-summaries/",
                     code.cur,
                     "-", run.suffix,
                     "-V", use.num, ".pdf")),
         width = 15, height = 20)
      ## map of abundance
      gp = ggarrange(out.abund$fig, gp.nfj.abund, ncol = 1)
      print(annotate_figure(gp, top = text_grob(plot.title, size = 18)))
      ## map of trends
      gp = ggarrange(out.trend$fig, gp.nfj.trend, ncol = 1)
      print(annotate_figure(gp, top = text_grob(plot.title, size = 18)))
      
      ## Turn knots.list into something human-readable. 
      knots.vec = "knots: "
      if(length(knots.list)==0){
        knots.vec = "knots: all automatically placed"
      }else{
        for(cur.name in names(knots.list)){
          knots.vec = paste0(knots.vec,
                             paste0(cur.name,": ", paste0(round(knots.list[[cur.name]],1), collapse = ", "), "\n       "))
        }
      }
      
      ## Page with summary information (in words)
      plot.words = plot.title
      plot.form = as.character(formula(fit))[3]
      plot.form = gsub("[+]", "+\n  ", plot.form)
      plot(NA, xlim = c(0,5), ylim = c(0,5), bty = "n", 
           xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      text(1,4, "Model fitting summary:", pos = 4, cex = 3)
      text(1,3, plot.words, pos = 4, cex = 2)
      text(1,2, plot.form, pos = 4, cex = 2)
      text(1,1, knots.vec, pos = 4)
      
      gp = ggarrange(gp.hist, gp.counts, ncol=1)
      print(annotate_figure(gp, top = text_grob(plot.title, size = 18)))
      
      ## activity curves for highest data denstiy, highest estimated abundance
      print(gp.activity.maxdata + ggtitle(
        paste0(gp.activity.maxdata$labels$title, "\n(Location of highest data density)")
      ))
      print(gp.activity.maxabund + ggtitle(
        paste0(gp.activity.maxabund$labels$title, "\n(Location of highest estimated abundance)")
      ))
      dev.off()
      ## save jpgs for faster viewing of maps etc.
      ggsave(here(paste0("4_res/fit-summaries/jpgs-fastviews/",
                         code.cur, "-", run.suffix,
                         "-V", use.num, "-abund-map.jpg")),
             out.abund$fig, 
             width = 15, height = 12)
      ggsave(here(paste0("4_res/fit-summaries/jpgs-fastviews/",
                         code.cur, "-", run.suffix,
                         "-V", use.num, "-trend-map.jpg")),
             out.trend$fig, 
             width = 15, height = 12)
      ggsave(here(paste0("4_res/fit-summaries/jpgs-fastviews/",
                         code.cur, "-", run.suffix,
                         "-V", use.num, "-trend-NFJ.jpg")),
             gp.nfj.trend, 
             width = 15, height = 12)
    }
  }
}
