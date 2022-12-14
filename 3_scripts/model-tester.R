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


## function to make site id for random effects.
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
  site.re = as.factor(site.re)
  site.re = relevel(site.re, ref = "light on data")
  return(as.factor(site.re))
}

## Function for specifying knots
## Knot placement:
## IF we're including inferred zeroes and not constraining geography, we may want
## to modify the default knot placement. mgcv:gam by default uses quantiles of the
## data to place knots (ie if data is more dense we want higher density of knots there).
## The inclusion of inferred zeroes can screw that up (tons of apparent data-density where
## nothing is actually happening).
## 
## We can specify knots based on the quantiles of the ACTUAL data 
## For cyclic splines (probably the best choice for day of year), we need to specify
## the endpoints that loop around to each other. I use day = 0.5 and day = 364.5,
## such that Dec 31 and Jan 1 are effectively 1 day apart (except for leap years)
##
## 
## Another option is to NOT specific knots. I've found specifying them works a bit
## better in my test cases, but feel free to experiment.
## 
## To make this easy to include in the loop, we define knot_maker as a function
## that takes in `dat` and returns a list. If the list is empty, knots will be generated
## automatically. Any item in the list that corresponds to a named variable in the model formula
## will rely on the provided knots

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

## Grab semi-USGS regions, designated by state based on this map: https://www.fws.gov/about/regions
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))

## specifying run name (to help identify/distinguish results files for different parameterizations)
run.suffix = "markdown-report-test" ## Change this for whatever you're trying out
## specify whether or not to use inferred 0s.
use.inferred = TRUE

## Specifying taxa:
specs.do.all = FALSE #If TRUE, try to apply model to all species with 400+data points
# If false, use the default four species defined here:

## should we constrain model to only region (convex hull) of non-inferred points?
geography.constrain = FALSE

## formula:
form.use = formula(count ~ te(lat, lon, by = year, k = c(10, 10), bs = c("cr", "cr")) +
                     sourcefac + 
                     s(site.refac, bs = 're')) ## can specify form listed above or use formula() to write it directly here.

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
specs.do = data.frame(code = c("PIERAP"),
                      specname = c("cabbage white")
)

# alternately, apply approach to ALL species
if(specs.do.all){
  dat = qread(paste0("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv"))
  code.tab = table(dat$code)
  code.tab = code.tab[code.tab >=400]
  code.tab = code.tab[names(code.tab != "?-SP")]
  code.tab = code.tab[names(code.tab != "TBD")]
  code.tab = code.tab[names(code.tab != "BFLY")]
  code.tab = code.tab[!grepl("-", names(code.tab))|
                        grepl("^S-", names(code.tab))]
  sum(dat$code %in% names(code.tab))
  specs.do = data.frame(code = names(code.tab))
  specs.do$specname = specs.do$code
}

## actual looping --------
summary.df = NULL
for(i.spec in 1:nrow(specs.do)){ 
  code.cur = specs.do$code[i.spec]
  out = model_runner(code.cur = code.cur,
                     form.use = form.use,
                     knot_maker = knot_maker,
                     sitere_maker = sitere_maker,
                     regions.dict = regions.dict,
                     use.inferred = use.inferred,
                     geography.constrain = geography.constrain,
                     use.only.source = use.only.source,
                     n.threads.use = n.threads.use)
  output.name = report_maker(out,
                            code.cur = code.cur,
                            run.suffix = run.suffix)
  summary.df = rbind(summary.df,
                     data.frame(
                       code = code.cur, 
                       abund.correlation = out$abund.cor,
                       filename = output.name
                     ))
}

## saving summary file ----------------

## to avoid accidental overwriting: automated numbering
cur.files = list.files(here(paste0("4_res/fit-summaries/")))
cur.file.code = cur.files[grepl(paste0("AAA-run-summary - ", run.suffix), cur.files)]
cur.file.code = cur.file.code[grepl("[.]txt", cur.file.code)]
if(length(cur.file.code)>0){
  cur.file.code = gsub("[.]txt", "", cur.file.code)
  cur.nums = as.numeric(gsub(paste0("AAA-run-summary - ", run.suffix, "-V"), "", cur.file.code))
  use.num = max(cur.nums)+1
}else{
  use.num = 1
}
path.use = here(paste0("4_res/fit-summaries/","AAA-run-summary - ", run.suffix, "-V", use.num, ".txt"))

parm_paster = function(parm.name){
  paste0(parm.name,":\n", "  ", as.character((get(parm.name))),"\n")
}

cat(c(paste0("Fit on: ", Sys.time(), "\n"),
      "Parameters:\n==================",
      parm_paster("form.use"),
      parm_paster("use.inferred"),
      parm_paster("geography.constrain"),
      parm_paster("use.only.source"),
      "\n"),
    sep = "\n",
    file = path.use)

cat("Correlation summaries:\n==================\n",
    append = TRUE,
    file = path.use)

write_delim(summary.df[,c("code","abund.correlation")], 
            file = path.use, col_names = TRUE, append = T,
            delim = "\t")

cat(c("\n\n Report files generated:\n==================",
      summary.df$filename),
    sep = "\n",
    file = path.use,
    append = TRUE)

cat("\n\n Human notes on model fit:\n==================\n[Fill in as desired]\n\n\n",
    file = path.use,
    append = TRUE)
cat(c("============================================================",
      "============================================================",
      "============================================================", "\n\n"),
    file = path.use,
    sep = "\n",
    append = TRUE)
cat("Extra details:  knot and site-making functions, region dictionary\n==================\n",
    file = path.use,
    append = TRUE)


sink(path.use, append = T)
cat("\n\nknot_maker():\n")
print(knot_maker)
cat("\n\nsitere_maker():\n")
print(sitere_maker)
cat("\n\nRegion dictionary:\n")
sink()
write_delim(regions.dict, 
            file = path.use, col_names = TRUE, append = T, delim = "\t\t\t")




