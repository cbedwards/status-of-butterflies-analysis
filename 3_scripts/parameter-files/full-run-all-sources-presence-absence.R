## parameter file to be sourced in 3_scripts/model-runner-parameterized.R
## 

## Functions ------------------

## These two functions allow tweaking of site identity for use in random effects, and
## for custom knot placement. Custom knot placement seemed more important with inferred
## zeroes and no geographic constraints. Can likely leave alone now. 
## The sitere_maker() allows us to avoid estimating separate random effects for
## sites without sufficient data. See the final "filter()" call in designating temp.


## ## function to make site id for random effects.
## treat any site with <3 observations as a generic site
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


## Parameters ----------
## Functional parameters
do.summary = TRUE #If true, create summary text file for the whole combination of runs. 
use.inferred = TRUE # specify whether or not to use inferred 0s.
fit.family = "binomial" #negative binomial fit
use.range = TRUE ## should we constrain data to within the range maps?
geography.constrain = FALSE ## should we constrain model to convex hull of non-zero counts? Superceded by use.range
use.only.source = NULL ## select a subset of the data sources to use. If left as NULL, will use all sources

## performance-tweaking parameters
copy.heatmaps = TRUE ## If true, save a copy of the abundance and trends heatmaps for each species. Takes extra time to do this.
n.threads.use = 4 ## How many threads to use when fitting? reduce if your computer is struggling.
do.pheno = FALSE ## !!IF!! the formula does not include phenology in the form of a doy-related term,
## setting this to FALSE will substantially speed up runs. WARNING: If do.pheno is FALSE and formula includes doy, results will be WRONG.

## formula: ------------
form.use = formula(presence ~ te(lat, lon, by = year, k = c(5, 5), bs = c("cr", "cr")) +
                     effort.universal:effort.universal.type + 
                     sourcefac +
                     s(site.refac, bs = 're')) ## can specify form listed above or use formula() to write it directly here.



## choosing species ---------------------
#Can manually create specs.do, as in the commented example here, and set specs.do.all = FALSE. 
#Alternately, can leave specs.do.all = TRUE, and fit each species that fits the given criterion.
#note that the only important column is $code. I include specname for ease of reading.

specs.do.all = TRUE #If TRUE, try to apply to most taxa (see code below)
# specs.do = data.frame(code = c("PIERAP", 
#                                "NYMVAU",
#                                "VANCAR",
#                                "EUPHPHA"),
#                       specname = c("cabbage white",
#                                    "Compton tortoiseshell",
#                                    "painted lady",
#                                    "Baltimore checkerspot")
# )

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
  cat("The following codes have no range map, and will be excluded:\n")
  print(codes.use[!(codes.use %in% maps.list)], collapse = ", ")
  codes.use = codes.use[codes.use %in% maps.list]
  specs.do = data.frame(code = codes.use)
  specs.do$specname = codes.use
}

problem.codes = c("CALLNIP", "CHLOGOR", "ERYZAR", "PIRPIR")

# Debugging: run everything after that
# specs.do = specs.do[-(1:max(which(specs.do$code %in% problem.codes))),]


## quick and dirty solution: skip those species.
## When doing final run
specs.do = specs.do[-((which(specs.do$code %in% problem.codes))),]

### Misc ------------
## Grab FWS regions, designated by state based on this map: https://www.fws.gov/about/regions
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))
