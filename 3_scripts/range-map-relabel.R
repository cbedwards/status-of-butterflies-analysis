# To relabel Eliza's range maps, from: https://github.com/elizagrames/butterflyblobs

## relabeling range maps -----------
library(here)
library(tidyverse)
source(here("3_scripts/funs.R"))

vec_rm = function(vec, entry.remove){
  vec[! vec %in% entry.remove]
}


## Identify species we're likely to work with:
dat = qread(paste0("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv"))
#count the number of years, sites, and observation events with > 0 butterflies reported (per species)
dat.sum = dat %>% 
  filter(presence == 1) %>% 
  group_by(code) %>% 
  summarise(nyear = length(unique(year)),
            nsite = length(unique(site)),
            nobs = n()
  ) %>% 
  ungroup()
dat.use = dat.sum %>% ## identifying thresholds for minimum years of data, minimum sites, minimum obs 
  filter(nobs >= 100) %>% 
  filter(nyear >= 5 ) %>% 
  filter(nsite >= 5 )
dat.use = dat.use[!(grepl("-", dat.use$code)),]
specs.do = dat.use$code
specs.do = vec_rm(specs.do, c("TBD", "BFLY", "NONE"))
specs.do = na.omit(specs.do)


## relabeling range maps
dict = read_csv(here("2_data_wrangling/dictionaries/dictionary-use.csv"))
dict = dict[dict$code %in% specs.do,]
dict$name = gsub(" ", "_", dict$name)
dict$name = tolower(dict$name)

file.names = list.files(here("2_data_wrangling/range-maps"))
file.names = grep("_outline.kml", file.names, value = T)
file.names = gsub("_outline.kml", "", file.names) %>% 
  tolower()

## identify codes that are missing maps. 
## NOTE: we have multiple species names for some codes.
## For each code, check if there is not exactly one map for a species name
## associated with that code.
matches = NULL
for(cur.code in unique(dict$code)){
  dict.cur = dict %>% filter(code == cur.code)
  res.cur =  data.frame(code = cur.code,
                        num.matches = sum(dict.cur$name %in% file.names),
                        files.matches = paste(dict.cur$name[dict.cur$name %in% file.names], collapse = " | "))
  matches = rbind(matches, res.cur)
  ## if possible, update 
  if(res.cur$num.matches == 1){
    map.name = paste0(file.names[file.names %in% dict.cur$name],"_outline.kml")
    new.name = paste0(cur.code, ".kml")
    file.copy(from = here("2_data_wrangling/range-maps/", map.name),
              to = here("2_data_wrangling/range-maps/relabeled/", new.name),)
  }
}

## Copy OCHSUL to OCHYUM, since they share a range map
file.copy(from = here("2_data_wrangling/range-maps/relabeled/OCHSYL"),
          to = here("2_data_wrangling/range-maps/relabeled/OCHYUM"))

# View(matches)
trouble = matches %>% filter(num.matches != 1)
trouble = trouble %>% 
  filter(code != "OCHYUM")

write_csv(trouble, here("2_data_wrangling/range-maps/range-map-mismatch.csv"))
