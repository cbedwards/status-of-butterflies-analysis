## visualizing data density.
## 
library(here)
library(tidyverse)
library(sp)
library(sf)
library(spData)
library(rworldmap)
library(ggpubr)
library(scatterpie)
source(here("3_scripts/funs.R"))

## set up map

state.map.data <- maps::map('state', fill = TRUE, plot = FALSE) %>%
  st_as_sf()

dat = qread("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv")
dat = dat[!is.na(dat$count),]

## going to skip 

code.cur = "CELLAD"


## examples for playing with settings -----------

## NFJ demo


dat_filt = function(dat, code.cur, dat.source, rnd.digits = 0){
  ## filter data to the nearest lat/lon (or finer scale for rnd.digits>0)
  ## and summarizer number of years of coverage (nyear), number of observations (nobs)
  ##  [includes reported 0s but not inferred 0s], and number of butterflies reported total (ncount)
  ##  Note that for nyear, we're not counting years in which only 0s were reported.
  if(dat.source %in% c("pollard", "Pollard")){
    dat.cur = dat %>% 
      filter(code == code.cur) %>% 
      filter(source != "NFJ") %>% 
      filter(source != "MASSBfly") %>% 
      mutate(lon.rnd = round(lon, digits= rnd.digits),
             lat.rnd = round(lat, digits= rnd.digits))
  }else{
  dat.cur = dat %>% 
    filter(code == code.cur) %>% 
    filter(source == dat.source) %>% 
    mutate(lon.rnd = round(lon, digits= rnd.digits),
           lat.rnd = round(lat, digits= rnd.digits))
  }
  dat.year = dat.cur %>%
    filter(count > 0) %>% 
    select(year, lon.rnd, lat.rnd) %>% 
    unique() %>% 
    group_by(lon.rnd, lat.rnd) %>% 
    summarize(nyear = n()) %>% 
    ungroup()
  dat.obs = dat.cur %>%
    select(lon.rnd, lat.rnd, date) %>% 
    unique() %>% 
    group_by(lon.rnd, lat.rnd) %>% 
    summarize(nobs = n()) %>% 
    ungroup()
  dat.count = dat.cur %>% 
    group_by(lon.rnd, lat.rnd) %>% 
    summarize(ncount = sum(count)) %>% 
    ungroup()
  dat.sum = inner_join(dat.year, dat.obs, by = c("lon.rnd", "lat.rnd"))
  dat.sum = inner_join(dat.sum, dat.count, by = c("lon.rnd", "lat.rnd"))
  return(dat.sum)
}


dat_filt_pie = function(dat, code.cur, dat.source, rnd.digits = 0){
  dat %>% 
    filter(code == code.cur) %>%
    filter(source == source) %>% 
    mutate(lon.rnd = round(lon, digits= rnd.digits),
           lat.rnd = round(lat, digits= rnd.digits)) %>% 
    group_by(name, lon.rnd, lat.rnd) %>% 
    summarise(value = sum(count)) %>% #scatterpie is pretty jank, values have to be in column called "value"
    ungroup()
}


dat.cur = dat_filt(dat, code.cur = code.cur, dat.source = "NFJ")
dat.spie = dat_filt_pie(dat, code.cur = code.cur,
                        dat.source = "NFJ", rnd.digits = 0)

gp.nfj = ggplot() +
  geom_sf(data = state.map.data) + 
  geom_point(data = dat.cur,
             aes(x = lon.rnd, y = lat.rnd, size = nyear, 
                 color = nyear >= 6))+
  ggtitle("Naba Fourth of July Count\nAcross-year coverage")+
  xlab("")+
  ylab("")
gp.nfj

## piecharts

ggplot() +
  geom_sf(data = state.map.data) + 
  geom_scatterpie(data = dat.spie,
                  aes(x = lon.rnd, y = lat.rnd), cols = "name",
                  long_format = TRUE, pie_scale = .3)+
  ggtitle("Naba Fourth of July Count, recorded names")+
  xlab("")+
  ylab("")+
  guides(fill=guide_legend("recorded name"))


## MA club demo
## 
dat.cur = dat_filt(dat, code = code.cur, dat.source = "MASSBfly", rnd.digits = 1)

gp.mab = ggplot() +
  geom_sf(data = state.map.data) + 
  geom_point(data = dat.cur,
             aes(x = lon.rnd, y = lat.rnd, size = nobs))+
  ggtitle("MA Butterfly club\nAcross-trip coverage")+
  xlab("")+
  ylab("")+
  xlim(c(-74, -70))+
  ylim(c(41.2, 42.8))
gp.mab

## Pollard demo
## 

dat.cur = dat_filt(dat, code = code.cur, dat.source = "Pollard", rnd.digits = 1)

gp.pol = ggplot() +
  geom_sf(data = state.map.data) + 
  geom_point(data = dat.cur,
             aes(x = lon.rnd, y = lat.rnd, size = nyear, 
                 color = nyear >= 6))+
  ggtitle("Pollard-type data\nAcross-year coverage")+
  xlab("")+
  ylab("")
gp.pol
#snag the most common "name" entry
gp.pointscale = scale_size_continuous(range = c(.5,4),
                                      limits = c(0,50))
name.common = names(rev(sort(table((dat.cur = dat %>% filter(code == code.cur))$name))))[[1]]
gg.multi = ggarrange( text_grob(paste0(code.cur,"\n",name.common), rot = 90, size = 14), 
                      gp.nfj+gp.pointscale, 
                      gp.pol+gp.pointscale, 
                      gp.mab, ncol = 4,
                      widths = c(.1, 1,1,1),
                      common.legend = F)
gg.multi

# gg.multimulti = ggarrange(gg.multi, 
#                           gg.multi, 
#                           gg.multi,
#                           gg.multi, 
#                           gg.multi, ncol = 1, nrow = 4)
# ggexport(gg.multimulti, filename = here("4_res/data-coverage/demo-maps.pdf"),
#          width = 25, height =20)


## Doing all codes---------------

## Looping
## Need to set a threshold for total number of observations
codes.plot = names(table(dat$code)[table(dat$code)>200])
#cut out codes with dashes, which don't correspond to individual species
codes.plot = codes.plot[!grepl("[-]", codes.plot)]
## there are some "TBD" code entries. not useful for us.
codes.plot = codes.plot[codes.plot != "TBD"]

#dataframe to store taxa names
taxa.df = NULL
plots.list = list()
for(i in 1:length(codes.plot)){
# for(i in 1:4){
  code.cur = codes.plot[i]
  
  ## NFJ 
  
  dat.cur = dat_filt(dat, code.cur = code.cur, dat.source = "NFJ")
  gp.nfj = ggplot() +
    geom_sf(data = state.map.data) + 
    geom_point(data = dat.cur,
               aes(x = lon.rnd, y = lat.rnd, size = nyear, 
                   color = nyear >= 6))+
    ggtitle("Naba Fourth of July Count\nAcross-year coverage")+
    xlab("")+
    ylab("")+
    guides(col = guide_legend("6+ Years"),
           size = guide_legend("Years\ncoverage"))
  # gp.nfj
  
  ## MA club 
  dat.cur = dat_filt(dat, code.cur = code.cur, dat.source = "MASSBfly", rnd.digits = 1)
  gp.mab = ggplot() +
    geom_sf(data = state.map.data) + 
    geom_point(data = dat.cur,
               aes(x = lon.rnd, y = lat.rnd, size = nobs))+
    ggtitle("MA Butterfly club\nAcross-trip coverage")+
    xlab("")+
    ylab("")+
    xlim(c(-74, -70))+
    ylim(c(41.2, 42.8))+
    guides(size = guide_legend("# observation\nevents"))
  gp.mab
  
  ## Pollard 
  dat.cur = dat_filt(dat, code.cur = code.cur, dat.source = "Pollard", rnd.digits = 1)

  gp.pol = ggplot() +
    geom_sf(data = state.map.data) + 
    geom_point(data = dat.cur,
               aes(x = lon.rnd, y = lat.rnd, size = nyear, 
                   color = nyear >= 6))+
    ggtitle("Pollard-type data\nAcross-year coverage")+
    xlab("")+
    ylab("")+
    guides(col = guide_legend("6+ Years"),
           size = guide_legend("Years\ncoverage"))
  
  #snag the most common "name" entry
  names.all = rev(sort(table((dat %>% filter(code == code.cur))$name)))
  names.string = names(names.all)
  names.string[1:length(names.string) %% 3 ==0] =  paste0(names.string[1:length(names.string) %% 3 ==0], "\n")
  
  taxa.df = rbind(taxa.df, data.frame(code = code.cur, 
                                      names.frequent = paste(names(names.all[names.all>10]),
                                                             collapse = " | ")))
  # join plots:
  # 
  gp.pointscale = scale_size_continuous(range = c(.5,2.5), limits = c(0,50))
  
  text.spec = paste0(dat$common[1], " (",code.cur,")","\n",
                     paste(names.string,
                           collapse = " | "), "\n", 
                     sum((dat %>% filter(code == code.cur))$count), " counted; ",
                     nrow(dat %>% filter(code == code.cur)), " observation events")
  plots.list[[i]] = ggarrange( text_grob(text.spec, rot = 90, size = 14), 
                               gp.nfj+gp.pointscale, 
                               gp.pol+gp.pointscale, 
                               gp.mab, ncol = 4,
                               widths = c(.1, 1,1,1))
}
write.csv(taxa.df,
          "4_res/data-coverage/taxa-list.csv",
          row.names = FALSE)

gg.all = ggarrange(plotlist = plots.list, nrow = 4, ncol = 1)
ggexport(gg.all, filename = here("4_res/data-coverage/maps.pdf"),
         width = 25, height =20, res = 150)
