## modification of analysis-v2-loops for the all-in-one-hat model.


## for working in base R:
setwd("G:/repos/status-of-butterflies-analysis")

## trying to figure out model running v2.0 - no spatial smooth

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
library(cowplot)
library(lubridate)
library(msm) #for delta method
library(patchwork)
library(terra)
library(tidyterra)
library(RColorBrewer)
source(here("3_scripts/funs.R"))
source(here("3_scripts/funs-fitting.R"))
source(here("3_scripts/model-fitting-function.R"))

gg_quickbar = function(vec, #vector with named entries. Typically table()
                       compress = NULL){ ## Must be null, or between 0 and 1. If not null, 
  ## entries that compress less than `compress` fraction of the total counts are collapsed
  ## into a single category. Useful to simplify viz with many categories
  if(!is.null(compress)){
    stopifnot(compress>=0 & compress<=1)
    vec.perc = vec/sum(vec)
    vec.new = sum(vec[vec.perc < compress])
    names(vec.new) = paste0("others comprising less than ", compress*100," percent")
    vec.new = c(vec.new, vec[vec.perc >= compress])
    vec = vec.new
  }
  dat.plot = data.frame(name = names(vec), 
                        val = as.numeric(vec))
  dat.plot$name = factor(dat.plot$name, levels = names(sort(vec)))
  dat.plot$dummy = ""
  ggplot(dat.plot, aes(fill = name, x = dummy, y = val))+
    geom_bar(position = "stack", stat="identity")+
    xlab("")
}

theme.larger =   theme(axis.title = element_text(size = rel(1.8)),
                       axis.text = element_text(size = rel(1.8)),
                       strip.text = element_text(size = rel(1.8)),
                       plot.title = element_text(size = rel(1.8)),
                       legend.text = element_text(size = rel(1.8)),
                       legend.title = element_text(size = rel(1.8)),
)

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

year.min = 1990 #cut off years before this
year.baseline = 2010
yr.region.min = 10 # need at least this many years represented in non-inferred data 
# to use a region.
n.threads.use = 6
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))
## read in usfws map
usfws = st_read(here("2_data_wrangling/FWS_National_Legacy_Regional_Boundaries/FWS_National_Legacy_Regional_Boundaries.shp"))
usfws = terra::vect(usfws)
usfws = usfws[usfws$REGNAME!="Alaska Region",]

## read in the range areas for combining into a species total
range.areas = read_csv(here("2_data_wrangling/range-area-by-regions.csv"))


## Reading in data, combining all species
dat = qread("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv")
# clear out non-community data, stuff that's too early
dat = dat %>% 
  filter(gather == "community") %>% 
  filter(year >= year.min) %>% 
  inner_join(regions.dict)

## make sure we have storage folder
if(!dir.exists(here("4_res/v2-plots"))){dir.create(here("4_res/v2-plots"))}
if(!dir.exists(here("4_res/v2-plots/one-hat"))){dir.create(here("4_res/v2-plots/one-hat"))}

dat.spec = dat %>% 
  filter(!is.na(count)) %>% 
  filter(year >= year.min) %>% 
  group_by(site, count, year, doy, source, event.type, event.id, duration, party.minutes, effort.universal, effort.universal.type, region) %>% 
  summarize(count = sum(count)) %>% 
  ungroup %>% 
  mutate(sourcefac = as.factor(source),
         effort.universal.type = as.factor(effort.universal.type))
dat.spec$site.refac = sitere_maker(dat.spec)
## handle regions
dat.spec = dat.spec %>% 
  filter(!is.na(region)) %>% 
  mutate(regionfac = as.factor(region))
dat.spec$year = dat.spec$year-year.min

## identify pheno windows - periods of time in which the majority of events occur
## so that we are not estimating abundance on inferred regions
dat.window = dat.spec %>% 
  group_by(region) %>% 
  filter(count > 0 ) %>% 
  summarize(pheno.start = quantile(doy, probs = 0.005),
            pheno.end = quantile(doy, probs = 0.995)) %>% 
  ungroup()

## Visualizations of the data -------------------
gg_quickbar(sort(table(dat.spec$source)), compress = .02)+
  ggtitle("Total representation in the one-hat data")+
  theme.larger+
  ylab("number of records")

## Fitting model  -----
cat("Fitting model\n")
init.time = proc.time()
out = bam(count ~ -1 + s(doy, by = regionfac, bs = "cc", k = 10) + 
            regionfac + year:regionfac + 
            sourcefac + effort.universal:effort.universal.type + 
            s(site.refac, bs = "re"),
          method="fREML", 
          knots = list(doy = c(0.5, 364.5)),
          family = "nb",
          discrete = TRUE,
          nthreads = n.threads.use,
          data = dat.spec)
proc.time()-init.time

## Phenology predictions, find density estimates (on a per-site-reporting basis)
cur.source = "NFJ"
cur.region = unique(dat.spec$regionfac)
h = 0.1
dat.pred = as.data.frame(expand.grid(year = 0:max(dat.spec$year),
                                     doy = seq(0, 365, by = h),
                                     regionfac = cur.region,
                                     sourcefac = cur.source,
                                     effort.universal = 0,
                                     event.typefac = "NFJ",
                                     effort.universal.type = "duration",
                                     site.refac = "totally new site"))
#filter to phenological windows
dat.pred1 = NULL
for(cur.regionfac in unique(dat.spec$regionfac)){
  cur.window = dat.window %>% filter(region == cur.regionfac)
  cur.dat = dat.pred %>% 
    filter(regionfac == cur.regionfac) %>% 
    filter(doy >= cur.window$pheno.start,
           doy <= cur.window$pheno.end)
  dat.pred1 = rbind(dat.pred1, cur.dat)
}
dat.pred = dat.pred1

dat.pred$count = predict(out, newdata = dat.pred, type = 'response',
                         discrete = FALSE)
dens.df = dat.pred %>% 
  group_by(regionfac, year) %>% 
  summarize(density = sum(count)*h) %>% 
  ungroup() %>% 
  rename(region = regionfac) %>% 
  mutate(region = as.character(region)) %>% 
  mutate(measure = "density") %>% 
  mutate(density = density/sum(density)) %>%  #rescale to sum to 1
  mutate(year = year + year.min)

dens.baseline = dens.df %>% 
  filter(year == year.baseline) %>% 
  mutate(density = density/sum(density))
dens.est = dens.baseline$density
names(dens.est) = dens.baseline$region

## Calculate "abundance", proportional to #s of butterflies per region per year.
# get area in square kilometer
cur.areas = range.areas %>%
  filter(code == "total region")
# vector form for ease of use
areas.vec = cur.areas$km2
names(areas.vec) = cur.areas$region
# we're going to scale to naba circle areas: 
area.nfj = pi*(25/2)^2
abund.df = left_join(dens.df, cur.areas) %>% 
  mutate(abund = density * km2 / area.nfj)
## total abundance across regions
abund.tot = abund.df %>% summarize(abund = sum(abund), .by = year)

cat("extracting coefficients\n")
## getting a clean "region" label for each row
out.sum = summary(out)
temp = out.sum$p.table
coefs.temp = data.frame(coef = rownames(temp), temp)
rownames(coefs.temp) = NULL
coefs.temp = coefs.temp %>% 
  filter(str_detect(.$coef, "year"))
coefs.temp$region = gsub("year", "", coefs.temp$coef)
coefs.temp$region = gsub("regionfac", "", coefs.temp$region)
coefs.temp$region = gsub(":", "", coefs.temp$region)
coefs.temp = coefs.temp %>% 
  rename(coefficient = coef,
         estimate = Estimate,
         SE = "Std..Error",
         tval = t.value,
         Pval = "Pr...t..")
coefs.temp

## when we have more than 1 region represented, we should use the delta method to 
## estimate the species trend
out.vcov = vcov(out, freq = TRUE)
out.vcov = out.vcov[grepl("year", rownames(out.vcov)), grepl("year", colnames(out.vcov))]
out.coef = coef(out)
out.coef = out.coef[grepl("year", names(out.coef))]
## extracting names for weighting
coef.names = gsub("year", "", names(out.coef))
coef.names = gsub("regionfac", "", coef.names)
coef.names = gsub(":", "", coef.names)

## looping over different years for alternate weighting schemes
for(cur.year in c(1990, 2000, 2010, 2020)){
  ##pulling out areas for weights, turning into relative representation
  abund.cur = abund.df %>% 
    filter(year == cur.year)
  weights.full = abund.cur$abund
  weights.full = weights.full/sum(weights.full)
  names(weights.full) = abund.cur$region
  
  ## delta method
  form.delta = paste0(paste0(paste0("x", 1:length(weights.full)), "*", weights.full), collapse = " + ")
  form.delta = paste0("~ ", form.delta)
  
  est.total = sum(weights.full * out.coef)  
  est.ses = deltamethod(as.formula(form.delta),
                        mean = out.coef, 
                        cov = out.vcov)
  coefs.temp = coefs.temp %>% 
    add_row(coefficient = "total", estimate = est.total, SE = est.ses, tval = NA, Pval = NA, region = as.character(cur.year))
}

cat("making figures\n")

##plotting results ----------------------

## Blobs vs region map
gp.map = ggplot() +
  geom_spatvector(data = usfws, aes(fill = REGNAME))+
  xlim(c(-130, -65))+
  ylim(c(24, 50))+
  ggtitle("FWS Regions")+
  theme.larger+
  scale_fill_brewer(palette = "Dark2")

## Abundance labels
abund.init = abund.tot$abund[abund.tot$year==1990]
abund.fin = abund.tot$abund[abund.tot$year==2020]
change.perc = round((abund.fin/abund.init-1)*100)
change.lab = "increase"
if(change.perc <0){change.lab = "decline"}

gp.abund = ggplot(abund.df, aes(x = year, y = abund, col = region))+
  geom_path(linewidth = 1)+
  geom_path(data = abund.tot, col = "black", linewidth = 1.5)+
  ggtitle(paste0(change.perc, "% ", change.lab, " in abundance between 1990 and 2020\n(black line is species total)"))+
  ylab("Abundance index (AUC) x region area (scaled to NFJ circles)")+
  xlab("")+
  theme.larger
gp.abund

# realized growth rate
gr.path = data.frame(gr = diff(log(abund.tot$abund)),
                     year = abund.tot$year[-1]-.5)
gp.grpath = ggplot(gr.path, aes(x = year, y = gr))+
  geom_path(linewidth = 0.8)+
  ggtitle("Realized growth rate\n(sum of exponentials, so gamma shape)")+
  ylab("Instantaneous growth rate\n(diff(log(abundance)))")+
  xlab("")+
  theme.larger

coefs.regional = coefs.temp %>% filter(coefficient != "total")
coefs.regional$region = gsub(" ", "\n", coefs.regional$region)
coefs.regional$region = gsub("-", "\n", coefs.regional$region)
coefs.regional$regionfac = factor(coefs.regional$region, levels = c("Pacific", "Pacific\nSouthwest", "Mountain\nPrarie",
                                                                    "Southwest", "Midwest", "Southeast", "Northeast")) 

## Estimated growth rates
gp.trends = ggplot(coefs.regional, aes(x = regionfac, y = estimate))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_point(size = 5)+
  geom_segment(aes(y = estimate - 1.96*SE,
                   yend = estimate + 1.96*SE,
                   xend = region), linewidth = 2)+
  ggtitle("All butterflies (lumping species)")+
  xlab("")+
  ylab("Growth rate")+
  scale_x_discrete(drop = FALSE)+
  theme.larger
gp.trends

## Overall growth rate
gp.continental = ggplot(coefs.temp %>% filter(coefficient == "total"), aes(x = region, y = estimate))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_point(size = 5)+
  geom_segment(aes(y = estimate - 1.96*SE,
                   yend = estimate + 1.96*SE,
                   xend = region), linewidth = 2)+
  ggtitle("Continental estimate estimate")+
  xlab("Year used for abundance weighting")+
  ylab("Growth rate")+
  ylim(layer_scales(gp.trends)$y$range$range)+ #let's grab the ylimits from gp.trends
  theme.larger
gp.continental

## Weightings

weights.df = abund.df %>% 
  filter(year == year.baseline) %>% 
  mutate(density = density/sum(density),
         area = km2/sum(km2),
         abund = abund/sum(abund)) %>% 
  select(-measure)
weights.plot = pivot_longer(weights.df,
                            cols = c("area", "density", "abund"),
                            names_to = "measure")
weights.plot$measure = gsub("abund", "total relative weighting",  
                            weights.plot$measure)
weights.plot$measure = gsub("area", "range area per\nregion (relative)",
                            weights.plot$measure)
weights.plot$measure = gsub("density", "species density\n(per site)",   
                            weights.plot$measure)

gp.weights = ggplot(weights.plot, aes(fill = region, y = value, x = measure))+
  geom_col(position = "stack")+
  xlab("")+
  ylab("Proportion")+
  ggtitle("Regional weights: area, density, total weight")+
  theme.larger





## Looking at phenology curves.

#add region and region label for plotting
dat.pred$region = as.character(dat.pred$regionfac)
dat.pred = inner_join(dat.pred, dat.window)
#create data with region label for geom_rug
dat.plot = inner_join(dat.spec, dat.window)

#add date for plotting
dat.pred$date = as.Date(dat.pred$doy+1, origin = "1990-01-01")
dat.plot$date = as.Date(dat.plot$doy+1, origin = "1990-01-01")

gp.phenos = ggplot(dat.pred %>% filter(year == year.baseline-year.min), aes(x = date, y = count))+
  # geom_point(data = dat.spec, aes(col = sourcefac))+
  geom_path(col = "black")+
  facet_wrap(. ~ region, scales = "free_y")+
  geom_rug(data = dat.plot,
           aes(x=date),
           sides = "b", inherit.aes = FALSE)+
  scale_x_date(date_labels = "%b %e")+
  theme.larger+
  xlab("")+
  ylab("estimated activity")+
  ggtitle(paste0("Activity curves in ",year.baseline))+
  theme(axis.text.x = element_text(size = rel(.7)))
# plot_grid(gp.trends, gp.phenos, ncol=1)


fig.comb = (gp.map /
              gp.abund / 
              (gp.trends + gp.continental + plot_layout(widths = c(4,1.2)))/
              gp.weights/
              gp.grpath /
              gp.phenos) + plot_layout(heights = c(1, .8, .8, .8, .8, .8))


write_csv(coefs.temp,
          file = here("4_res/v2-plots/one-hat/all-butterflies-trend.csv"))
write_csv(abund.df,
          file = here("4_res/v2-plots/one-hat/all-butterflies-regional-abund.csv"))
ggsave(here(paste0("4_res/v2-plots/one-hat/all-butterflies-trends.pdf")),
       fig.comb, width = 20, height = 45)
saveRDS(fig.comb, file = here("4_res/v2-plots/one-hat/all-butterflies-fig.RDS"))
