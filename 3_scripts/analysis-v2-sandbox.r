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
source(here("3_scripts/funs.R"))
source(here("3_scripts/funs-fitting.R"))
source(here("3_scripts/model-fitting-function.R"))

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

cur.code = "ATACAM"
year.min = 1990 #cut off years before this
yr.region.min = 10 # need at least this many years represented in non-inferred data 
# to use a region.

n.threads.use = 2

if(FALSE){
  make_dataset(cur.code, use.range = T,
               infer.messy.levels = c("GENUS", "SUBFAMILY", "FAMILY", "COMPLEX"),
               name.pretty = NULL)
}
dat.spec = read_csv(paste0("2_data_wrangling/cleaned by code/", cur.code, ".csv"))
dat.spec = dat.spec %>% 
  filter(!is.na(count)) %>% 
  filter(year >= year.min)
dat.spec = dat.spec %>% 
  mutate(sourcefac = as.factor(source),
         effort.universal.type = as.factor(effort.universal.type))
dat.spec$site.refac = sitere_maker(dat.spec)
## handle regions
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))
dat.spec = left_join(dat.spec, regions.dict, by = "state")
dat.spec = dat.spec %>% 
  filter(!is.na(region)) %>% 
  mutate(regionfac = as.factor(region))
dat.spec$year = dat.spec$year-1990

dat.filt = dat.spec %>% 
  group_by(region, year) %>% 
  summarize(n.noninferred = sum(inferred == FALSE)) %>% 
  group_by(region) %>% 
  summarize(nyear.good = sum(n.noninferred>0)) %>% 
  ungroup() %>% 
  filter(nyear.good >= 10)

dat.spec = dat.spec %>% 
  filter(region %in% dat.filt$region) %>% 
  mutate(event.typefac = as.factor(event.type))

#### identify pheno windows
dat.window = dat.spec %>% 
  group_by(region) %>% 
  filter(count > 0 ) %>% 
  summarize(pheno.start = quantile(doy, probs = 0.005),
            pheno.end = quantile(doy, probs = 0.995),
            n.noninferred = sum(inferred == FALSE)) %>% 
  ungroup()

hist((dat.spec %>% filter(inferred == FALSE, region == "Pacific Southwest"))$doy, breaks=40)

ggplot(dat.spec %>% filter(inferred == FALSE, region == "Pacific Southwest"),
       aes(x = doy, y = count)) +
  geom_point()+
  geom_smooth()


## baseline function -----

init.time = proc.time()
out = bam(count ~ s(doy, by = regionfac, 
                        bs = "cc", k = 10) + #phenology
            regionfac + year:regionfac + #variation across regions
            event.typefac + effort.universal:effort.universal.type + #accounting for effort and source
            s(site.refac, bs = "re"), #random effect to account for effort, biology, pseudoreplication
          method="fREML", 
          knots = list(doy = c(0.5, 364.5)),
          family = "nb",
          discrete = TRUE,
          nthreads = n.threads.use,
          data = dat.spec)
proc.time()-init.time

## simple version of the model
init.time = proc.time()
out.simple = bam(count ~ 0 + s(doy, bs = "cc", k = 10) + #phenology
            regionfac + year:regionfac + #variation across regions
              event.typefac + effort.universal:effort.universal.type + #accounting for effort and source
            s(site.refac, bs = "re"), #random effect to account for effort, biology, pseudoreplication
          method="fREML", 
          knots = list(doy = c(0.5, 364.5)),
          family = "nb",
          discrete = TRUE,
          nthreads = n.threads.use,
          data = dat.spec)
proc.time()-init.time


## extract coefs
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

##plot of results
gp.trends = ggplot(coefs.temp, aes(x = region, y = estimate))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_point(size = 5)+
  geom_segment(aes(y = estimate - 1.96*SE,
                   yend = estimate + 1.96*SE,
                   xend = region), linewidth = 2)+
  ggtitle(paste0(cur.code, ": ", dat.spec$common[1]))+
  xlab("")+
  ylab("Growth rate")+
  theme.larger
gp.trends
## what if we use the simple version

## extract coefs
out.simple.sum = summary(out.simple)
temp = out.simple.sum$p.table
coefs.simple.temp = data.frame(coef = rownames(temp), temp)
rownames(coefs.simple.temp) = NULL
coefs.simple.temp = coefs.simple.temp %>% 
  filter(str_detect(.$coef, "year"))
coefs.simple.temp$region = gsub("year", "", coefs.simple.temp$coef)
coefs.simple.temp$region = gsub("regionfac", "", coefs.simple.temp$region)
coefs.simple.temp$region = gsub(":", "", coefs.simple.temp$region)
coefs.simple.temp = coefs.simple.temp %>% 
  rename(coefficient = coef,
         estimate = Estimate,
         SE = "Std..Error",
         tval = t.value,
         Pval = "Pr...t..")
coefs.simple.temp

##plot of results
gp.trends.simple = ggplot(coefs.simple.temp, aes(x = region, y = estimate))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_point(size = 5)+
  geom_segment(aes(y = estimate - 1.96*SE,
                   yend = estimate + 1.96*SE,
                   xend = region), linewidth = 2)+
  ggtitle(paste0(cur.code, ": ", dat.spec$common[1]))+
  xlab("")+
  ylab("Growth rate")+
  theme.larger
gp.trends.simple

## Looking at phenology curves.

cur.source = "NFJ"
cur.region = unique(dat.spec$regionfac)
dat.pred = as.data.frame(expand.grid(year = 10,
                                     doy = seq(0, 365, by = .1),
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
# ggplot(dat.pred, aes(x = doy, y = count))+
#   geom_point(data = dat.spec, aes(col = sourcefac))+
#   geom_path(col = "black")+
#   facet_wrap(. ~ regionfac, scales = "free")+
#   theme.larger

gp.phenos = ggplot(dat.pred, aes(x = doy, y = count))+
  # geom_point(data = dat.spec, aes(col = sourcefac))+
  geom_path(col = "black")+
  facet_wrap(. ~ regionfac, scales = "free")+
  theme.larger
gp.phenos

##pheno for simples
cur.source = "NFJ"
cur.region = dat.spec$regionfac[1]
dat.pred.simple = as.data.frame(expand.grid(year = 10,
                                     doy = seq(0, 365, by = .1),
                                     regionfac = cur.region,
                                     sourcefac = cur.source,
                                     effort.universal = 0,
                                     event.typefac = "NFJ",
                                     effort.universal.type = "duration",
                                     site.refac = "totally new site"))
dat.pred.simple$count = predict(out.simple, newdata = dat.pred, type = 'response',
                         discrete = FALSE)
gp.pheno.simple =  ggplot(dat.pred.simple, aes(x = doy, y = count))+
  # geom_point(data = dat.spec, aes(col = sourcefac))+
  geom_path(col = "black")+
  theme.larger
gp.pheno.simple

plot_grid(gp.trends, gp.trends.simple+ggtitle("with only single activity curve"), ncol=1)

plot_grid(gp.phenos, gp.pheno.simple, ncol = 1, rel_heights = c(2,1))

anova.gam(out, out.simple)