## looping over species for creating diagnostics for 3/3/2023

## for working in base R:
setwd("C:/repos/status-of-butterflies-analysis")

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

year.min = 1990 #cut off years before this
yr.region.min = 10 # need at least this many years represented in non-inferred data 
# to use a region.
n.threads.use = 3
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))


# cur.code = "POAZAB"

## identify species
dat = qread(paste0("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv"))
#count the number of years, sites, and observation events with > 0 butterflies reported (per species)
dat = dat %>% 
  filter(code != "TBD") %>% 
  filter(code != "BFLY") %>% 
  filter(code != "NONE")
dat.sum = dat %>% 
  filter(presence == 1) %>% 
  group_by(code) %>% 
  summarise(nyear = length(unique(year)),
            nsite = length(unique(site)),
            nsource = length(unique(source)),
            neffort.type = length(unique(effort.universal.type)),
            neffort.type.rich =  sum(table(effort.universal.type)>10),
            nobs = n()
  ) %>% 
  ungroup()
dat.use = dat.sum %>% ## identifying thresholds for minimum years of data, minimum sites, minimum obs 
  filter(nobs >= 100) %>% 
  filter(nsite >= 10 ) %>% 
  filter(nsource > 1 ) %>% 
  filter(neffort.type > 1)
dat.use = dat.use[!(grepl("-", dat.use$code)),]
codes.use = dat.use$code
maps.list = list.files(here("2_data_wrangling/range-maps/relabeled/"))
maps.list = grep("*.kml", maps.list, value = TRUE)
maps.list = gsub(".kml", "", maps.list)
cat("The following codes have no range map, and will be excluded:\n")
print(codes.use[!(codes.use %in% maps.list)], collapse = ", ")
codes.use = codes.use[codes.use %in% maps.list]
## checking for at least 1 region with 10+years
dat.reg = dat %>% 
  filter(code %in% codes.use) %>% 
  inner_join(regions.dict) %>% 
  filter(count > 0 ) %>% 
  group_by(code, region) %>% 
  summarize(nyear = length(unique(year))) %>% 
  group_by(code) %>% 
  summarize(nyear.best = max(nyear)) %>% 
  filter(nyear.best >= yr.region.min)
codes.use = dat.reg$code
dat = dat.use = NULL

## if we break, list them here
code.prob = c("BOLOCHA", "CELNEG", "CHLOLAC", 
              "CHLOLEA", "COLEU1", ##Everything after this failed 
              #even with the new single-region toggling ability
              "ENOCRE", "ENOPOR","HESLEO", "LERACC", "SATFAV", "URBPR3"
              )
codes.use = codes.use[-(1:max(which(codes.use %in% code.prob)))]

### Loop ------------
for(cur.code in codes.use){
  if(TRUE){
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
  
  if(length(unique(dat.spec$source))>1 &
     length(unique(dat.spec$event.type))>1){
  
  #### identify pheno windows
  dat.window = dat.spec %>% 
    group_by(region) %>% 
    filter(count > 0 ) %>% 
    summarize(pheno.start = quantile(doy, probs = 0.005),
              pheno.end = quantile(doy, probs = 0.995),
              n.noninferred = sum(inferred == FALSE)) %>% 
    ungroup()
  
  # hist((dat.spec %>% filter(inferred == FALSE, region == "Pacific Southwest"))$doy, breaks=40)
  
  # ggplot(dat.spec %>% filter(inferred == FALSE, region == "Pacific Southwest"),
  # aes(x = doy, y = count)) +
  # geom_point()+
  # geom_smooth()
  
  
  ## baseline function -----
  cat(paste0("Fitting ", cur.code, "\n"))
  init.time = proc.time()
  if(length(unique(dat.spec$region))>1){
  out = bam(count ~ s(doy, by = regionfac, 
                      bs = "cc", k = 10) + #phenology
              regionfac + year:regionfac + #variation across regions
              sourcefac + effort.universal:effort.universal.type + #accounting for effort and source
              s(site.refac, bs = "re"), #random effect to account for effort, biology, pseudoreplication
            method="fREML", 
            knots = list(doy = c(0.5, 364.5)),
            family = "nb",
            discrete = TRUE,
            nthreads = n.threads.use,
            data = dat.spec)
  proc.time()-init.time
  }else{
    out = bam(count ~ s(doy, bs = "cc", k = 10) + #phenology
                year + #variation across regions
                sourcefac + effort.universal:effort.universal.type + #accounting for effort and source
                s(site.refac, bs = "re"), #random effect to account for effort, biology, pseudoreplication
              method="fREML", 
              knots = list(doy = c(0.5, 364.5)),
              family = "nb",
              discrete = TRUE,
              nthreads = n.threads.use,
              data = dat.spec)
    print(proc.time()-init.time)
  }
  cat("extracting coefficients\n")
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
  
  cat("making figures\n")
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
  
  
  #add region and region label for plotting
  dat.pred$region = as.character(dat.pred$regionfac)
  dat.pred = inner_join(dat.pred, dat.window)
  dat.pred$regionlab = paste0(dat.pred$regionfac, " (N =", dat.pred$n.noninferred, ")")
  #create data with region label for geom_rug
  dat.plot = inner_join(dat.spec, dat.window)
  dat.plot$regionlab = paste0(dat.plot$regionfac, " (N =", dat.plot$n.noninferred, ")")
  
  #add date for plotting
  dat.pred$date = as.Date(dat.pred$doy+1, origin = "1990-01-01")
  dat.plot$date = as.Date(dat.plot$doy+1, origin = "1990-01-01")
  
  gp.phenos = ggplot(dat.pred, aes(x = date, y = count))+
    # geom_point(data = dat.spec, aes(col = sourcefac))+
    geom_path(col = "black")+
    facet_wrap(. ~ regionlab, scales = "free")+
    geom_rug(data = dat.plot %>% filter(inferred == FALSE),
             aes(x=date),
             sides = "b", inherit.aes = FALSE)+
    scale_x_date(date_labels = "%b %e")+
    theme.larger+
    xlab("")+
    ylab("estimated activity")+
    theme(axis.text.x = element_text(size = rel(.7)))
  plot_grid(gp.trends, gp.phenos, ncol=1)
  
  dat.window = dat.window %>% 
    select(region, n.noninferred)
  
  coefs.temp = full_join(dat.window, coefs.temp)
  fig.comb = plot_grid(gp.trends, gp.phenos, ncol=1)
  
  write_csv(coefs.temp,
            file = here(paste0("4_res/v2-plots/", cur.file, "-trend.csv")))
  ggsave(here(paste0("4_res/v2-plots/", cur.code,"-trends.pdf")),
         fig.comb, width = 16, height = 16)
  saveRDS(fig.comb, file = here(paste0("4_res/v2-plots/", cur.code,"-fig.RDS")))
  }else{
    cat("insufficient data after filtering. Skipping. /n")
  }
  
}

### collect and save single file; make regional histograms
res.files = list.files(here("4_res/v2-plots/"))
res.files = res.files[grepl(".*-trends[.]csv", res.files, )]

coefs.all = NULL
for(cur.file in res.files){
  cur.coef = read_csv(here("4_res/v2-plots/", cur.file))
  cur.code = gsub("-trends[.]csv", "", cur.file)
  cur.coef$code = cur.code
  coefs.all = rbind(coefs.all, cur.coef)
}
dict.common = read_csv(here("2_data_wrangling/dictionaries/code-to-common-name.csv"))

coefs.all = left_join(coefs.all, dict.common)

write_csv(coefs.all, 
          here("4_res/v2-plots/summary/version-2-regional-trends.csv"))

coefs.all = coefs.all %>% 
  filter(!is.na(region))

region.means = coefs.all %>% 
  group_by(region) %>% 
  summarize(trend.mean = mean(estimate, na.rm = T))

gp = ggplot(coefs.all)+
  geom_histogram(aes(x = estimate), bins = 50)+
  facet_wrap(. ~ region) + 
  geom_vline(data=region.means, aes(xintercept = trend.mean), col = 'blue')+
  geom_vline(xintercept = 0, linetype=2)+
  theme_bw()+
  theme.larger
ggsave(here("4_res/v2-plots/summary/version-2-regional-trends.pdf"),
       gp,
       width = 16, height = 12)
