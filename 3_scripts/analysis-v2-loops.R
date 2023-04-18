## looping over species for creating diagnostics for 3/3/2023

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
n.threads.use = 7
regions.dict = read.csv(here("2_data_wrangling/FWS-regions-by-state.csv"))

## read in the range areas for combining into a species total
range.areas = read_csv(here("2_data_wrangling/range-area-by-regions.csv"))


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
              "CHLOLEA", "COLEU1", "ENOCRE", "ENOPOR","HESLEO", "LERACC", "SATFAV", "URBPR3"
)
codes.use = codes.use[-(1:max(which(codes.use %in% code.prob)))]

if(!dir.exists(here("4_res/v2-plots"))){dir.create(here("4_res/v2-plots"))}
if(!dir.exists(here("4_res/v2-plots/summary"))){dir.create(here("4_res/v2-plots/summary"))}


#### Loop ------------
for(cur.code in codes.use){
# for(cur.code in codes.use[1:3]){
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
  
  ## summarize information by region to help avoid fitting regions we have less than 10 years of non-inferred data for.
  dat.filt = dat.spec %>% 
    group_by(region, year) %>% 
    summarize(n.noninferred = sum(inferred == FALSE)) %>% 
    group_by(region) %>% 
    summarize(nyear.good = sum(n.noninferred>0)) %>% 
    ungroup() %>% 
    filter(nyear.good >= 10)
  
  ## remove those years with insufficient data
  dat.spec = dat.spec %>% 
    filter(region %in% dat.filt$region) %>% 
    mutate(event.typefac = as.factor(event.type))
  
  ## if we have more than 1 data source and more than 1 event type
  
  ## identify pheno windows - periods of time in which the majority of events occur
  ## so that we are not estimating abundance on inferred regions
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
  
  
  
  # 
  ## create baseline formula: include regional effects if more than 1 region represented
  if(length(unique(dat.spec$region))>1){
    form = "count ~-1 +  s(doy, by = regionfac, bs = \"cc\", k = 10) + regionfac + year:regionfac"
  }else{
    form = "count ~ -1 + s(doy, bs = \"cc\", k = 10) + year"
  }
  ## Iff multiple sources, include sourcefac
  if(length(unique(dat.spec$source))>1){
    form = paste0(form, " + sourcefac")
  }
  ## Iff effort.universal.type is not JUST site-based (in which case effort.universal is all 0s), include effort.universal
  if(!all(dat.spec$effort.universal.type == "site-based")){
    ## Iff multiple event types, include effort.universal.type interaction, otherwise no
    if(length(unique(dat.spec$effort.universal.type))==1){
      form = paste0(form, " + effort.universal")  
    }else{
      form = paste0(form, " + effort.universal:effort.universal.type")  
    }
  }
  
  ## Add random effect
  form = paste0(form, " + s(site.refac, bs = \"re\")")
  
  
  ##fit the model, time it
  cat(paste0("Fitting ", cur.code, "\n"))
  init.time = proc.time()
  out = bam(as.formula(form), 
            method="fREML", 
            knots = list(doy = c(0.5, 364.5)),
            family = "nb",
            discrete = TRUE,
            nthreads = n.threads.use,
            data = dat.spec)
  proc.time()-init.time
  
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
  ## if we only had 1 region used, region will show as blank
  ## let's update to list the region used
  if(nrow(coefs.temp)==1 & all(coefs.temp$region=="")){coefs.temp$region = unique(dat.spec$region)}
  ## Note that this is intentionally fragile - if we fiddle with model construction, this will need
  ## to be changed, and I'd rather an error than misleading results.
  
  ## extract regional densities
  coef(out)[1]
  dens.coef = coef(out)
  dens.coef = dens.coef[grepl("^regionfac", names(dens.coef))]
  dens.coef = dens.coef[!grepl("[:]", names(dens.coef))]
  ##NOTE: IF WE HAVE AN INTERCEPT, WE NEED TO ADD IT HERE!!
  dens.est = exp(dens.coef)
  dens.est = dens.est/sum(dens.est)
  names(dens.est) = gsub("regionfac", "", names(dens.est))
  dens.df = data.frame(density = dens.est, 
                       region = names(dens.est), row.names = NULL)
  dens.df$measure = "density"
  
  
  ## when we have more than 1 region represented, we should use the delta method to 
  ## estimate the species trend
  if(nrow(coefs.temp)>1){
    out.vcov = vcov(out, freq = TRUE)
    out.vcov = out.vcov[grepl("year", rownames(out.vcov)), grepl("year", colnames(out.vcov))]
    out.coef = coef(out)
    out.coef = out.coef[grepl("year", names(out.coef))]
    ## extracting names for weighting
    coef.names = gsub("year", "", names(out.coef))
    coef.names = gsub("regionfac", "", coef.names)
    coef.names = gsub(":", "", coef.names)
    ##pulling out areas for weights, turning into relative representation
    cur.areas = range.areas %>% 
      filter(code == cur.code)
    ## vector form for ease of use
    areas.vec = cur.areas$km2
    names(areas.vec) = cur.areas$region
    coef.weights = areas.vec[coef.names]
    coef.weights = coef.weights/sum(coef.weights)
    weights.df = data.frame(weight = coef.weights,
                            region = names(coef.weights))
    row.names(weights.df) = NULL
    weights.df$measure = "area"
    
    ## turn dens.est into correct order for weight multiplication
    dens.est = dens.est[coef.names]
    
    weights.full = coef.weights*dens.est/sum(coef.weights*dens.est)
    
    ## delta method
    form.delta = paste0(paste0(paste0("x", 1:length(coef.weights)), "*", weights.full), collapse = " + ")
    form.delta = paste0("~ ", form.delta)
    
    est.total = sum(weights.full * out.coef)  
    est.ses = deltamethod(as.formula(form.delta),
                          mean = out.coef, 
                          cov = out.vcov)
  }else{
    est.total = coefs.temp$estimate
    est.ses = coefs.temp$SE
    weights.df = data.frame(weight = 1, region = dat.spec$region[1], dummy = "dummy")
  }
  coefs.temp = coefs.temp %>% 
    add_row(coefficient = "total", estimate = est.total, SE = est.ses, tval = NA, Pval = NA, region = "continental")
  
  cat("making figures\n")
  ##plot of results
  gp.trends = ggplot(coefs.temp %>% filter(region != "continental"), aes(x = region, y = estimate))+
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
  gp.continental = ggplot(coefs.temp %>% filter(region == "continental"), aes(x = region, y = estimate))+
    geom_hline(yintercept = 0, linetype = 2)+
    geom_point(size = 5)+
    geom_segment(aes(y = estimate - 1.96*SE,
                     yend = estimate + 1.96*SE,
                     xend = region), linewidth = 2)+
    ggtitle("Specieswide estimate")+
    xlab("")+
    ylab("Growth rate")+
    ylim(layer_scales(gp.trends)$y$range$range)+ #let's grab the ylimits from gp.trends
    theme.larger
  gp.continental
  
  fullweights.df = full_join(weights.df %>% 
                               select(-measure) %>% 
                               rename(rel.area = weight),
                             dens.df %>% 
                               select(-measure))
  fullweights.df$weight = fullweights.df$rel.area * fullweights.df$density
  fullweights.df$weight = fullweights.df$weight/sum(fullweights.df$weight)
  weights.plot = pivot_longer(fullweights.df,
                              cols = c("rel.area", "density", "weight"),
                              names_to = "measure")
  weights.plot$measure = gsub("weight", "total relative weighting",  
                              weights.plot$measure)
  weights.plot$measure = gsub("rel.area", "range area per\nregion (relative)",
                              weights.plot$measure)
  weights.plot$measure = gsub("density", "species density\n(per site)",   
                              weights.plot$measure)
  
  gp.weights = ggplot(weights.plot, aes(fill = region, y = value, x = measure))+
    geom_col(position = "stack")+
    xlab("")+
    ylab("Proportion")+
    ggtitle("Regional weights: area, density, total weight")+
    theme.larger
  
  
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
  # plot_grid(gp.trends, gp.phenos, ncol=1)
  
  dat.window = dat.window %>% 
    select(region, n.noninferred)
  
  coefs.temp = full_join(dat.window, coefs.temp)
  fig.comb = (gp.trends + gp.continental + plot_layout(widths = c(4,1)))/
    gp.weights/
    gp.phenos
  
  
  write_csv(coefs.temp,
            file = here(paste0("4_res/v2-plots/", cur.code, "-trend.csv")))
  ggsave(here(paste0("4_res/v2-plots/", cur.code,"-trends.pdf")),
         fig.comb, width = 22, height = 26)
  saveRDS(fig.comb, file = here(paste0("4_res/v2-plots/", cur.code,"-fig.RDS")))
  
}

### collect and save single file; make regional histograms
res.files = list.files(here("4_res/v2-plots/"))
res.files = res.files[grepl(".*-trend[.]csv", res.files, )]

coefs.all = NULL
for(cur.file in res.files){
  cur.coef = read_csv(here("4_res/v2-plots/", cur.file))
  cur.code = gsub("-trend[.]csv", "", cur.file)
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
