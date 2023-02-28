library(sf)
library(tidyverse)
library(MASS, include.only = "glm.nb")
library(viridis)
library(here)
source(here("3_scripts/funs.R"))

## Notes for convex hulls:
# ## identifying convex hull of locations -- don't want to predict in weird spaces because SPLINE
# library(geometry)
# shape.hull = convhulln(dat[,c("lon", "lat")])
# 
# loc.plot = as.data.frame(expand.grid(lon = seq(min(dat$lon), max(dat$lon), by = 1),
#                                      lat = seq(min(dat$lat), max(dat$lat), by = 1),
#                                      doy = seq(0,365, by = 1),
#                                      year = unique(dat$year),
#                                      sourcefac = "NFJ"))
# loc.plot = loc.plot[inhulln(shape.hull, as.matrix(loc.plot[, c("lon","lat")])),]
# pts.rect = as.matrix(expand.grid(lon = seq(min(dat$lon), max(dat$lon), by = 1),
#                                  lat = seq(min(dat$lat), max(dat$lat), by = 1)))
# pts.test = pts.rect[inhulln(shape.hull, pts.rect),]

grid_plot_allyr = function(dat, 
                           regions.dict, 
                           dat.constrain,
                           use.range = TRUE, 
                           do.pheno = TRUE){#if TRUE, predict on each day of the year. Otherwise, predict on one day of hte year.
  ## for easy calculating abundance
  if(do.pheno){ ## model includes phenology curve, so predict along that.
    doy.vec = seq(0,365, by = 1)
  } else{
    doy.vec = 180
  }
  loc.pts = as.data.frame(expand.grid(lon = seq(round(min(dat$lon)), round(max(dat$lon)), by = 1),
                                      lat = seq(round(min(dat$lat)), round(max(dat$lat)), by = 1)))
  if(use.range){
    loc.pts = use_range(dat = loc.pts, code.cur = dat$code[1])
  }else if(dat.constrain){
    shape.hull = convhulln(loc.pts[,c("lon", "lat")])
    loc.pts = loc.pts[inhulln(shape.hull, as.matrix(loc.pts[, c("lon","lat")])),]
  }
  
  loc.pts$country = coords_2country(long = loc.pts$lon, lat = loc.pts$lat)
  loc.pts = loc.pts[!is.na(loc.pts$country),]
  loc.pts = loc.pts[loc.pts$country=="USA",]
  loc.pts$state = coords_2state(long = loc.pts$lon, lat = loc.pts$lat)
  loc.pts = left_join(loc.pts, regions.dict, by = "state")
  loc.pts = loc.pts %>% 
    filter(!is.na(region)) %>% 
    rename(regionfac = region)
  loc.pts = loc.pts[loc.pts$regionfac %in% dat$regionfac, ]
  grid.plot = expand_grid(loc.pts[, c("lat","lon", "regionfac")],
                          doy = doy.vec,
                          year = unique(dat$year),
                          effort.universal = 0,
                          effort.universal.type = factor("party.minutes", levels = levels(dat$effort.universal.type)))
  if(all(unique(dat$effort.universal.type == "site-based"))){ ## accounting for case when site-based is the only option
    grid.plot$effort.universal = 1
    effort.universal.type = factor("site-based", levels = levels(dat$effort.universal.type))
  }
  if("site.refac" %in% names(dat)){
    grid.plot$site.refac = factor("dummy.site", levels = levels(dat$site.refac))
  }
  return(grid.plot)
}

grid_plot_oneyr = function(dat, regions.dict, do.pheno = TRUE, use.range = TRUE, dat.constrain = FALSE){
  ## for calculating trends
  loc.pts = as.data.frame(expand.grid(lon = seq(min(dat$lon), max(dat$lon), by = .5),
                                      lat = seq(min(dat$lat), max(dat$lat), by = .5)))
  if(do.pheno){ ## model includes phenology curve, so predict along that.
    doy.vec = seq(0,365, by = 1)
  } else{
    doy.vec = 180
  }
  if(use.range){
    loc.pts = use_range(loc.pts, code.cur = dat$code[1])
  }else if(dat.constrain){
    shape.hull = convhulln(loc.pts[,c("lon", "lat")])
    loc.pts = loc.pts[inhulln(shape.hull, as.matrix(loc.pts[, c("lon","lat")])),]
  }
  loc.pts$country = coords_2country(long = loc.pts$lon, lat = loc.pts$lat)
  loc.pts = loc.pts[!is.na(loc.pts$country),]
  loc.pts = loc.pts[loc.pts$country=="USA",]
  loc.pts$state = coords_2state(long = loc.pts$lon, lat = loc.pts$lat)
  loc.pts = left_join(loc.pts, regions.dict, by = "state")
  loc.pts = loc.pts %>% 
    filter(!is.na(region)) %>% 
    rename(regionfac = region)
  loc.pts = loc.pts[loc.pts$regionfac %in% dat$regionfac, ]
  grid.plot = expand_grid(loc.pts[, c("lat","lon", "regionfac")],
                          doy = doy.vec,
                          effort.universal = 0,
                          effort.universal.type = factor(dat$effort.universal.type[1], 
                                                         #note: we're setting effort to 0 (remember: normalized), 
                                                         #  so we can use any effort type. Just want one that's present. 
                                                         levels = levels(dat$effort.universal.type)))
  if("site.refac" %in% names(dat)){
    grid.plot$site.refac = factor("dummy.site", levels = c(levels(dat$site.refac), "dummy.site"))
  }
  return(grid.plot)
}



##Update: this is the past abundance mapper. It struggles with larger data sets, I ran into issues with RAM.
##It's also inefficient when we know our population trends are linear. 



#Calculates raw trend, and abundance.
#Note: mean abundance for a single site across years, when growth is exactly exponential (ie our models here):
#1/(N+1) \sum_{t = 0}^n N_0 e^{rt}
#(just tex writing for the sum N_0 + N_0*e^r + N_0*e^(r*2)...)
#This has a closed-form solution that doesn't require calculating individual years
#(valuable given memory limitations, size of our data)
#Solution is \frac{N_0}{n+1} \frac{e^(rt + r)-1}{e^r-1}
## Test case:
#n0 = 1000 #initial size of 1000
# r = .05 #growth rate
# n = 50 #number of years
# ## calculate each year and average
# yrs = 0:n
# pop = n0*exp(r*yrs)
# mean(pop)
# ## closed form
# n0/(n+1)*(exp(n*r + r)-1)/(exp(r)-1)

trend_and_abund_calc = function(dat, fit, regions.dict, 
                                fit.family,
                                use.range = TRUE,
                                dat.constrain = FALSE, 
                                do.pheno = TRUE, #if no doy term in the model, set to FALSE for faster calculations.
                                sourcefac = "NFJ"){
  ## Prepping two sequential years for predicting
  loc.plot = grid_plot_oneyr(dat = dat, regions.dict = regions.dict,
                             dat.constrain = dat.constrain, 
                             do.pheno = do.pheno, use.range = use.range)
  loc.plot$sourcefac = sourcefac
  loc.plot.y1 = loc.plot.y2 = loc.plot
  # loc.plot.y1$year  = min(dat$year)
  # loc.plot.y2$year = loc.plot.y1$year + 1
  loc.plot.y1$year  = 2010
  loc.plot.y2$year = loc.plot.y1$year + 1
  
  #predict
  loc.plot.y1$count = predict(fit, newdata = loc.plot.y1, type = "response", discrete = FALSE)
  loc.plot.y2$count = predict(fit, newdata = loc.plot.y2, type = "response", discrete = FALSE)
  
  ## combining
  loc.plot$gr = log(loc.plot.y2$count)-log(loc.plot.y1$count)
  nyear = max(dat$year)-min(dat$year)
  # loc.plot$abund.index = loc.plot.y1$count/(nyear+1) * (exp(nyear * loc.plot$gr)-1)/(exp(loc.plot$gr)-1)
  loc.plot$abund.index = loc.plot.y1$count
  if(do.pheno){
    loc.plot = loc.plot %>% 
      group_by(lon, lat) %>% 
      summarize(gr.med = median(gr),
                abund.index = sum(abund.index)) %>% 
      ungroup()
  }else{
    loc.plot = loc.plot %>% 
      rename(gr.med = gr) %>% 
      mutate(abund.index = abund.index * 365)%>% #making it butterfly-days, equiv to when phenology is measured
      select(lon, lat, gr.med, abund.index)
  }
  abund.tot = sum(loc.plot$abund.index)
  gr.tot = sum(loc.plot$gr.med * (loc.plot$abund.index/abund.tot))
  return(list(loc.plot = loc.plot, global.gr = gr.tot))
}


## new version: takes loc.plot dataframe from trend_and_abund_calc(), data, plots it
trend_plotter = function(dat, loc.plot){ 
  state.map.data <- maps::map('state', fill = TRUE, plot = FALSE) %>%
    st_as_sf()
  dat.pts = viz_filter(dat %>%  select(lon,lat, inferred), reso = 2)
  gp = ggplot()+
    geom_tile(data = loc.plot, 
              aes(x = lon, y = lat, fill = gr.med))+
    geom_jitter(data = dat.pts %>% filter(inferred == TRUE), aes(x = lon, y = lat), col = "coral", size = .2, width = 0.1, height = 0.1)+
    geom_jitter(data = dat.pts %>% filter(inferred == FALSE), aes(x = lon, y = lat), col = "springgreen3", size = .2, width = 0.1, height = 0.1)+
    geom_sf(data = state.map.data, fill = NA)+
    theme_minimal()+
    labs(fill = "Growth rate")+
    ggtitle(paste0(dat$code[1], ": ", dat$common[1]," growth rates\ncoral points: \"inferred\" zeroes; green points: direct observations"))+
    scale_fill_gradient2(low = muted("blue"),
                         mid = "white",
                         high = muted("green"),
                         midpoint = 0)
  return(list(fig = gp, data  = loc.plot))
}


abund_plotter = function(dat, loc.plot){
  state.map.data <- maps::map('state', fill = TRUE, plot = FALSE) %>%
    st_as_sf()
  
  ## visualizing
  dat.pts = viz_filter(dat %>%  select(lon,lat, inferred), reso = 2)
  gp = ggplot()+
    geom_tile(data = loc.plot, aes(x = lon, y = lat, fill = log10(abund.index)))+
    geom_jitter(data = dat.pts %>% filter(inferred == TRUE), aes(x = lon, y = lat), col = "coral", size = .2, width = 0.1, height = 0.1)+
    geom_jitter(data = dat.pts %>% filter(inferred == FALSE), aes(x = lon, y = lat), col = "springgreen3", size = .2, width = 0.1, height = 0.1)+
    geom_sf(data = state.map.data, fill = NA)+
    scale_fill_viridis()+
    theme_minimal()+
    ggtitle(paste0(dat$code[1], ": ", dat$common[1], 
                   " estimated abundance (averaged across all years)\ncoral points: \"inferred\" zeroes; green points: direct observations"))
  res = list(fig = gp)
  return(res)
}


demo_activity_plots = function(dat, fit, regions.dict){ 
  ## identify area of interest - where are butterflies most commonly reported?
  dat.freq = dat %>%
    filter(count > 0) %>%
    mutate(lat.rnd = round(lat),
           lon.rnd = round(lon),
           loq.unq = paste(lat.rnd, lon.rnd)) %>%
    select(loq.unq)
  loq.freq = rev(sort(table(dat.freq$loq.unq)))
  loq.plot = as.numeric(strsplit(names(loq.freq)[1], " ")[[1]])
  names(loq.plot) = c("lat", "lon")
  activity_plotter(dat, fit, regions.dict, lat.plot = loq.plot[["lat"]], lon.plot = loq.plot[["lon"]])
}


activity_plotter = function(dat, fit, regions.dict, lat.plot, lon.plot, 
                            allyears = FALSE,#allyears = TRUE will plot all years in the overall data
                            #  instead of all years with at least 1 data point. Good for plotting at point of max abund.
                            source.adaptive = TRUE #if TRUE, use the source best represented in the locaiton.
                            #  if FALSE, use NJF. Recommend FALSE if allyears = TRUE
){ 
  ## snag observations in that area
  dat.plot = dat %>%
    filter(round(lon) == lon.plot) %>%
    filter(round(lat) == lat.plot)
  state.plot = coords_2state(long = lon.plot, lat = lat.plot)
  region.plot = regions.dict$region[regions.dict$state == state.plot]
  ## predict densities there
  if(allyears){ 
    years.plot = unique(dat$year)
  }else{
    years.plot = unique(dat.plot$year)
  }
  if(source.adaptive){
    sourcefac = names(rev(sort(table(dat.plot$sourcefac))))[1]
  }else{
    sourcefac = "NFJ"
  }
  dat.pred = as.data.frame(expand.grid(doy = seq(0,365, by = .1),
                                       year = years.plot,
                                       lat = lat.plot,
                                       lon = lon.plot,
                                       sourcefac = sourcefac,
                                       regionfac = region.plot))
  if("site.refac" %in% names(dat)){
    dat.pred$site.refac = as.factor("light on data")
  }
  dat.pred$count = predict(fit, newdata = dat.pred, type = "response")
  ## plot
  dat.pt = dat %>%  
    select(doy, count) %>% 
    unique()
  ggplot(dat.pred, aes(x = doy, y = count)) +
    geom_point(data = dat.pt, aes(col = sourcefac))+
    geom_line(col = 'black') +
    facet_wrap( ~ year)+
    ggtitle(paste0(dat$code[1],": ", dat$sommon[1], " at (", lon.plot, ", ", lat.plot,") predicted for ", names(rev(sort(table(dat.plot$sourcefac))))[1]))+
    theme_minimal()
}

NFJ_compare = function(dat, fit, regions.dict, fit.family,
                       use.range = FALSE,
                       nyears = 5,
                       across.doy = TRUE,#if FALSE, will just look at points estimates of doy for predictions. If TRUE, calculates activity-days in that year
                       across.year = TRUE #if FALSE, will treat each year separately. If TRUE, will report the average across all years. TRUE is probably more
                       #faithful as a model validation, as the model isn't structured to capture year-to-year fluctations.
){ 
  dat.comp = dat %>% 
    filter(sourcefac == "NFJ") %>% 
    filter(max(year)-year < nyears)
  if(use.range){
    dat.comp = use_range(dat.comp, code.cur = dat.comp$code[1])
  }
  n.sites = length(unique(dat.comp$site))
  if(across.doy){
    dat.pred = expand_grid(nesting(dat.comp[, c("lat", "lon", "year", "regionfac")]),
                           doy = 0:365)
    if("site.refac" %in% names(dat)){
      dat.pred$site.refac = as.factor("light on data")
    }
    dat.pred$sourcefac = "NFJ"
    dat.pred$count.pred = predict(fit, newdata = dat.pred, type = "response")
    dat.pred = dat.pred %>% 
      group_by(site, year) %>% 
      summarize(count.pred = mean(count.pred)*365) %>% 
      ungroup()
    dat.full = inner_join(dat.comp, dat.pred)
  }else{
    dat.full = dat.comp
    dat.full$count.pred = predict(fit, newdata = dat.full, type = "response")
  }
  if(across.year){
    dat.full = dat.full %>% 
      group_by(site, region, code) %>% #code is here to maintain it for plot title later
      summarize(count = mean(count),
                presence = mean(presence),
                count.pred = mean(count.pred)) %>% 
      ungroup()
  }
  if(fit.family == "binomial"){ #quick handling of count vs presence/absence
    dat.cor = cor.test(dat.full$presence, dat.full$count.pred)
  }else{
    dat.cor = cor.test(dat.full$count, dat.full$count.pred)
  }
  p_pretty = function(p){
    if(is.na(p) | is.null(p)){
      res = "P undefined"
    }else if(p<0.001){
      res = "P < 0.001"
    }else(
      res = paste0("P = ",round(p, 3))
    )
    return(res)
  }
  if(fit.family == "binomial"){
    gp = ggplot(dat.full, aes(x = presence, y = count.pred))+
      geom_point(aes(col = region))+
      geom_smooth(method = lm)+
      # geom_smooth(method = lm, formula = dat.full$count.pred ~ dat.full$count:dat.full$region)+
      ggtitle(paste0(dat.full$code[1],": ", dat$common[1], ", NFJ presence vs predictions, last ", nyears, " years\n",
                     "Correlation of ", round(dat.cor$estimate,4),", ", p_pretty(dat.cor$p.value)))+
      xlab("Actual presence (avg across year)")+
      ylab("Predicted presence (avg across year)")
  }else{
    gp = ggplot(dat.full, aes(x = count, y = count.pred))+
      geom_point(aes(col = region))+
      geom_smooth(method = lm)+
      # geom_smooth(method = lm, formula = dat.full$count.pred ~ dat.full$count:dat.full$region)+
      ggtitle(paste0(dat.full$code[1],": ", dat$common[1], ", NFJ counts vs predictions, last ", nyears, " years\n",
                     "Point = Site (averaged across years), ", n.sites, " total sites\n",
                     "Correlation of ", round(dat.cor$estimate,4),", ", p_pretty(dat.cor$p.value)))+
      xlab("Actual count")+
      ylab("Prediction")
  }
  return(list(fig = gp, cor = round(dat.cor$estimate,4)))
}

NFJ_regional_trends = function(dat, regions.dict, use.range = FALSE){
  ## calculate states as reminder for each region
  regions.df = data.frame(region = unique(dat$region))
  regions.df$label = ""
  for(i in 1:nrow(regions.df)){
    states = regions.dict[regions.dict$region==regions.df$region[i], "state"]
    states[-length(states)] = paste0(states[-length(states)], ",")
    states[1:length(states) %% 3 ==0] = paste0(states[1:length(states) %% 3 ==0], "\n")
    regions.df$regionlabel[i]=paste0(regions.df$region[i],"\n", paste0(states, collapse = " "))
  }
  dat.plot = dat %>% 
    filter(source=="NFJ") %>% 
    left_join(regions.df, by = "region")
  ggplot(dat.plot, aes(x = year, y = count))+
    geom_point(aes(y = count +1), position=position_jitter(width=.3, height=0),
               shape = 1)+
    facet_wrap(~ regionlabel, scale = "free_y")+
    ylab("Count")+
    geom_smooth(method = "glm.nb")+
    scale_y_log10()+
    ggtitle(paste0(dat$code[1],": ", dat$sommon[1]," Trends in NFJ data, log scale, by region\nline: glm.nb; points: count+1"))
}

# #### Old versions of trend + abund calculations from before optimization.
# 
# trend_plotter_old = function(dat, fit, regions.dict, 
#                              dat.constrain = FALSE, #if TRUE, only predict in the convex hull of the lat/lon of the data 
#                              do.pheno = TRUE, #if no doy term in the model, set to FALSE for faster calculations.
#                              sourcefac = "NFJ"){ #if FALSE, set up 2-color gradient moving outward from 0. 
#   #                                                  If TRUE, use a single color gradient adapted to observe GR.
#   ## using median year and median year + 1 to help minimize numerical weirdness.
#   ## 
#   state.map.data <- maps::map('state', fill = TRUE, plot = FALSE) %>%
#     st_as_sf()
#   
#   ## Prepping two sequential years for predicting
#   loc.plot = grid_plot_oneyr(dat, regions.dict, dat.constrain)
#   loc.plot$sourcefac = sourcefac
#   loc.plot.y1 = loc.plot.y2 = loc.plot
#   loc.plot.y1$year  = round(median(dat$year))
#   loc.plot.y2$year = loc.plot.y1$year + 1
#   
#   #predict
#   loc.plot.y1$count = predict(fit, newdata = loc.plot.y1, type = "response")
#   loc.plot.y2$count = predict(fit, newdata = loc.plot.y2, type = "response")
#   
#   ## combining
#   loc.plot$gr = log(loc.plot.y2$count)-log(loc.plot.y1$count)
#   
#   loc.plot = loc.plot %>% 
#     group_by(lon, lat) %>% 
#     summarize(gr.med = median(gr)) %>% 
#     ungroup()
#   
#   dat.pts = viz_filter(dat %>%  select(lon,lat, inferred), reso = 2)
#   gp = ggplot()+
#     geom_tile(data = loc.plot, 
#               aes(x = lon, y = lat, fill = gr.med))+
#     geom_point(data = dat.pts %>% filter(inferred == TRUE), aes(x = lon, y = lat), col = "coral", size = .2)+
#     geom_point(data = dat.pts %>% filter(inferred == FALSE), aes(x = lon, y = lat), col = "springgreen3", size = .2)+
#     geom_sf(data = state.map.data, fill = NA)+
#     theme_minimal()+
#     labs(fill = "Growth rate")+
#     ggtitle(paste0(dat$code[1], ": ", dat$sommon[1]," growth rates\ncoral points: inferred observations; green points: direct observations"))+
#     scale_fill_gradient2(low = muted("blue"),
#                          mid = "white",
#                          high = muted("green"),
#                          midpoint = 0)
#   return(list(fig = gp, data  = loc.plot))
# }
# 
# abund_mapper_old = function(dat, fit, regions.dict, sourcefac = "NFJ", 
#                             use.range = FALSE, #if TRUE, constrain to species range KML
#                             dat.constrain = FALSE,#if constraining geography to convex hull, set this to TRUE to only predict there.
#                             do.confidence = FALSE, #if true, calculate upper and lower confidence limits. Possibly. Interpretation seems tricky
#                             do.pheno = TRUE # set to FALSE if we're not fitting a phenology curve. if TRUE, predicts activity every day of the year. 
#                             # if FALSE, predicts for a single day in each year.
# ){
#   state.map.data <- maps::map('state', fill = TRUE, plot = FALSE) %>%
#     st_as_sf()
#   grid.plot = grid_plot_allyr(dat, regions.dict, use.range = use.range, 
#                               dat.constrain = dat.constrain, do.pheno = do.pheno)
#   grid.plot$sourcefac = sourcefac
#   ## confidence stuff is messy. Not keeping this updated.
#   if(do.confidence){
#     mod.pred = predict(fit, newdata = grid.plot, type = "response", se = TRUE)
#     grid.plot$count = mod.pred$fit
#     grid.plot$count.lower = mod.pred$fit - 1.96*mod.pred$se
#     grid.plot$count.upper = mod.pred$fit + 1.96*mod.pred$se
#     loc.sum = grid.plot %>% 
#       group_by(lat, lon) %>% 
#       summarize(abund.index = mean(count)*365,
#                 abund.lower = mean(count.lower)*365,
#                 abund.upper = mean(count.upper)*365)
#   }else{
#     grid.plot$count = predict(fit, newdata = grid.plot, type = "response")
#     grid.yearly = grid.plot %>% 
#       group_by(lat, lon, year) %>% 
#       summarize(abund.index = mean(count)*365) %>% 
#       ungroup()
#     loc.sum = grid.yearly %>% 
#       group_by(lat, lon) %>% 
#       summarize(abund.index = mean(abund.index)) %>% 
#       ungroup()
#   }
#   ## calculating abundance metrics
#   
#   ## yearly total abundance
#   abund.species = grid.yearly %>% 
#     group_by(year) %>% 
#     summarize(abund.index = sum(abund.index)) %>% 
#     ungroup()
#   ## yearly abundance at location of highest predicted density.
#   loc.highdense = loc.sum[which.max(loc.sum$abund.index),]
#   abund.highdense = grid.yearly %>% 
#     filter(lat == loc.highdense$lat, 
#            lon == loc.highdense$lon)
#   ## yearly abundance at location with most NFJ sightings. 
#   # reminder: predictions are happening at the rounded lat/lon, so let's find the rounded lat/lon with
#   # the most NFJ obs
#   dat.nfj = dat %>%
#     filter(source == "NFJ") %>% 
#     mutate(lon = round(lon),
#            lat = round(lat)) %>% 
#     group_by(lon, lat) %>% 
#     summarize(nfj.abund = mean(count)) %>% 
#     ungroup()
#   dat.nfj = inner_join(grid.yearly, dat.nfj)
#   abund.bestnfj = dat.nfj %>% 
#     filter(nfj.abund == max(nfj.abund))
#   ## visualizing
#   dat.pts = viz_filter(dat %>%  select(lon,lat, inferred), reso = 2)
#   gp = ggplot()+
#     geom_tile(data = loc.sum, aes(x = lon, y = lat, fill = log10(abund.index)))+
#     geom_point(data = dat.pts %>% filter(inferred == TRUE), aes(x = lon, y = lat), col = "coral", size = .2)+
#     geom_point(data = dat.pts %>% filter(inferred == FALSE), aes(x = lon, y = lat), col = "springgreen3", size = .2)+
#     geom_sf(data = state.map.data, fill = NA)+
#     scale_fill_viridis()+
#     theme_minimal()+
#     ggtitle(paste0(dat$code[1], ": ", dat$sommon[1], " estimated abundance (averaged across all years)\ncoral points: inferred observations; green points: direct observations"))
#   ## TROUBLE ADDING RANGE MAP. UPDATE HERE WITH ELIZA INPUT
#   # if(use.range){
#   #   range.map = rgdal::readOGR(here(paste0("2_data_wrangling/range-maps/", code.cur, ".kml")))
#   #   gp = gp + 
#   #     geom_path(data = range.map, aes(x = long, y = lat), fill = NA, color = "black")
#   # }
#   ## confidence stuff is messy. Not keeping this updated.
#   if(do.confidence){
#     gp.lower = ggplot()+
#       geom_tile(data = loc.sum, aes(x = lon, y = lat, fill = log10(abund.lower)))+
#       geom_point(data = dat.pts, aes(x = lon, y = lat), size = .2)+
#       geom_sf(data = state.map.data, fill = NA)+
#       scale_fill_viridis()+
#       theme_minimal()+
#       ggtitle(paste0(dat$code[1], " estimated abundance, LOWER CONFIDENCE LIMIT"))
#     gp.upper= ggplot()+
#       geom_tile(data = loc.sum, aes(x = lon, y = lat, fill = log10(abund.upper)))+
#       geom_point(data = dat.pts, aes(x = lon, y = lat), size = .2)+
#       geom_sf(data = state.map.data, fill = NA)+
#       scale_fill_viridis()+
#       theme_minimal()+
#       ggtitle(paste0(dat$code[1], " estimated abundance, UPPER CONFIDENCE LIMIT"))
#   }
#   ## not keeping do.confidence updated
#   if(do.confidence){
#     res = list(fig = gp, fig.lower = gp.lower, fig.upper = gp.upper, data  = loc.sum)
#   }else{
#     res = list(fig = gp, data = loc.sum, 
#                abund.species = abund.species, abund.highdense = abund.highdense, 
#                abund.bestnfj = abund.bestnfj)
#   }
#   return(res)
# }