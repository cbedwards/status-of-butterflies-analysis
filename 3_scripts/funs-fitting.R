library(sf)
library(tidyverse)
library(MASS, include.only = "glm.nb")
library(viridis)
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

grid_plot_allyr = function(dat, regions.dict, dat.constrain){
  ## for easy calculating abundance
  loc.pts = as.data.frame(expand.grid(lon = seq(min(dat$lon), max(dat$lon), by = 1),
                                      lat = seq(min(dat$lat), max(dat$lat), by = 1)))
  if(dat.constrain){
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
                          doy = seq(0,365, by = 1),
                          year = unique(dat$year))
  if("site.refac" %in% names(dat)){
    grid.plot$site.refac = as.factor("light on data")
  }
  return(grid.plot)
}

grid_plot_oneyr = function(dat, regions.dict, dat.constrain){
  ## for calculating trends
  loc.pts = as.data.frame(expand.grid(lon = seq(min(dat$lon), max(dat$lon), by = .2),
                                      lat = seq(min(dat$lat), max(dat$lat), by = .2)))
  if(dat.constrain){
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
                          doy = seq(0,365, by = 1))
  if("site.refac" %in% names(dat)){
    grid.plot$site.refac = as.factor("light on data")
  }
  return(grid.plot)
}

abund_mapper = function(dat, fit, regions.dict, sourcefac = "NFJ", 
                        dat.constrain = FALSE,#if constraining geography to convex hull, set this to TRUE to only predict there.
                        do.confidence = FALSE){ #if true, calculate upper and lower confidence limits. Possibly. Interpretation seems tricky
  state.map.data <- maps::map('state', fill = TRUE, plot = FALSE) %>%
    st_as_sf()
  grid.plot = grid_plot_allyr(dat, regions.dict, dat.constrain)
  grid.plot$sourcefac = sourcefac
  if(do.confidence){
    mod.pred = predict(fit, newdata = grid.plot, type = "response", se = TRUE)
    grid.plot$count = mod.pred$fit
    grid.plot$count.lower = mod.pred$fit - 1.96*mod.pred$se
    grid.plot$count.upper = mod.pred$fit + 1.96*mod.pred$se
    loc.sum = grid.plot %>% 
      group_by(lat, lon) %>% 
      summarize(abund.index = mean(count)*365,
                abund.lower = mean(count.lower)*365,
                abund.upper = mean(count.upper)*365)
  }else{
    grid.plot$count = predict(fit, newdata = grid.plot, type = "response")
    loc.sum = grid.plot %>% 
      group_by(lat, lon) %>% 
      summarize(abund.index = mean(count)*365)
  }
  gp = ggplot()+
    geom_tile(data = loc.sum, aes(x = lon, y = lat, fill = log10(abund.index)))+
    geom_point(data = dat %>% filter(inferred == TRUE), aes(x = lon, y = lat), col = "coral", size = .2)+
    geom_point(data = dat %>% filter(inferred == FALSE), aes(x = lon, y = lat), col = "springgreen3", size = .2)+
    geom_sf(data = state.map.data, fill = NA)+
    scale_fill_viridis()+
    theme_minimal()+
    ggtitle(paste0(dat$code[1], ": ", dat$sommon[1], " estimated abundance (averaged across all years)\ncoral points: inferred observations; green points: direct observations"))
  if(do.confidence){
    gp.lower = ggplot()+
      geom_tile(data = loc.sum, aes(x = lon, y = lat, fill = log10(abund.lower)))+
      geom_point(data = dat, aes(x = lon, y = lat), size = .2)+
      geom_sf(data = state.map.data, fill = NA)+
      scale_fill_viridis()+
      theme_minimal()+
      ggtitle(paste0(dat$code[1], " estimated abundance, LOWER CONFIDENCE LIMIT"))
    gp.upper= ggplot()+
      geom_tile(data = loc.sum, aes(x = lon, y = lat, fill = log10(abund.upper)))+
      geom_point(data = dat, aes(x = lon, y = lat), size = .2)+
      geom_sf(data = state.map.data, fill = NA)+
      scale_fill_viridis()+
      theme_minimal()+
      ggtitle(paste0(dat$code[1], " estimated abundance, UPPER CONFIDENCE LIMIT"))
  }
  if(do.confidence){
    res = list(fig = gp, fig.lower = gp.lower, fig.upper = gp.upper, data  = loc.sum)
  }else{
    res = list(fig = gp, data = loc.sum)
  }
  return(res)
}

trend_plotter = function(dat, fit, regions.dict, 
                         dat.constrain = FALSE, #if TRUE, only predict in the convex hull of the lat/lon of the data 
                         sourcefac = "NFJ"){ #if FALSE, set up 2-color gradient moving outward from 0. 
  #                                                  If TRUE, use a single color gradient adapted to observe GR.
  ## using median year and median year + 1 to help minimize numerical weirdness.
  ## 
  state.map.data <- maps::map('state', fill = TRUE, plot = FALSE) %>%
    st_as_sf()
  
  ## Prepping two sequential years for predicting
  loc.plot = grid_plot_oneyr(dat, regions.dict, dat.constrain)
  loc.plot$sourcefac = sourcefac
  loc.plot.y1 = loc.plot.y2 = loc.plot
  loc.plot.y1$year  = round(median(dat$year))
  loc.plot.y2$year = loc.plot.y1$year + 1
  
  #predict
  loc.plot.y1$count = predict(fit, newdata = loc.plot.y1, type = "response")
  loc.plot.y2$count = predict(fit, newdata = loc.plot.y2, type = "response")
  
  ## combining
  loc.plot$gr = log(loc.plot.y2$count)-log(loc.plot.y1$count)
  
  loc.plot = loc.plot %>% 
    group_by(lon, lat) %>% 
    summarize(gr.med = median(gr)) %>% 
    ungroup()
  
  gp = ggplot()+
    geom_tile(data = loc.plot, 
              aes(x = lon, y = lat, fill = gr.med))+
    geom_point(data = dat %>% filter(inferred == TRUE), aes(x = lon, y = lat), col = "coral", size = .2)+
    geom_point(data = dat %>% filter(inferred == FALSE), aes(x = lon, y = lat), col = "springgreen3", size = .2)+
    geom_sf(data = state.map.data, fill = NA)+
    theme_minimal()+
    labs(fill = "Growth rate")+
    ggtitle(paste0(dat$code[1], ": ", dat$sommon[1]," growth rates\ncoral points: inferred observations; green points: direct observations"))+
    scale_fill_gradient2(low = muted("blue"),
                         mid = "white",
                         high = muted("green"),
                         midpoint = 0)
  return(list(fig = gp, data  = loc.plot))
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
  ggplot(dat.pred, aes(x = doy, y = count)) +
    geom_point(data = dat.plot, aes(col = sourcefac))+
    geom_line(col = 'black') +
    facet_wrap( ~ year)+
    ggtitle(paste0(dat$code[1],": ", dat$sommon[1], " at (", lon.plot, ", ", lat.plot,") predicted for ", names(rev(sort(table(dat.plot$sourcefac))))[1]))+
    theme_minimal()
}

NFJ_compare = function(dat, fit, regions.dict, nyears = 5,
                       across.doy = TRUE,#if FALSE, will just look at points estimates of doy for predictions. If TRUE, calculates activity-days in that year
                       across.year = TRUE #if FALSE, will treat each year separately. If TRUE, will report the average across all years. TRUE is probably more
                       #faithful as a model validation, as the model isn't structured to capture year-to-year fluctations.
){ 
  dat.comp = dat %>% 
    filter(sourcefac == "NFJ") %>% 
    filter(max(year)-year < nyears)
  if(across.doy){
    dat.pred = expand_grid(nesting(dat.comp[, c("lat", "lon", "year", "regionfac")]),
                           doy = 0:365)
    if("site.refac" %in% names(dat)){
      dat.pred$site.refac = as.factor("light on data")
    }
    dat.pred$sourcefac = "NFJ"
    dat.pred$count.pred = predict(fit, newdata = dat.pred, type = "response")
    dat.pred = dat.pred %>% 
      group_by(lat, lon, year) %>% 
      summarize(count.pred = mean(count.pred)*365) %>% 
      ungroup()
    dat.full = inner_join(dat.comp, dat.pred)
  }else{
    dat.full = dat.comp
    dat.full$count.pred = predict(fit, newdata = dat.full, type = "response")
  }
  if(across.year){
    dat.full = dat.full %>% 
      group_by(lat, lon, code, region) %>% 
      summarize(count = mean(count),
                count.pred = mean(count.pred)) %>% 
      ungroup()
  }
  dat.cor = cor.test(dat.full$count, dat.full$count.pred)
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
  gp = ggplot(dat.full, aes(x = count, y = count.pred))+
    geom_point(aes(col = region))+
    geom_smooth(method = lm)+
    # geom_smooth(method = lm, formula = dat.full$count.pred ~ dat.full$count:dat.full$region)+
    ggtitle(paste0(dat.full$code[1],": ", dat$sommon[1], ", NFJ counts vs predictions, last ", nyears, " years\n",
                   "Correlation of ", round(dat.cor$estimate,4),", ", p_pretty(dat.cor$p.value)))+
    xlab("Actual count")+
    ylab("Prediction")
  return(list(fig = gp, cor = round(dat.cor$estimate,4)))
}

NFJ_regional_trends = function(dat, regions.dict){
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



