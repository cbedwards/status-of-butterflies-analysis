## Code to plot phenology curves from fitted model_runner obect
## 
## Assuming you've already run something along the lines of 
## out <- model_runner(STUFF HERE)

## To look at individual source or region (also, for cases when doy doesn't 
## vary across those factors)
cur.region = "Midwest"
cur.source = "NFJ"
## if we're allowing our DOY smooth to vary between source or region,
## change cur.region or cur.source to span all options
## example:
# cur.region = unique(out$data$regionfac)
lat.cur = median(out$data$lat)
lon.cur = median(out$data$lon)
dat.pred = as.data.frame(expand.grid(year = 2010,
                                     doy = seq(0, 365, by = .1),
                                     lat = lat.cur,
                                     lon = lon.cur,
                                     regionfac = cur.region,
                                     sourcefac = cur.source,
                                     effort.universal = 0,
                                     effort.universal.type = "duration",
                                     site.refac = "totally new site"))
dat.pred$count = predict(out$fitted.model, newdata = dat.pred, type = 'response',
                         discrete = FALSE)
ggplot(dat.pred, aes(x = doy, y = count))+
  geom_path()+
  facet_wrap(. ~ sourcefac + regionfac)
