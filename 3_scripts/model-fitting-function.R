## Function for fitting a single model and returning various results
library(tidyverse)
library(formula.tools)
library(rmarkdown)
library(sp)
library(here)
source(here("3_scripts/funs.R"))
source(here("3_scripts/funs-fitting.R"))

model_runner = function(code.cur, #GU code for taxa of interest
                        form.use, #gam formula to fit
                        knot_maker, #function taking as sole argument `dat`, handles custom knot placement. 
                        #            for default knot placement, knot_maker = function(dat){return(list())}
                        sitere_maker, # function to create site.refac, factors corresponding to sites for random effects.
                        #               if not using site-level random effects in your model, 
                        #               sitere_maker = function(dat){return(NA)}
                        #               if using each site as a separate level for random effects,
                        #               sitere_maker = function(dat){return(as.factor(dat$site))}
                        regions.dict, #dictionary to map states to regions
                        fit.family = "nb",
                        use.range = FALSE, #if TRUE, reads in associated range map from 2_data_wrangling and cuts observations outside it
                        suitability.cutoff = NA, #if NA, use kml polygons of range map. If specified, estimated habitat suitability from .tiff, and include points with the specified suitability
                        use.inferred = TRUE, #if TRUE, infer zeros from community surveys that didn't report focal species
                        infer.messy.levels = c("GENUS", "SUBFAMILY", "FAMILY", "COMPLEX"), #for infering zeroes, what level of "unidentified" do we not infer zeroes
                        ## at the default, if a trip did not report the focal species, but did report an unknown in the family, subfamily, genus, or species complex
                        ## we do NOT infer zeroes. FOr a more selective approach, try setting to just c("GENUS", "COMPLEX") - this will infer more zeroes.
                        do.pheno = TRUE, #If model doesn't include doy term, set to FALSE to speed up calculations. MAKE SURE THIS IS TRUE WHEN THERE IS A DOY TERM IN THE MODEL!
                        geography.constrain = FALSE, #if TRUE, restrict data to only the observations within the convex hull (in lat/lon)
                        #                               of non-inferred data. NOTE: this is overridden if using the range maps (if use.range == TRUE)
                        pheno.window = NULL, #if NULL, include all data. Otherwise, specify vector of probabilities for quantiles, and data will be clipped to the doy defined by those quanitles.
                        #we recommend c(0.001, 0.999) as a good window choice.
                        regions.use = NULL, #if specified, limit analyses to only the provided FWS region(s)
                        min.year = -9999, #only use years after this point. For fast check of Wayne's question
                        use.only.source= NULL, #if NULL, use all sources. If a vector of characters, use only the specified sources.
                        n.threads.use = 4){ #how many threads to let mgcv::bam() use. 
  #                                          Decrease if you're running into computer performance issues.
  
  if(!do.pheno){cat("do.pheno set to FALSE, calculations of trends + abundance calculated accordingly.\nCheck that there is no `doy` term in the model!\n")}
  if(use.inferred){
    events.missed.messy = make_dataset(code = code.cur, 
                 use.range = use.range, 
                 infer.messy.levels = infer.messy.levels)
    dat = qread(paste0("2_data_wrangling/cleaned by code/",code.cur, ".csv")) 
  }else{
    ## curious if the inferred 0s are messing things up
    dat = qread(paste0("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv"))
    dat = dat %>%
      filter(code == code.cur)
    events.missed.messy = NA
  }
  dat = dat %>% 
    filter(!is.na(count))
  
  ## filter by year:
  dat = dat[dat$year >= min.year,]
  
  ## filter by doy
  if(!is.null(pheno.window)){
    dat.nzero = dat %>% 
      filter(count > 0)
    doy.window = quantile(dat.nzero$doy, pheno.window, na.rm = T)
    cat("Clipping data to the the window between these days:\n")
    print(doy.window)
    dat = dat[dat$doy >= doy.window[[1]],]
    dat = dat[dat$doy <= doy.window[[2]],]
  }
  

  
  ## 
  ## UPDATE: this is being applied in the data generation step (make_dataset) to reduce size of data files.
  ## Crop to range, if that option is taken
  ## 
  # if(use.range){
  #   dat.sp = dat # separate object for spatial coordinates, because I'm paranoid about interactions with other code #use polygon range map
  #   range.map = rgdal::readOGR(here(paste0("2_data_wrangling/range-maps/", code.cur, ".kml")))
  #   coordinates(dat.sp) <- ~ lon + lat
  #   sp::proj4string(dat.sp) <- sp::proj4string(range.map)
  #   if(is.na(suitability.cutoff)){ ## just use polygon of range map
  #     overlaid <- sp::over(dat.sp, as(range.map, "SpatialPolygons"))
  #     ind.within <- which(!is.na(overlaid))
  #   }else{ #use gradient of habitat suitability and specified cutoff.
  #     hsm <- raster::raster(here(paste0("2_data_wrangling/range-maps/", code.cur, ".tiff")))
  #     site.suitability <- raster::extract(hsm, dat.sp)
  #     ind.within <- which(site.suitability > suitability.cutoff)
  #   }
  #   dat = dat[ind.within,]
  # }
  
  ## make factor version of source
  dat$sourcefac = as.factor(dat$source)
  ## make factor version of universal effort type
  dat$effort.universal.type = as.factor(dat$effort.universal.type)
  
  ## adding customized site levels for random effects of sites
  dat$site.refac = sitere_maker(dat)
  
  ## add regions and a factor version of regions to the data  
  dat = left_join(dat, regions.dict, by = "state")
  
  if(!is.null(regions.use)){
    dat = dat %>% 
      filter(region %in% regions.use)
  }
  

  
  ## FOR TESTING PURPOSES let's remove region = NA -- there are a few sites that are 
  ## at boundaries of coastlines and are being mistakenly identified as state == ""
  dat = dat %>% 
    filter(!is.na(dat$region))
  dat$regionfac = as.factor(dat$region)
  
  cat(paste0("working with a total of ", nrow(dat), " records.\n"))
  
  ##round counts to nearest whole number - some data are measured using distance sampling
  ##or other methods that can give decimals.
  dat$count = round(dat$count)
  
  # make simple data/data coverage plots
  gp.counts = ggplot(dat) + 
    geom_point(aes(x = doy, y = count, col = inferred), shape = 1)+
    facet_wrap(~regionfac)+ggtitle("Counts across regions")+
    xlab("day of year")+
    theme_minimal()
  gp.hist = ggplot(dat) +
    geom_histogram(aes(x=doy, fill = count>0)) +
    facet_wrap(~ inferred)+
    ggtitle("Histogram of # records across doy\nseparated by whether inferred or not (TRUE = inferred)")+
    xlab("day of year")+
    theme_minimal()
  
  ## Get species names
  names.vec = unique(dat$name)
  ##remove "absence for" entries in list of species names -- these are to denote inferences
  names.vec = names.vec[!grepl("absence for ", names.vec)]
  plot.title = paste0(dat$common[1], " (", code.cur,")", ": " , paste0(names.vec, collapse = ", "))
  cat(paste0("\n******************************************\nFitting ", plot.title,
             "\n******************************************\n"))
  
  ## constraining geography if called for - this version is clipping to the convex hull
  ## My thinking here is that splines can misbehave if we have many observations (of 0)
  ## in biologically uninteresting areas (ie lots of 0s on the east coast for a 
  ## west-coast-only species). This can happen because of our conferred zeroes
  ## 
  ## Constraining data might be best with range maps, but a quick-and-dirty solution
  ## is to only use observations within the convex hull of our non-inferred data.
  ## 
  ## Potential elaboration: add a buffer of ~X lat/lon around each non-inferred data point.
  ## 
  if(geography.constrain == TRUE){
    shape.hull = convhulln(dat[dat$inferred==FALSE,c("lon", "lat")])
    dat = dat[inhulln(shape.hull, as.matrix(dat[, c("lon","lat")])),]
  }
  
  ## Filtering source, if any
  if(!is.null(use.only.source)){
    dat = dat %>% 
      filter(source %in% use.only.source)
  }
  
  #
  cat("fitting gam\n")
  time.start = proc.time()
  fit = bam(form.use,
            data = dat,
            method="fREML", 
            family=fit.family,
            knots=knot_maker(),
            discrete = TRUE,
            nthreads = n.threads.use)
  time.end = proc.time()
  print(time.end-time.start)
  print(fit$formula)
  
  if(FALSE){
    gam.check(fit)
  }
  
  ## generate plots --------------
  ## plot abundance map
  cat("calculating abundance and trends\n")
  
  out.abund.and.trends = trend_and_abund_calc(dat=dat, fit = fit, 
                                              regions.dict = regions.dict, 
                                              fit.family = fit.family,
                                              use.range = use.range, 
                                              dat.constrain = geography.constrain, 
                                              do.pheno = do.pheno)
  loc.plot = out.abund.and.trends$loc.plot ## for simplicity
  cat("making heatmap plots\n")
  #plot abundance
  out.abund = abund_plotter(dat = dat, loc.plot = loc.plot)
  ## plot trends
  out.trend = trend_plotter(dat = dat, loc.plot = loc.plot)
  cat("identifying trends from key locations\n")
  trend.species = out.abund.and.trends$global.gr
  
  trend.highabund = loc.plot$gr.med[which.max(loc.plot$abund.index)]
  
  
  
  if(sum(dat$source == "NFJ")>0){
    dat.nfj = dat %>%
      filter(source == "NFJ") %>% 
      mutate(lon = round(lon),
             lat = round(lat)) %>% 
      group_by(lon, lat) %>% 
      summarize(nfj.abund = mean(count)) %>% 
      ungroup()
    loc.plot.merge = loc.plot %>% 
      mutate(lon = round(lon),
             lat = round(lat)) %>% 
      group_by(lon, lat) %>% 
      summarize(gr.med = mean(gr.med),
                abund.index = mean(abund.index)) %>% 
      ungroup()
    dat.nfj = inner_join(loc.plot.merge, dat.nfj, by = c("lon", "lat"))
    dat.nfj = dat.nfj %>% 
      filter(nfj.abund == max(nfj.abund))
    trend.mostnfj = dat.nfj$gr.med[1];
  }
  else{
    trend.mostnfj = NA
  }
  
  ## Plot NFJ abundance by region
  ## IF NFJ is in the data
  if(sum(dat$source == "NFJ")>0){
    cat("Comparing to NFJ abundance\n")
    compare.abund = NFJ_compare(dat, fit, regions.dict, fit.family = fit.family, use.range = use.range, nyears = 10, across.doy = FALSE)
  }else{
    cat("No NFJ data, so no comparison to NFJ abundance")
    compare.abund = list(fig = ggplot()+ggtitle("No NFJ data"), cor = NA)
  }
  
  ## Plot NFJ trends by region
  ## This is no longer being used in the report. uncomment if being used again.
  # cat("Comparing to NFJ trends")
  # gp.nfj.trend = NFJ_regional_trends(dat, regions.dict, use.range = use.range)
  
  ## calculate some diagnostics for info-plot
  
  ## Turn knots.list into something human-readable. 
  knots.vec = "knots: "
  knots.list = knot_maker(dat)
  if(length(knots.list)==0){
    knots.vec = "knots: all automatically placed"
  }else{
    for(cur.name in names(knots.list)){
      knots.vec = paste0(knots.vec,
                         paste0(cur.name,": ", paste0(round(knots.list[[cur.name]],1), collapse = ", "), "\n       "))
    }
  }
  ## Page with summary information (in words)
  # plot.form = as.character(formula(fit))[3]
  # plot.form = gsub("[+]", "+\n  ", plot.form)
  diagnostics.text = c(heading = "Model fitting summary:",
                       title = plot.title,
                       formula = as.character(form.use),
                       sources = paste0("Sources limited to:  ", paste0(use.only.source, collapse=", ")),
                       knots = knots.vec
  )
  if(is.null(use.only.source)){diagnostics.text$sources = "All sources used"}
  cat("\n")
  return(list(plot.title = plot.title,
              loc.plot = loc.plot,
              abund.cor = compare.abund$cor,
              fig.abund = out.abund$fig,
              fig.trend = out.trend$fig,
              fig.nfj.abund = compare.abund$fig,
              fig.hist =gp.hist,
              fig.counts = gp.counts,
              diagnostics.text = diagnostics.text,
              fitted.model = fit,
              trend.species = trend.species,
              trend.highabund = trend.highabund,
              trend.mostnfj = trend.mostnfj,
              n = nrow(dat),
              data = dat,
              events.missed.messy = events.missed.messy))
}

report_maker = function(modelfit.output,
                        code.cur,
                        res.path, 
                        copy.heatmaps = FALSE,
                        fig.width = 7.5, 
                        fig.height = 5){
  run.name = gsub(".*/", "", res.path)
  filename.use = paste0(code.cur,
                        "-", run.name,".html")
  path.use = here(paste0(res.path ,"/", filename.use))
  
  ##fitting individual models is only taking a few minutes at most.
  ##Probably easier to re-fit than to save and load
  # saveRDS(modelfit.output$fitted.model,
  #         here(paste0(res.path, "/", code.cur, "-fitted-model.RDS")))
  ## saving jpgs for making the report
  dir.create(here("4_res/temp-files"))
  ggsave(here("4_res/temp-files/cur.abund.jpg"),
         modelfit.output$fig.abund, 
         width = fig.width, height = fig.height)
  ggsave(here("4_res/temp-files/cur.nfj.abund.jpg"),
         modelfit.output$fig.nfj.abund, 
         width = fig.width, height = fig.height)
  ggsave(here("4_res/temp-files/cur.trend.jpg"),
         modelfit.output$fig.trend, 
         width = fig.width, height = fig.height)
  ## copy heatmaps for later use in Rshiny app. Slows process
  if(copy.heatmaps){
    file.copy(from = here("4_res/temp-files/cur.abund.jpg"),
              to = here(paste0(res.path,"/heatmaps/", code.cur, "-abund-map.jpg")))
    file.copy(from = here("4_res/temp-files/cur.trend.jpg"),
              to = here(paste0(res.path,"/heatmaps/", code.cur, "-trend-map.jpg")))
  }
  # ggsave(here("4_res/fit-summaries/temp-files/cur.nfj.trend.jpg"),
  #        modelfit.output$fig.nfj.trend, 
  #        width = fig.width, height = fig.height)
  saveRDS(modelfit.output$diagnostics.text,
          here("4_res/temp-files/cur.diagnostics.RDS"))
  ## Note: I've cut the following from the reports for the sake of speed. 
  ## These plots didn't seem particularly useful, and plot-saving + report-making
  ## is ~ 1/2 of the total runtime. If these plots are returned in the report,
  ## uncomment these lines.
  # ggsave(here("4_res/temp-files/cur.counts.jpg"),
  #        modelfit.output$fig.counts, 
  #        width = fig.width, height = fig.height)
  # ggsave(here("4_res/ftemp-files/cur.hist.jpg"),
  #        modelfit.output$fig.hist, 
  #        width = fig.width, height = fig.height)
  # ggsave(here("4_res/temp-files/cur.activity.maxdata.jpg"),
  #        modelfit.output$fig.activity.maxdata, 
  #        width = fig.width, height = fig.height)
  # ggsave(here("4_res/temp-files/cur.activity.maxabund.jpg"),
  #        modelfit.output$fig.activity.maxabund, 
  #        width = fig.width, height = fig.height)
  render(input = here("3_scripts/species-fitting-report.Rmd"),
         output_dir = here(res.path),
         params = list(title = out$diagnostics.text[["title"]]),
         output_file = filename.use)
  return(path.use)
  files.temp = list.files(here("4_res/temp-files"))
  file.remove(paste(here("4_res/temp-files/", files.temp, collapse = "")))
}


# ## Note: this makes a big pdf that's hard to load/view.
# model_saver = function(modelfit.output,
#                        code.cur,
#                        run.suffix){
#   ## identify next available version number (to avoid overwriting)
#   cur.files = list.files(here(paste0("4_res/fit-summaries/")))
#   cur.file.code = cur.files[grepl(paste0(code.cur,"-", run.suffix), cur.files)]
#   cur.file.code = cur.file.code[grepl("[.]pdf", cur.file.code)]
#   if(length(cur.file.code)>0){
#     cur.file.code = gsub("[.]pdf", "", cur.file.code)
#     cur.nums = as.numeric(gsub(paste0(code.cur, "-", run.suffix,"-V"), "", cur.file.code))
#     use.num = max(cur.nums)+1
#   }else{
#     use.num = 1
#   }
#   ## saving PDF
#   path.use = here(paste0("4_res/fit-summaries/",
#                          code.cur,
#                          "-", run.suffix,
#                          "-V", use.num, ".pdf"))
#   cat(paste0("Saving ", modelfit.output$plot.title," to\n",
#              path.use))
#   pdf(path.use,
#       width = 15, height = 20)
#   ## map of abundance
#   gp = ggarrange(modelfit.output$fig.abund,
#                  modelfit.output$fig.njf.abund, ncol = 1)
#   print(annotate_figure(gp, top = text_grob(modelfit.output$plot.title, size = 18)))
#   ## map of trends
#   gp = ggarrange(modelfit.output$fig.trend, 
#                  modelfit.output$fig.nfj.trend, ncol = 1)
#   print(annotate_figure(gp, top = text_grob(modelfit.output$plot.title, size = 18)))
#   
#   plot(NA, xlim = c(0,5), ylim = c(0,5), bty = "n", 
#        xaxt = "n", yaxt = "n", xlab = "", ylab = "")
#   text(1,4, modelfit.output$diagnostics.text["heading"], pos = 4, cex = 3)
#   text(1,3, modelfit.output$diagnostics.text["title"], pos = 4, cex = 2)
#   text(1,2, modelfit.output$diagnostics.text["formula"], pos = 4, cex = 2)
#   text(1,1.5, modelfit.output$diagnostics.text["sources"], pos = 4)
#   text(1,1, modelfit.output$diagnostics.text["knots"], pos = 4)
#   
#   gp = ggarrange(modelfit.output$fig.hist, modelfit.output$fig.counts, ncol=1)
#   print(annotate_figure(gp, top = text_grob(modelfit.output$plot.title, size = 18)))
#   
#   ## activity curves for highest data denstiy, highest estimated abundance
#   print(modelfit.output$fig.activity.maxdata + ggtitle(
#     paste0(modelfit.output$fig.activity.maxdata$labels$title, "\n(Location of highest data density)")
#   ))
#   print(modelfit.output$fig.activity.maxabund + ggtitle(
#     paste0(modelfit.output$fig.activity.maxabund$labels$title, "\n(Location of highest estimated abundance)")
#   ))
#   dev.off()
#   ## save jpgs for faster viewing of maps etc.
#   ggsave(here(paste0("4_res/fit-summaries/jpgs-fastviews/",
#                      code.cur, "-", run.suffix,
#                      "-V", use.num, "-abund-map.jpg")),
#          modelfit.output$fig.abund, 
#          width = 15, height = 12)
#   ggsave(here(paste0("4_res/fit-summaries/jpgs-fastviews/",
#                      code.cur, "-", run.suffix,
#                      "-V", use.num, "-trend-map.jpg")),
#          modelfit.output$fig.trend, 
#          width = 15, height = 12)
#   ggsave(here(paste0("4_res/fit-summaries/jpgs-fastviews/",
#                      code.cur, "-", run.suffix,
#                      "-V", use.num, "-trend-NFJ.jpg")),
#          modelfit.output$fig.nfj.abund, 
#          width = 15, height = 12)
#   cat("\n")
#   return(path.use)
# }
# 
