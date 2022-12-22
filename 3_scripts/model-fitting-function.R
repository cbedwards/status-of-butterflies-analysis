## Function for fitting a single model and returning various results
library(tidyverse)
library(formula.tools)
library(rmarkdown)

model_runner = function(code.cur, #GU code for taxa of interest
                        form.use, #gam formula to fit
                        knot_maker, #function taking as sole argument `dat`, handlescustom knot placement. 
                        #            for default knot placement, knot_maker = function(dat){return(list())}
                        sitere_maker, # function to create site.refac, factors corresponding to sites for random effects.
                        #               if not using site-level random effects in your model, 
                        #               sitere_maker = function(dat){return(NA)}
                        #               if using each site as a separate level for random effects,
                        #               sitere_maker = function(dat){return(as.factor(dat$site))}
                        regions.dict, #dictionary to map states to regions
                        use.inferred = TRUE, #if TRUE, infer zeros from community surveys that didn't report focal species
                        geography.constrain = FALSE, #if TRUE, restrict data to only the observations within the convex hull (in lat/lon)
                        #                               of non-inferred data
                        use.only.source= NULL, #if NULL, use all sources. If a vector of characters, use only the specified sources.
                        n.threads.use = 4){ #how many threads to let mgcv::bam() use. 
  #                                          Decrease if you're running into computer performance issues.
  if(use.inferred){
    make_dataset(code = code.cur)
    dat = qread(paste0("2_data_wrangling/cleaned by code/",code.cur, ".csv"))
  }else{
    ## curious if the inferred 0s are messing things up
    dat = qread(paste0("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv"))
    dat = dat %>%
      filter(code == code.cur)
  }
  dat = dat %>% 
    filter(!is.na(count))
  # dim(dat)
  # head(dat)
  
  ## make factor version of source
  dat$sourcefac = as.factor(dat$source)
  
  ## adding customized site levels for random effects of sites
  dat$site.refac = sitere_maker(dat)
  
  ## add regions and a factor version of regions to the data  
  dat = left_join(dat, regions.dict, by = "state")
  ## FOR TESTING PURPOSES let's remove region = NA -- there are a few sites that are 
  ## at boundaries of coastlines and are being mistakenly identified as state == ""
  dat = dat %>% 
    filter(!is.na(dat$region))
  dat$regionfac = as.factor(dat$region)
  
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
  
  ## constraining geography if called for
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
  
  print(paste0("Fitting ", plot.title))
  time.start = proc.time()
  fit = bam(form.use,
            data = dat,
            method="fREML", 
            family="nb",
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
  print("calculating abundance")
  out.abund = abund_mapper(dat, fit, regions.dict, dat.constrain = geography.constrain)
  out.abund$fig
  
  print("calculating trends")
  ## plot trends
  out.trend = trend_plotter(dat, fit, regions.dict, 
                            dat.constrain = geography.constrain)
  out.trend$fig #+ scale_fill_viridis()
  
  cat("calculating activity curve at data-dense region")
  ## plot activity curves for point of max data
  gp.activity.maxdata = demo_activity_plots(dat, fit, regions.dict)
  ## plot activity curves for point of max estimated density (diangostic for unreasonable activity curves)
  pt.maxabund = out.abund$data[which.max(out.abund$data$abund.index),]
  cat("calculating activity curve at high estimated abundance region")
  gp.activity.maxabund = activity_plotter(dat, fit, regions.dict, 
                                          lat.plot = pt.maxabund$lat, 
                                          lon.plot = pt.maxabund$lon,
                                          allyears = TRUE,
                                          source.adaptive = FALSE)
  
  ## Plot NFJ abundance by region
  cat("Comparing to NFJ abundance")
  compare.abund = NFJ_compare(dat, fit, regions.dict, nyears = 10)
  
  ## Plot NFJ trends by region
  cat("Comparing to NFJ trends")
  gp.nfj.trend = NFJ_regional_trends(dat, regions.dict)
  
  ## calculate some diagnostics for info-plot
  
  ## Turn knots.list into something human-readable. 
  knots.vec = "knots: "
  knots.list = knot_maker()
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
  cat("\n")
  return(list(plot.title = plot.title,
              dat.abund = out.abund$data,
              dat.trend = out.trend$data,
              abund.cor = compare.abund$cor,
              fig.abund = out.abund$fig,
              fig.nfj.abund = compare.abund$fig,
              fig.trend = out.trend$fig,
              fig.nfj.trend = gp.nfj.trend,
              fig.activity.maxdata = gp.activity.maxdata,
              fig.activity.maxabund = gp.activity.maxabund,
              fig.hist =gp.hist,
              fig.counts = gp.counts,
              diagnostics.text = diagnostics.text))
}

report_maker = function(modelfit.output,
                        code.cur,
                        run.suffix, 
                        fig.width = 7.5, 
                        fig.height = 5){
  cur.files = list.files(here(paste0("4_res/fit-summaries/")))
  cur.file.code = cur.files[grepl(paste0(code.cur,"-", run.suffix), cur.files)]
  cur.file.code = cur.file.code[grepl("[.]pdf", cur.file.code)]
  if(length(cur.file.code)>0){
    cur.file.code = gsub("[.]pdf", "", cur.file.code)
    cur.nums = as.numeric(gsub(paste0(code.cur, "-", run.suffix,"-V"), "", cur.file.code))
    use.num = max(cur.nums)+1
  }else{
    use.num = 1
  }
  filename.use = paste0(code.cur,
                        "-", run.suffix,
                        "-V", use.num, ".html")
  path.use = here(paste0("4_res/fit-summaries/", filename.use))
  
  ## saving jpgs for making the report
  ggsave(here("4_res/fit-summaries/temp-files/cur.abund.jpg"),
         modelfit.output$fig.abund, 
         width = fig.width, height = fig.height)
  ggsave(here("4_res/fit-summaries/temp-files/cur.nfj.abund.jpg"),
         modelfit.output$fig.nfj.abund, 
         width = fig.width, height = fig.height)
  ggsave(here("4_res/fit-summaries/temp-files/cur.trend.jpg"),
         modelfit.output$fig.trend, 
         width = fig.width, height = fig.height)
  ggsave(here("4_res/fit-summaries/temp-files/cur.nfj.trend.jpg"),
         modelfit.output$fig.nfj.trend, 
         width = fig.width, height = fig.height)
  saveRDS(modelfit.output$diagnostics.text,
          here("4_res/fit-summaries/temp-files/cur.diagnostics.RDS"))
  ggsave(here("4_res/fit-summaries/temp-files/cur.counts.jpg"),
         modelfit.output$fig.counts, 
         width = fig.width, height = fig.height)
  ggsave(here("4_res/fit-summaries/temp-files/cur.hist.jpg"),
         modelfit.output$fig.hist, 
         width = fig.width, height = fig.height)
  ggsave(here("4_res/fit-summaries/temp-files/cur.activity.maxdata.jpg"),
         modelfit.output$fig.activity.maxdata, 
         width = fig.width, height = fig.height)
  ggsave(here("4_res/fit-summaries/temp-files/cur.activity.maxabund.jpg"),
         modelfit.output$fig.activity.maxabund, 
         width = fig.width, height = fig.height)
  render(input = here("3_scripts/species-fitting-report.Rmd"),
         output_dir = here("4_res/fit-summaries/"),
         params = list(title = out$diagnostics.text[["title"]]),
         output_file = filename.use)
  return(path.use)
}


## Note: this makes a big pdf that's hard to load/view.
model_saver = function(modelfit.output,
                       code.cur,
                       run.suffix){
  ## identify next available version number (to avoid overwriting)
  cur.files = list.files(here(paste0("4_res/fit-summaries/")))
  cur.file.code = cur.files[grepl(paste0(code.cur,"-", run.suffix), cur.files)]
  cur.file.code = cur.file.code[grepl("[.]pdf", cur.file.code)]
  if(length(cur.file.code)>0){
    cur.file.code = gsub("[.]pdf", "", cur.file.code)
    cur.nums = as.numeric(gsub(paste0(code.cur, "-", run.suffix,"-V"), "", cur.file.code))
    use.num = max(cur.nums)+1
  }else{
    use.num = 1
  }
  ## saving PDF
  path.use = here(paste0("4_res/fit-summaries/",
                         code.cur,
                         "-", run.suffix,
                         "-V", use.num, ".pdf"))
  cat(paste0("Saving ", modelfit.output$plot.title," to\n",
             path.use))
  pdf(path.use,
      width = 15, height = 20)
  ## map of abundance
  gp = ggarrange(modelfit.output$fig.abund,
                 modelfit.output$fig.njf.abund, ncol = 1)
  print(annotate_figure(gp, top = text_grob(modelfit.output$plot.title, size = 18)))
  ## map of trends
  gp = ggarrange(modelfit.output$fig.trend, 
                 modelfit.output$fig.nfj.trend, ncol = 1)
  print(annotate_figure(gp, top = text_grob(modelfit.output$plot.title, size = 18)))
  
  plot(NA, xlim = c(0,5), ylim = c(0,5), bty = "n", 
       xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  text(1,4, modelfit.output$diagnostics.text["heading"], pos = 4, cex = 3)
  text(1,3, modelfit.output$diagnostics.text["title"], pos = 4, cex = 2)
  text(1,2, modelfit.output$diagnostics.text["formula"], pos = 4, cex = 2)
  text(1,1.5, modelfit.output$diagnostics.text["sources"], pos = 4)
  text(1,1, modelfit.output$diagnostics.text["knots"], pos = 4)
  
  gp = ggarrange(modelfit.output$fig.hist, modelfit.output$fig.counts, ncol=1)
  print(annotate_figure(gp, top = text_grob(modelfit.output$plot.title, size = 18)))
  
  ## activity curves for highest data denstiy, highest estimated abundance
  print(modelfit.output$fig.activity.maxdata + ggtitle(
    paste0(modelfit.output$fig.activity.maxdata$labels$title, "\n(Location of highest data density)")
  ))
  print(modelfit.output$fig.activity.maxabund + ggtitle(
    paste0(modelfit.output$fig.activity.maxabund$labels$title, "\n(Location of highest estimated abundance)")
  ))
  dev.off()
  ## save jpgs for faster viewing of maps etc.
  ggsave(here(paste0("4_res/fit-summaries/jpgs-fastviews/",
                     code.cur, "-", run.suffix,
                     "-V", use.num, "-abund-map.jpg")),
         modelfit.output$fig.abund, 
         width = 15, height = 12)
  ggsave(here(paste0("4_res/fit-summaries/jpgs-fastviews/",
                     code.cur, "-", run.suffix,
                     "-V", use.num, "-trend-map.jpg")),
         modelfit.output$fig.trend, 
         width = 15, height = 12)
  ggsave(here(paste0("4_res/fit-summaries/jpgs-fastviews/",
                     code.cur, "-", run.suffix,
                     "-V", use.num, "-trend-NFJ.jpg")),
         modelfit.output$fig.nfj.abund, 
         width = 15, height = 12)
  cat("\n")
  return(path.use)
}

