# QUICK NOTE: Noam Ross suggested that gamm4 package may be a bigger speedup when 
# the limiting factor is random effects. I don't have currently ahve the time to try rewriting
# using a new package (2 days to flying to CO), but I should look at this next chance I get.
# 
# https://peerj.com/articles/6876/#fig-19

## To facilitate running multiple versions of the model, this will read in a parameter file
## from 3_scripts/parameter-files and run it. Results will be saved to 4_res/fits-parameterized
## 
## ## Collin will be using this as his main workspace. model-tester will be a more accessible tutorial space.

rm(list = ls(all.names = TRUE))
gc() 

### Parameter file name should be the only piece that needs in THIS script. Save parameters in an R file in 3_scripts/parameter-files
### and provide the name (without ".R") here

run.name = "full-run-all-sources-presence-absence"

## If running on Collin's desktop
# setwd("G:/repos/status-of-butterflies-analysis") ## for running out of rstudio, specific to Collin's desktop


## QUICK NOTE FOR DEBUGGING FRUSTRATING PROBLEM!!
## If you want to run this in base R rather than Rstudio, you will need 
## pandoc the stand-alone program, not just pandoc the R package.
## If you repeatedly try to install and uninstally pandoc the package, it will not solve your problems.
## No amount of R package debugging will fix this. 
## On the other hand, the free (GNU license free) program is available here: https://pandoc.org/

## Libraries and sourcing functions ------------

library(here)
library(tidyverse)
library(mgcv)
library(sf)
library(sp)
library(geometry)
library(viridis)
library(scales)
library(ggpubr)
library(geometry)
library(formula.tools)
source(here("3_scripts/funs.R"))
source(here("3_scripts/funs-fitting.R"))
source(here("3_scripts/model-fitting-function.R"))

### reading parameter file ----------

sim.start = proc.time()

source(here(paste0("3_scripts/parameter-files/", run.name, ".R")))

### creating results directory ----------

## Identify run number to avoid accidental over-writing with duplicated parameter file names
## Quick and dirty - update results species by species, given crashes
cur.files = list.files(here(paste0("4_res/fits-parameterized/")))
cur.file.code = cur.files[grepl(run.name, cur.files)]
if(length(cur.file.code)>0){
  cur.nums = as.numeric(gsub(paste0(run.name, "-V"), "", cur.file.code))
  use.num = max(cur.nums)+1
}else{
  use.num = 1
}
res.path = here(paste0("4_res/fits-parameterized/", run.name, "-V",use.num))
dir.create(res.path)
dir.create(paste0(res.path,"/heatmaps"))
file.copy(from = here(paste0("3_scripts/parameter-files/", run.name, ".R")),
          to = here(paste0(res.path, "/parameters-used.R")))
## actual model fitting and saving --------
summary.df = NULL
for(i.spec in 1:nrow(specs.do)){
    code.cur = specs.do$code[i.spec]
    # try({
    out = model_runner(code.cur = code.cur,
                       form.use = form.use,
                       knot_maker = knot_maker,
                       sitere_maker = sitere_maker,
                       regions.dict = regions.dict,
                       fit.family = fit.family,
                       use.inferred = use.inferred,
                       geography.constrain = geography.constrain,
                       use.only.source = use.only.source,
                       n.threads.use = n.threads.use,
                       use.range = use.range,
                       do.pheno = do.pheno)
    ## NOTE: I've updated the report_maker() to reduce the number of plots it creates.
    ## This step is a sizeable chunk of the run-time, and creating things we don't use
    ## wastes time.
    cat("creating report\n")
    output.name = report_maker(out,
                               code.cur = code.cur,
                               res.path = res.path,
                               copy.heatmaps = copy.heatmaps)
    res.cur =  data.frame(
      code = code.cur, 
      abund.correlation = out$abund.cor,
      species.trend = out$trend.species,
      trend.at.highest.predictions = out$trend.highabund,
      trend.at.best.nfj.data = out$trend.mostnfj,
      filename = output.name,
      n.data = out$n
    )
    summary.df = rbind(summary.df,
                       res.cur)
    ## Quick and dirty - update results species by species, given crashes
    write.csv(summary.df,
              here(paste0(res.path,"/AAA-run-summary - trends - ongoing.csv")),
              row.names = FALSE, append = TRUE)
  # },
  # outFile = file.create(paste0("res.path/","error-log.txt"))
  # ) 
}

## saving summary file ----------------
if(do.summary){
  write.csv(summary.df,
            here(paste0(res.path,"/AAA-run-summary-trends.csv")),
            row.names = FALSE)
  ## save text here:
  path.use = here(paste0(res.path,"/AAA-run-summary.txt"))
  
  parm_paster = function(parm.name){
    paste0(parm.name,":\n", "  ", as.character((get(parm.name))),"\n")
  }
  
  cat(c(paste0("Fit on: ", Sys.time(), "\n"),
        "Parameters:\n==================",
        parm_paster("form.use"),
        parm_paster("use.inferred"),
        parm_paster("geography.constrain"),
        parm_paster("use.only.source"),
        "\n"),
      sep = "\n",
      file = path.use)
  
  cat("Correlation summaries:\n==================\n",
      append = TRUE,
      file = path.use)
  
  write_delim(summary.df[,c("code","abund.correlation")], 
              file = path.use, col_names = TRUE, append = T,
              delim = "\t")
  
  cat(c("\n\n Report files generated:\n==================",
        summary.df$filename),
      sep = "\n",
      file = path.use,
      append = TRUE)
  
  cat("\n\n Human notes on model fit:\n==================\n[Fill in as desired]\n\n\n",
      file = path.use,
      append = TRUE)
  cat(c("============================================================",
        "============================================================",
        "============================================================", "\n\n"),
      file = path.use,
      sep = "\n",
      append = TRUE)
  cat("Extra details:  knot and site-making functions, region dictionary\n==================\n",
      file = path.use,
      append = TRUE)
  
  
  sink(path.use, append = T)
  cat("\n\nknot_maker():\n")
  print(knot_maker)
  cat("\n\nsitere_maker():\n")
  print(sitere_maker)
  cat("\n\nRegion dictionary:\n")
  sink()
  write_delim(regions.dict, 
              file = path.use, col_names = TRUE, append = T, delim = "\t\t\t")
}

cat("\nOverall runtime (minutes):\n")
print((proc.time()-sim.start)/60)

