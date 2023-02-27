# status-of-butterflies-analysis

## Overview
This Rproject contains code to work with the SBUS data and fit individual taxa (based on taxa code) with a variety of generalized additive models (GAMs). 
Note that while github is handy for sharing (and updating) code, our data is protected by a range of MOUs, and we CANNOT host it in a publicly available space.
Additionally, the full data file is half a Gb, which is not github-friendly. 

## Setup

Because of MOUs, we cannot share the raw data online. If you don't already have access to the dropbox folder with the cleaned data, contact me (edwards.evoeco [ at ] gmail.com). 
After pulling or downloading + unzipping this github repository, 

 - copy `cleaned-data-aggregated.csv` from the shared dropbox folder into `2_data_wrangling/cleaned-data/`
 - copy `trip-codes.csv` and `trip-templates.csv` from the shared dropbox folder into `2_data_wrangling/trip-wrangling/`

You should now be good to go. 

## Usage 

I have recently (2/16/2023) updated the structure of the code to simplify carrying out multiple runs. Typical workflow is now to specify parameters, model formula, any complications in knot placement + random effect variables, and the list of species to be fit all in a separate file, placed in `3_scripts/parameter-files/` (I highly recommend creating a new file for each time of run you do, rather than using a single file and changing the parameters etc.).

Then in `scripts/model-runner-parameterized.R`, change the line defining `run.name` to match the name of your parameter file (without the `.R` suffix). At this point, you can run the file, and it will fit and save results for all species. I may further modify this to simplify running multiple parameter files at once, but I'm wary of running multiple simulations in the same environment. See the `parameters` section below for a breakdown of all the pieces.

### Results

Results from a run of the model are saved in `fits-parameterized`, in a folder based on the `run.name` (and so the original parameter file). A "-V#" will be tacked onto the end, specifying a number that hasn't been used yet. This will prevent accidental overwriting of simulations if a parameter file is changed but not renamed. Within a run's folder, 

- `parameters-used.R` is a copy of the parameter file used
- each taxa has a separate .html report with heatmaps of abundance and trends and some additional diagnostics + information, a separate .RDS file with the fitted model for that taxa
- AAA-run-summary-trends.csv contains information about the estimated trends for each taxa - correlation between the fitted model predictions and actual Naba NFJ sites, estimated trends for the species, estimated trends in the location with the highest estimated density (when the model goes wonky, it often estimates very high density, so this location is important for diagnostics), and the estimated trends in the location with the highest density of NFJ data.
- AAA-run-summary - trends - ongoing.csv saves the current version of AAA-run-summary-trends.csv with each new taxa, as a failsafe in case of R crashing or computer losing power.
-  AAA-run-summary.txt contains summary information for the run overall, including taxa used, regions used, parameters, etc.
- `heatmaps/` contains a copy of the trend and abundance heatmaps for each taxa, for use in Rshiny apps or presentations. By default, these are not generated (it adds to the runtime). Set `copy.heatmaps = TRUE` in the parameter file to generate them.

I have an example using a single species in `3_scripts/model-fitting-example.R` to look at the process in a more granular way. 

### Parameters

In your parameter file, you should set the following parameters:
- `do.summary` - boolean, if `TRUE`, `model-runner-paramterized.R` saves a text file with summary information after fitting all specified species
- `use.inferred` - boolean. If `TRUE`, infer zeroes from the absence of reporting from programs that report all butterfly species seen (and only butterfly species seen). Unless you have a strong reason to do otherwise, use `TRUE`
- `fit.family` - character string. What distribution should be used when fitting the data. While this technically allows us to fit a presence/absence model by specifying "binomial", several components of the modeling process currently (as of 2/26/23) rely on a log-link function. Unless you have strong reason to do otherwise, use `nb` (for negative binomial).
- `use.range` - boolian. Filter our data to only use observations within the ranges provided by Eliza Grames. Because inferring zeroes creates zeroes across the continent, this is *very* useful! Note that there are a few taxa (3) that do not currently have range maps. We're working to resolve this.
- `geography.constrain` - boolian. Filter our data to only include zero-count observations within the convex hull defined by non-zero points. This has been eclipsed by the much more elegant range maps. Use `use.range` instead, unless you have a very clear reason not to.
- `use.only.source` - NULL or vector of character strings. Specify which source or sources of data to use in the analysis. If NULL, uses data from all sources (a good default).
- `copy.heatmaps` - boolean. Should the heatmaps of abundance trend be saved outside of the individual reports? Mostly useful for creating Rshiny apps or other reports that need these images.
- `n.threads.use` - numeric. how many cores should be used when fitting the gams? More threads = faster, if the computer has the cores and ram for it.
- `do.pheno` - boolean. Calculations of abundance and trends generally require prediction of counts across a fine grid of day of year, which is computationally intensive. However, this is relevant only in cases when the `doy` term is included in the model. For models that do not include phenology in this way, abundance is identical across all days, and a prediction per year is sufficient. Set `do.pheno` to `TRUE` in these cases (and ONLY in these cases) to substantially speed up analyses. 

The model formula itself should be stored as `form.use`, and created using the `formula()` function. This will be run by `mgcv::bam()`, an optimized version of `mgcv::gam()`, and should use mgcv nomenclature. Any variable in the cleaned data can be included. Additionally you can use `site.refac`, which is not in the cleaned data and instead is created by the provided `sitere_maker()` (described below).

You should also specify two functions:
- `sitere_maker()` - this function takes the provided data set, and is used to create the variable `$site.refac` to serve as a random effect. To make `site.refac` serve as a normal site-level random effect, merely return site as a factor. To speed up model fitting and potentially avoid numerical issues, `sitere_maker()` can first lump together sites with only a few observations - these sites are unlikely to be assigned estimates that are far from zero, but drammatically inflate our count of sites, complicating model fitting procedures. See `full-run-all-sources.R` for an example of this.
- `knot_maker()` - this function takes the provided data set, and allows for specification of non-default knot locations. Typically you would want to do this using quantiles (rather than absolute lat/lon). This is likely unnecessary now that we can clip our data to range maps. If unneeded, have the function return an empty list; see `full-run-all-sources.R` for an example of this.

All species codes given in `specs.do$code` will be run. For simple test cases, you might create specs.do with a few manually identified species. For full runs, I recommend reading in the full data set, and identifying all species codes that meet desired data specifications (e.g. number of observations, number of sites, etc). See `3_scripts/parameter-files/full-run-all-sources.R` for an example. 

## Function descriptions

(*in progress*)

### `model_runner` 
Runs the model and returns various results for a given species. This includes creating and saving the data for that species, fitting the provided model, and predicting trends and abundance. This function takes the parameters listed above, as well as infer.messy.levels, which specifies how to handle "unknown in X genus" or "X family". For whatever levels are specified, during the "inferring zeroes" step in the code, trips that didn't report the focal species but *did* report an unknown at the specified level will be skipped over.

Results structure:
- `plot.title`: pretty name for the plot
- `loc.plot`: data frame of lon, lat, estimated abundance averaged across all years (`abund.index`), estimated growth rate `gr.med`
- `abund.cor`: correlation in abundance between the average NFJ counts for each site across the last 10 years, and the average of the model predictions for those counts.
- `fig.abund`: ggplot map of estimated abundance
- `fig.trend`: same, but for trends
- `fig.nfj.abund`: plot showing a more detailed breakdown of the abund.cor value. Each point is a site.
- `fig.hist`: histogram of observtions across day of year, broken down by whether or not they were inferreed, and colored by whether or not they were non-zero. Potentially useful for diagnostics.
- `fig.counts`: breakdown of counts ~ doy for each region. Potentially useful for diagnostics.
- `diagnostics.text`: vector of key aspects of model fit, for ease of report creation.
- `fitted.model`: actual fitted gam object.
- `trend.species`: estimate of the species-level growth rate. Weighted average across the map (weighted by estimated abundance)
- `trend.highabund`: estimated trend at the location of highest estimated abundance. Potentially useful for diagnostics.
- `trend.mostnfj`: estimated trend at the location of highest average NFJ counts. Potentially useful for diagnostics.
- `n`: total number of data points used.
- `events.missed.messy`: number of zeroes NOT inferred because of the given infer.messy.levels.

## Git ignore details

The .gitignore files should ensure that the raw data and the species-filtered data (generated at runtime) in 2_data_wrangling won't be commited to git and pushed to github. As a reminder, cleaned-data-aggregated.csv SHOULD NOT BE PUSHED TO GITHUB OR SHARED, and the individual data files generated in `2_data_wrangling/cleaned by code/` similarly SHOULD NOT BE PUSHED TO GITHUB OR SHARED.

I've included \*.pdf and \*.jpg in the overall .gitignore file because those files will rapidly clog up the repository, and you don't want github trying save the history of binary files like pdfs.

## Visualizing data coverage

To facilitate looking at data coverage across space and time, I've made a shiny app: https://edwards-evoeco.shinyapps.io/sbus-data-coverage/. This provides summaries of data coverage at relatively broad scales (ie nearest whole lat/lon, etc), and so does not violate our MOUs. `data-density-mapping.R` contains the code to generate a pdf equivalent of this for all taxa. This is probably not something that will need to be rerun in the future, but it's convenient to store the code here with everything else. 
