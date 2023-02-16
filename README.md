# status-of-butterflies-analysis

## Overview
This Rproject contains code to work with the SBUS data and fit individual taxa (based on taxa code) with a variety of generalized additive models (GAMs). 
Note that while github is handy for sharing (and updating) code, our data is protected by a range of MOUs, and we CANNOT host it in a publicly available space.
Additionally, the full data file is half a Gb, which is not github-friendly. 

## Usage

If you don't already have access to the dropbox folder with the cleaned data, contact me (edwards.evoeco [ at ] gmail.com). 
After pulling or downloading + unzipping this github repository, 

 - copy `cleaned-data-aggregated.csv` from the shared dropbox folder into `2_data_wrangling/cleaned-data/`
 - copy `trip-codes.csv` and `trip-templates.csv` from the shared dropbox folder into `2_data_wrangling/trip-wrangling/`

You should now be good to go. 

I have recently (2/16/2023) updated the structure of the code to simplify carrying out multiple runs. Parameters (including model formula!) are specified in an R file in `3_scripts/parameter-files/` (create a new file for each run you do).

`model-runner-parameterized.R` runs the specified model. Update the `run.name` parameter to be the name of the parameter file used (without the `.R` suffix), and then run `model-runner-parameterized.R` to carry out the model. I may further modify this to simplify running multiple parameter files at once, but I'm wary of running multiple simulations in the same environment.

Results from a run of the model are saved in `fits-parameterized`, in a folder based on the `run.name` (and so the original parameter file). A "-V#" will be tacked onto the end, specifying a number that hasn't been used yet. This will prevent accidental overwriting of simulations if a parameter file is changed but not renamed. Within a run's folder, 

- `parameters-used.R` is a copy of the parameter file used
- each taxa has a separate .html report with heatmaps of abundance and trends and some additional diagnostics + information, a separate .RDS file with the fitted model for that taxa
- AAA-run-summary-trends.csv contains information about the estimated trends for each taxa - correlation between the fitted model predictions and actual Naba NFJ sites, estimated trends for the species, estimated trends in the location with the highest estimated density (when the model goes wonky, it often estimates very high density, so this location is important for diagnostics), and the estimated trends in the location with the highest density of NFJ data.
- AAA-run-summary - trends - ongoing.csv saves the current version of AAA-run-summary-trends.csv with each new taxa, as a failsafe in case of R crashing or computer losing power.
-  AAA-run-summary.txt contains summary information for the run overall, including taxa used, regions used, parameters, etc.
- `heatmaps/` contains a copy of the trend and abundance heatmaps for each taxa, for use in Rshiny apps or presentations. By default, these are not generated (it adds to the runtime). Set `copy.heatmaps = TRUE` in the parameter file to generate them.

I will update an example using a single species in `3_scripts/model-fitting-example.R` to look at the individual pieces of the model-fitting and results-creation process. 


## Git ignore details

The .gitignore files should ensure that the raw data and the species-filtered data (generated at runtime) in 2_data_wrangling won't be commited to git and pushed to github. As a reminder, cleaned-data-aggregated.csv SHOULD NOT BE PUSHED TO GITHUB OR SHARED, and the individual data files generated in `2_data_wrangling/cleaned by code/` similarly SHOULD NOT BE PUSHED TO GITHUB OR SHARED.

I've included \*.pdf and \*.jpg in the overall .gitignore file because those files will rapidly clog up the repository, and you don't want github trying save the history of binary files like pdfs.

## Visualizing data coverage

To facilitate looking at data coverage across space and time, I've made a shiny app: https://edwards-evoeco.shinyapps.io/sbus-data-coverage/. This provides summaries of data coverage at relatively broad scales (ie nearest whole lat/lon, etc), and so does not violate our MOUs. `data-density-mapping.R` contains the code to generate a pdf equivalent of this for all taxa. This is probably not something that will need to be rerun in the future, but it's convenient to store the code here with everything else. 
