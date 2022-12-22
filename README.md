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

You should now be good to go. `3_scripts/model-fitting-example.R` walks through fitting a single taxa (cabbage white), looking at individual pieces of the model-fitting and results-creation process. `3_scripts/model-tester.R` fits an any given set of taxa with user-provided mgcv formula and data filtering parameters, leaning heavily on the functions in `model-fitting-function.R`. `model-fitting-function.R` also saves a summary text file with the parameters and species run, for easy tracking of which parameter combinations have been tested.

The .gitignore files should ensure that the raw data and the species-filtered data (generated at runtime) in 2_data_wrangling won't be commited to git and pushed to github. As a reminder, cleaned-data-aggregated.csv SHOULD NOT BE PUSHED TO GITHUB OR SHARED, and the individual data files generated in `2_data_wrangling/cleaned by code/` similarly SHOULD NOT BE PUSHED TO GITHUB OR SHARED.

I've included *.pdf and *.jpg in the overall .gitignore file because those files will rapidly clog up the repository, and you don't want github trying save the history of binary files like pdfs.

## Visualizing data coverage

To facilitate looking at data coverage across space and time, I've made a shiny app: https://edwards-evoeco.shinyapps.io/sbus-data-coverage/. This provides summaries of data coverage at relatively broad scales (ie nearest whole lat/lon, etc), and so does not violate our MOUs. `data-density-mapping.R` contains the code to generate a pdf equivalent of this for all taxa. This is probably not something that will need to be rerun in the future, but it's convenient to store the code here with everything else. 
