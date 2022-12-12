# status-of-butterflies-analysis
Shareable code for SBUS analysis team.

## Overview
This Rproject contains code to work with the SBUS data and fit individual taxa (based on taxa code) with a variety of generalized additive models (GAMs). 
Note that while github is handy for sharing (and updating) code, our data is protected by a range of MOUs, and we CANNOT host it in a publicly available space.
Additionally, the full data file is half a Gb, which is not github-friendly. 

## Usage

If you don't already have access to the dropbox folder with the cleaned data, contact me (edwards.evoeco [ at ] gmail.com). 
After pulling or downloading + unzipping this github repository, 
 - copy `cleaned-data-aggregated.csv` from the shared dropbox folder into 2_data_wrangling/cleaned-data
 - copy `trip-codes.csv` and `trip-templates.csv` from the shared dropbox folder into 2_data_wrangling/trip-wrangling

You should now be good to go. `3_scripts/species-loop-prototype.R` walks through fitting a single taxa (cabbage white),
and has code to loop over any given taxa codes and generate pdf reports of the results in 4_res/fit-summaries. `species-loop-prototype.R`
is flexible, and should accomodate experimenting with a wide range of model formulation and data filtering.

