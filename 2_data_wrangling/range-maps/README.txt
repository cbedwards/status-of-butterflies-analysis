These are the range maps created by Dr. Eliza Grames and posted here: https://github.com/elizagrames/butterflyblobs.
I last updated these files on 2/23/2023, from commit 584ca60b86b636438b96a430311703bf793f6073 [584ca60]


To simplify use in analysis scripts, 3_scripts/range-maps-relabel relabels these files to their GU codes. Any GU code for which there isn't a 1:1 mapping between map name and GU dictionary entry are reported in range-map-mismatch.csv

===================================

Notes: 
Phoesta currently needs a conversation, as there are two species maps that match for this.

For PARRMAL, we have a typo in the dictionary, with a hyphen in malbum. I've added a new entry to the dictionary, but since the dictionary is a derived product, this may need to be updated upstream at some point. Check this point of failure if PARRYMAL 

PYROIL and PYRPHI are reported together in GBIF, so share a range map. Presently it does not appear that this joint blob is on the repo. 

OCHYUM uses the same range map as OCHSYL; I've added a line of code in the range map script to make a second copy for this code.