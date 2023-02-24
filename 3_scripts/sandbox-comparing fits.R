#comparing model fits when fit with ALL sources, vs when fit with just NFJ.
# Note that these are AAA-run-summary - full-run-only-NFJ-V1.trends.csv
# and   AAA-run-summary - full-run-all-sources-V2-trends.csv
library(here)
library(tidyverse)

dat.allsource = read_csv(here("4_res/fit-summaries/AAA-run-summary - full-run-all-sources-V2-trends.csv"))
dat = dat.allsource %>% 
  select(code, abund.correlation, species.trend)
dat.NFJ = read_csv(here("4_res/fit-summaries/AAA-run-summary - full-run-only-NFJ-V1-trends.csv")) %>% 
  select(code, abund.correlation, species.trend) %>% 
  rename(nfjonly.abund.correlation = "abund.correlation",
         nfjonly.species.trend = "species.trend")
dat = inner_join(dat, dat.NFJ)

ggplot(dat, aes (x = species.trend, y = nfjonly.species.trend))+
  geom_abline(yintercept = 0, slope = 1, linetype = 2) +
  geom_point()

## Create a table of species by error
dat = dat %>% 
  mutate(error = nfjonly.species.trend - species.trend) %>% 
  mutate(abserror = abs(error)) %>% 
  arrange(desc(abserror))
## want to integrate dictionary
dict.spec = read_csv(here("2_data_wrangling/dictionaries/dictionary-use.csv")) %>% 
  group_by(code) %>% 
  summarize(sci.name = paste0(name, collapse = " | "))
dict.com = read_csv(here("2_data_wrangling/dictionaries/code-to-common-name.csv")) %>% 
  select(code, common) %>% 
  unique()

## quick and dirty

dat.clean = left_join(dat, dict.com)
dat.clean = left_join(dat.clean, dict.spec)

write_csv(dat.clean, here("4_res/NFJ-vs-all-sources-comparison-2-15-22.csv"))

## what are the species we COULDN'T fit? 
problem.taxa = data.frame(code = c("CALLDUM", "CALLPOL", "CHLOPAL", "ERYZAR", "EUPYBIM"))
problem.taxa = left_join(problem.taxa, dict.com)
problem.taxa = left_join(problem.taxa, dict.spec)

write_csv(dat.clean, here("4_res/NFJ-vs-all-sources-problem-taxa-2-15-22.csv"))
