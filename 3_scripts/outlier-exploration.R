## looking at outliers
library(tidyverse)
library(here)
dat = read_csv(here("2_data_wrangling/cleaned-data/cleaned-data-aggregated.csv"))
tail(sort(dat$count),100)


hist(tail(sort(dat$count),1000))
plot(tail(sort(dat$count),500))

sum(dat$count>=5000, na.rm=T)


## let's take a quick look at everything 5k+

temp = dat %>% 
  filter(dat$count >= 5000) %>% 
  select(count, name, date, source, state, effort.universal) %>% 
  arrange(count)

write_csv(temp, "Outliers.csv")
