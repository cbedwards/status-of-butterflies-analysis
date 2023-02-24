## aggregating and comparing results from specific analyses
library(here)
library(tidyverse)
library(ggforce)
res.full = read_csv(here("4_res/fits-parameterized/full-run-all-sources-V1/AAA-run-summary-trends.csv"))
res.full$run = "all sources"
res.ohio = read_csv(here("4_res/fits-parameterized/full-run-only-OhioLeps-V1/AAA-run-summary-trends.csv"))
res.ohio$run = "OhioLeps"
res.NFJ = read_csv(here("4_res/fits-parameterized/full-run-only-NFJ-V1/AAA-run-summary-trends.csv"))
res.NFJ$run = "NFJ"
res.MASSBfly = read_csv(here("4_res/fits-parameterized/full-run-only-MASSBfly-V1/AAA-run-summary-trends.csv"))
res.MASSBfly$run = "MASSBfly"
res.ill = read_csv(here("4_res/fits-parameterized/full-run-only-Illinois-BMN-V1/AAA-run-summary-trends.csv"))
res.ill$run = "Illinois BMN"

res.all = bind_rows(res.full, res.ohio, res.NFJ, res.MASSBfly, res.ill)
## give easy labels here
res.all$nlabel = paste0("(N = ", res.all$n.data, ")")
## add common names:
dict = read_csv(here("2_data_wrangling/dictionaries/code-to-common-name.csv"))
dict = dict %>% select(code, common) %>%  unique()
res.all = left_join(res.all, dict) %>% 
  mutate(label = paste0(code, ": ", common))


range(res.all$species.trend, na.rm=T) #Looks like trend estimates are within +/- 0.3

theme.larger =   theme(axis.title = element_text(size = rel(1.8)),
                       axis.text = element_text(size = rel(1.6)),
                       strip.text = element_text(size = rel(1.8)),
                       plot.title = element_text(size = rel(1.8)),
                       legend.text = element_text(size = rel(1.8)),
                       legend.title = element_text(size = rel(1.8))
)

theme.pdf =   theme(axis.title = element_text(size = rel(1)),
                    axis.text = element_text(size = rel(1)),
                    strip.text = element_text(size = rel(1)),
                    plot.title = element_text(size = rel(1)),
                    legend.text = element_text(size = rel(1)),
                    legend.title = element_text(size = rel(1)),
                    legend.position = "top"
)



nrow.use = 5
tot.plot = length(unique(res.all$label))
pdf(here("4_res/trends/trends-by-source.pdf"))
for(i in 1:ceiling(tot.plot/nrow.use)){ #total set
# for(i in 1:2){ #1:2 to fiddle with plotting options
  gp = ggplot(data = res.all, aes(x=species.trend, y = run, col = run, shape = run))+
    geom_vline(xintercept = 0, linetype=2)+
    geom_point(size = 3)+
    geom_text(aes(x = species.trend + 0.06*sign(species.trend), y = run, label = nlabel),
              show.legend = FALSE)+
    xlim(c(-0.3, 0.3))+
    xlab("Estimated trend")+
    ylab("")+
    facet_wrap_paginate(.~label, ncol = 1, nrow = nrow.use, page = i)+
    theme.pdf
  print(gp)
}
dev.off()
write_csv(res.all, here("4_res/trends/trends-by-source.csv"))
