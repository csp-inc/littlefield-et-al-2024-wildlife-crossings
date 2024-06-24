library(tidyverse)
library(glmmTMB)
source('code/utils/crossval-utils.R')

# Get seasonal RSF top models and migration top model
load(file = "output/winter-summer-rsf-top-models.rda", verbose = TRUE)
load(file = "output/migration-ssf-top-model-v2.rda", verbose = TRUE)

summer_cv <- kfoldRSF(summer_top, k = 5, nrepet = 30,nbins = 10, 
                      random = TRUE, jitter = TRUE, reproducible = FALSE)
winter_cv <- kfoldRSF(winter_top, k = 5, nrepet = 30,nbins = 10, 
                      random = TRUE, jitter = TRUE, reproducible = FALSE)
mig_cv <- kfoldSSF(mod = mig_top, strata_name = "step_id", k = 5, 
                   nrepet = 30, reproducible = FALSE, details = TRUE)

save(list = c('summer_cv','winter_cv','mig_cv'), file = "output/rsf-ssf-cross-val.rds")

cv_summarize <- function(dat, mod_name){
  d <- dat %>% 
    filter(kfold != 0) %>% 
    group_by(type) %>% 
    summarise(mean = mean(kfold), sd = sd(kfold))
  
  d2 <- data.frame('Season' = mod_name, 
                   'Corr Observed' = paste0(round(d$mean[d$type == 'obs'], 2),
                                            " (", round(d$sd[d$type == 'obs'], 2), ")"),
                   'Corr Rand' = paste0(round(d$mean[d$type == 'rand'], 2),
                                        " (", round(d$sd[d$type == 'rand'], 2), ")"))
  return(d2)
}

out<- rbind(cv_summarize(winter_cv, "Winter"),
            cv_summarize(summer_cv, "Summer"),
            cv_summarize(mig_cv, "Migration"))
out %>% write_csv(file = 'output/cross-val-results.csv')
