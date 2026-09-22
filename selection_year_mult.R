library(tidyverse)
library(broom)
load("data/selection_gradient_dfs.rda")

# function to calculate standardized selection gradients in each year
# only includes subplots with unmanipulated snowmelt and precipitation
# response: the fitness metric (survived, flowered, RGR), observed the following year
# set: lt1/lt2 = leaf traits (trichomes, SLA, water content), collection round 1 or 2
# leaves from 2019 round 2 were not weighed, so this round is excluded
# set: pt = physiology traits (conductance, photosynthesis, WUE)
# pt and fitness measured on less than 10 plants in 2019 and 2020, excluded
# the multiple regression includes only the three traits from that set,
# since each set of traits was measured on different plants
# fitness and traits are standardized within each year
# uses datasets with the mean trait for each round and year (lt.plantyr, pt.plantyr)
selection_year_mult_mod <- function(response, set, .year) {
  
  traits <- list(lt1=leaftraits, lt2=leaftraits, pt=phystraits)[[set]]
  
  t.plantyr <- list(lt1 = filter(lt.plantyr, round=="1"), 
                    lt2 = filter(lt.plantyr, round=="2"),
                    pt = pt.plantyr)[[set]] %>% 
    filter(year == .year, water == "Control", snow == "Normal") %>% 
    rename("response"=response) %>% select(response, !!!traits) %>% 
    drop_na(all_of(traits), response) %>% 
    mutate(across(all_of(traits), \(x) x/sd(x)), # divide traits by standard deviation
           relfitness = response/mean(response)) # divide fitness by mean
  
  if(nrow(t.plantyr) < 10) return(NA) # not enough data for that round
  
  lm(formula(paste("relfitness ~ ", paste(traits, collapse="+"))), data=t.plantyr)
}

# the combinations of fitness trait, trait set, and year to run models on
selection.combos <- expand_grid(response = fitnesstraits, 
                                set = c("lt1", "lt2", "pt"),
                                .year = levels(treatments$year))

# run the models 
coef.selection <- selection.combos %>%
  mutate(model = pmap(., selection_year_mult_mod)) %>% 
  filter(!is.na(model)) %>% rename(year = .year) %>% 
  mutate(coefs = map(model, tidy), 
         n_obs = map_int(model, nobs),
         round = c(lt1=1, lt2=2, pt=1)[set]) %>% 
  select(-model) %>% unnest(coefs) %>% filter(term!="(Intercept)")

#plot the gradients
ggplot(coef.selection, aes(x=paste(year, round, sep="."), y=estimate, color=year, size=n_obs,
                           ymax=estimate+std.error, ymin=estimate-std.error)) +
  facet_grid(vars(response), vars(factor(term, levels = names(traitnames))), 
             labeller = as_labeller(c(traitnames, fitnessnames)),
             scales = "free_y") + 
  geom_hline(yintercept=0) + geom_pointrange() + 
  scale_color_manual(values=year_pal) + scale_size_area(max_size=2) + 
  labs(y="SD-standardized selection gradient", x="Year, round", color="Year", size="Sample size") +
  theme_minimal() + theme(legend.position = "top", axis.text.x = element_text(angle=90))
