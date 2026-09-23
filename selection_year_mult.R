library(tidyverse)
library(broom)

# save(traitnames, year_pal,
#      fitnessnames_leaf, leaftraits, phystraits, lt.plantyr, pt.plantyr, #leaf data from maxfield_data.R
#      fitnessnames_floral, floraltraits, meta.pca.nosepal, mnps.plantyr.traits, #floral data from maxfield_volatiles_2020.Rmd
#      file="data/selection_gradient_inputs.rda")
load("data/selection_gradient_dfs.rda")

# Leaf traits -------------------------------------------------------------

# function to calculate standardized selection gradients in each year
# only includes subplots with unmanipulated snowmelt and precipitation
# fitnesstrait: the fitness metric (survived, flowered, RGR), observed the following year
# set: lt1/lt2 = leaf traits (trichomes, SLA, water content), collection round 1 or 2
# leaves from 2019 round 2 were not weighed, so this round is excluded
# set: pt = physiology traits (conductance, photosynthesis, WUE)
# pt and fitness measured on less than 10 plants in 2019 and 2020, excluded
# the multiple regression includes only the three traits from that set,
# since each set of traits was measured on different plants
# fitness and traits are standardized within each year
# uses datasets with the mean trait for each round and year (lt.plantyr, pt.plantyr)
selection_year_leaf <- function(fitnesstrait, set, .year) {
  
  traits <- list(lt1=leaftraits, lt2=leaftraits, pt=phystraits)[[set]]
  
  t.plantyr <- list(lt1 = filter(lt.plantyr, round=="1"), 
                    lt2 = filter(lt.plantyr, round=="2"),
                    pt = pt.plantyr)[[set]] %>% 
    filter(year == .year, water == "Control", snow == "Normal") %>% 
    rename("fitness"=fitnesstrait) %>% select(fitness, !!!traits) %>% 
    drop_na(all_of(traits), fitness) %>% 
    mutate(across(all_of(traits), \(x) x/sd(x)), # divide traits by standard deviation
           relfitness = fitness/mean(fitness)) # divide fitness by mean
  
  if(nrow(t.plantyr) < 10) return(NA) # not enough data for that round
  
  lm(formula(paste("relfitness ~ ", paste(traits, collapse="+"))), data=t.plantyr)
}

# the combinations of fitness trait, trait set, and year to run models on
selection.combos.leaf <- expand_grid(fitnesstrait = names(fitnessnames_leaf), 
                                set = c("lt1", "lt2", "pt"),
                                .year = names(year_pal))

# run the models 
coef.selection.leaf <- selection.combos.leaf %>%
  mutate(model = pmap(., selection_year_leaf)) %>% 
  filter(!is.na(model)) %>% rename(year = .year) %>% 
  mutate(coefs = map(model, tidy), 
         n_obs = map_int(model, nobs),
         round = c(lt1=1, lt2=2, pt=1)[set]) %>% 
  select(-model) %>% unnest(coefs) %>% filter(term!="(Intercept)")

#plot the gradients
ggplot(coef.selection.leaf, aes(x=paste(year, round, sep="."), y=estimate, color=year, size=n_obs,
                           ymax=estimate+std.error, ymin=estimate-std.error)) +
  facet_grid(vars(fitnesstrait), vars(factor(term, levels = names(traitnames))), 
             labeller = as_labeller(c(traitnames, fitnessnames_leaf)),
             scales = "free_y") + 
  geom_hline(yintercept=0) + geom_pointrange() + 
  scale_color_manual(values=year_pal) + scale_size_area(max_size=2) + 
  labs(y="SD-standardized selection gradient", x="Year, round", color="Year", size="Sample size") +
  theme_minimal() + theme(legend.position = "top", axis.text.x = element_text(angle=90))

# Floral traits -----------------------------------------------------------

# function to calculate standardized selection gradients in each year
# only includes subplots with unmanipulated snowmelt and precipitation
# fitnesstrait: the fitness metric, observed the same year
# uses datasets with the mean trait for each year (meta.pca.nosepal, mnps.plantyr.traits)
selection_year_floral <- function(fitnesstrait, regtraits, .year) {
  t.plantyr <- meta.pca.nosepal %>% 
    filter(year == .year,  water == "Control", snow == "Normal") %>% 
    select(-any_of(fitnesstrait)) %>% 
    left_join(mnps.plantyr.traits %>% select(year, plotid, plantid, all_of(fitnesstrait)), 
              by=c("year", "plotid", "plantid")) %>% 
    select("fitness"=all_of(fitnesstrait), all_of(regtraits)) %>% 
    drop_na(fitness, all_of(regtraits)) %>% 
    mutate(across(all_of(regtraits), \(x) x/sd(x)),
           relfitness = fitness/mean(fitness)) %>% 
    select(where(~sum(is.na(.))==0))
  
  if(nrow(t.plantyr) < 10) return(NA) # not enough data for that round
  
  lm(formula(paste("relfitness ~", paste0(regtraits, collapse="+"))), data = t.plantyr) 
}

# the combinations of fitness trait and year to run models on
selection.combos.floral <- expand_grid(fitnesstrait=names(fitnessnames_floral), .year=names(year_pal))

# run the models 
coef.selection.floral <- selection.combos.floral %>%
  mutate(model = pmap(., selection_year_floral, 
                      regtraits=floraltraits[-c(4,8)])) %>% #exclude sepal width and flower number
  filter(!is.na(model)) %>% rename(year = .year) %>% 
  mutate(coefs = map(model, tidy), 
         n_obs = map_int(model, nobs)) %>% 
  select(-model) %>% unnest(coefs) %>% filter(term!="(Intercept)")

#plot the gradients
ggplot(coef.selection.floral, aes(x=year, y=estimate, color=year, size=n_obs,
                                  ymax=estimate+std.error, ymin=estimate-std.error)) +
  facet_grid(vars(fitnesstrait), vars(factor(term, levels = names(traitnames))), 
             labeller = as_labeller(c(traitnames, fitnessnames_floral))) + 
  geom_hline(yintercept=0) + geom_pointrange() + 
  scale_color_manual(values=year_pal) + scale_size_area(max_size=1.5) + 
  labs(y="SD-standardized selection gradient", x="Year", color="Year", size="Sample size") +
  theme_minimal() + theme(legend.position = "top", axis.text.x = element_text(angle=90))


bind_rows(coef.selection.leaf, coef.selection.floral) %>% 
  write_csv("data/selection_gradients_controls.csv")
  
