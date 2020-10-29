# Load in all the Maxfield Meadow Ipomopsis data
# Includes treatments, loggers, census, phenology, floral and leaf traits, and seeds

# TODO write metadata  to /metadata for the final data files
# TODO add the following datasets:
#       + census 2018 - loaded but needs to be cleaned up and merged with other years
#       + phenology 2018, 2019 - loaded but need to be cleaned up
#       + floral volatiles 2018, 2019
#       + leaf physiology measures (LICOR) 2018, 2019, 2020
#       + seeds 2019, 2020

# setup -------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(RColorBrewer)

# The data is sourced from Google Drive.
# You will need to authenticate interactively with the gs4_auth() and drive_auth() commands to make this work locally.
# This creates an Oauth token in ~/.R/gargle for future use in any project (and avoids syncing this token to the web).
library(googlesheets4)
gs4_auth(email = T)
library(googledrive)
drive_auth(email = T)

# Get file listing for "2020 RMBL Campbell Lab"
datasheets <- drive_ls(as_id("1xVG466gMSCTfsMbvXudhq-Ruc2TKOZJB"), recursive=TRUE) 

# The experimental treatments. "water4" keeps all four treatments, "water" lumps control and mock rainout together.
treatments <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="treatments") %>%
  mutate(
    water4 = factor(recode(water, 
                           "control"="Control", "mock rainout"="Mock rainout", "addition"="Addition", "rainout"="Reduction")),
    water = factor(recode(water4, "Mock rainout"="Control")),
    snow = factor(recode(snow, "early"="Early","normal"="Normal"))) %>% 
  separate(plotid, c("plot","subplot"), sep=1, remove=F)

# Color palettes for water and snow treatments
water4_pal <- setNames(brewer.pal(9,name="Set1")[c(2,9,8,1)], levels(treatments$water4))
water_pal <- setNames(brewer.pal(9,name="Set1")[c(2,9,1)], levels(treatments$water))
snow_pal <- setNames(brewer.pal(3, name="Dark2")[c(2,1)], levels(treatments$snow))

# hobo --------------------------------------------------------------------

# Tarps were placed on plots 2, 4 and 5 in the spring to make the snow melt faster.
# HOBO pendant temperature/light loggers (model UA-002-xx?) were placed at the corners of each plot.
# HOBO manual: https://www.onsetcomp.com/files/manual_pdfs/9556-M%20UA-002%20Manual.pdf
# When the snow melts, the temperature rises above freezing and the light intensifies.
# The day of first bare ground is when these measurements exceed the following thresholds.
# Since sunlight may reach the logger through the snow, we will use the warm_date

melt_threshold_sun <- 15000 # light units (probably lux = lumen / m2)
melt_threshold_warm <- 1    # 1 degree Celsius

hobo <- drive_download(filter(datasheets, name=="maxfield_hobo_data.csv"), overwrite = T)$local_path %>% read_csv() %>% 
  mutate(time=with_tz(time, "America/Denver"),
         sun = intensity > melt_threshold_sun, 
         warm = temp > melt_threshold_warm) %>% 
  mutate_at(c("year","plot"), as.factor)

meltdates <- full_join(
  hobo %>% filter(month(time)%in%4:6, sun==TRUE) %>% 
    group_by(year, plot) %>% summarize_at("time", min) %>% 
    mutate(time=yday(time)) %>%  rename(sun_date=time),
  hobo %>% filter(month(time)%in%4:6, warm==TRUE) %>% 
    group_by(year, plot) %>% summarize_at("time", min) %>% 
    mutate(time=yday(time)) %>%  rename(warm_date=time))
meltdates <- left_join(meltdates, treatments %>% select(plot, snow) %>% distinct)

# weather -----------------------------------------------------------------

# Weather data from billy barr's RMBL station in Gothic, CO, USA provided by the Western Regional Climate Center
# Original (not FPA) data from 2017-2020 in metric units
# Downloaded as xls (csv download broken) from https://wrcc.dri.edu/cgi-bin/rawMAIN.pl?corbil
# Columns renamed from WRCC names according to /metadata/weather_wrcc_metadata.csv

wx <- drive_download(filter(datasheets, name=="billy_rmbl_wrcc_weather_2017_2020.csv"), overwrite = T)$local_path %>% 
  read_tsv(col_types="Dtnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn")
daily_precip <- wx %>% group_by(date) %>% summarize_at("precip_mm", sum) %>% mutate(year=year(date), julian=yday(date))
daily_temp <- wx %>% group_by(date) %>% summarize_at("av_temp_2m_C", mean, na.rm=T) %>% mutate(year=year(date), julian=yday(date))

# soil --------------------------------------------------------------------

# Soil moisture readings were taken with a Campbell Scientific HydroSense II https://www.campbellsci.com/hs2
# For each subplot the corners and middle of the plot were tested with 12 cm (or 20 cm?) rods
# The units are volumetric water content expressed as a percentage
# Data from 2018-06-06 and 2018-06-07 did not include all plots and are excluded

sm_sheets <- filter(datasheets, name=="2020 Maxfield Soil Moisture")
sm <- bind_rows(read_sheet(sm_sheets, sheet="2020"),
                read_sheet(sm_sheets, sheet="2019"),
                bind_rows(
                  read_sheet(sm_sheets, sheet="2018 - long"),
                  read_sheet(sm_sheets, sheet="2018 - wide") %>% 
                    pivot_longer(UR:M, names_to="location", values_to="VWC"))) %>% 
  mutate(plot = as.character(plot), plotid = paste0(plot, subplot), year=factor(year(date))) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor) %>% 
  filter(!(format(date) %in% c("2018-06-06", "2018-06-07")))

# Means within each subplot across days, months, and years
sm.subplot <- sm %>% group_by(year, date, plot, subplot, plotid, water, snow) %>% 
  summarize(VWC = mean(VWC, na.rm=T), .groups="drop") 

sm.subplotmonth <- sm.subplot %>% mutate(mo = month(date)) %>% group_by(year, mo, plot, subplot, plotid, water, snow) %>% 
  summarize(VWC = mean(VWC, na.rm=T), .groups="drop")

sm.subplotyear <- sm.subplot %>% group_by(year, plot, subplot, plotid, water, snow) %>% 
  summarize(VWC = mean(VWC, na.rm=T), .groups="drop")

# census ------------------------------------------------------------------

# Census of all plants at the beginning of each season to mark if they are dead, vegetative, or flowering
# For vegetative plants counted number of leaves on each rosette and longest leaf

# After loading the censuses, check the allowed states of notes/rosettes/flowering for each year and 
# [translate](https://docs.google.com/spreadsheets/d/17qSYLvl4aO3TNIHwz1YW_XdfjX_grSqoy38m4mKYXe8/edit#gid=816285957) 
# each combination into a simplified plant "status". If there is a problem, status is flagged "recheck". 

allowed_nrf <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="status")
status_priority <- c("flowering", "vegetative", "recheck",  "dead_nf", "tagnf") # determines order for deduplicating
translate_nrf_19 <- with(filter(allowed_nrf, year==2019),  setNames(status, nrf))
translate_nrf_20 <- with(filter(allowed_nrf, year==2020),  setNames(status, nrf))
status_ok <- c("flowering", "vegetative", "dead_nf")
status_pal <- setNames(brewer.pal(3, name="Set1")[c(1,3,2)], status_ok)

cen18 <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="2018") %>% 
  mutate_at(vars(ends_with(c("longest","leaves"))), ~as.integer(as.character(.))) %>% 
  mutate(plot = as.character(plot),
         plotid = paste0(plot, subplot),
         plant = as.character(plant),
         notes = other_notes %>% str_match("tagnf|nf|dead|chewed|eaten")  %>% factor %>% 
           fct_explicit_na("blank") %>% recode(eaten="chewed"),
         flowering = flowering %>% tolower %>% factor %>% fct_explicit_na("blank"),
         n_rosettes = as.integer(as.character(rosettes)),
         rosettes = recode(as.character(rosettes), ID = "indist", `0` = "zero" , `NULL`="blank",.default="one_or_more", ),
         r1_longest	= ifelse(rosettes=="zero", NA, r1_longest),
         r1_leaves = ifelse(rosettes=="zero", NA, r1_leaves)) %>% 
  left_join(treatments, by="plotid") %>%
  mutate_if(is.character, as.factor)

cen19 <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="2019foranalysis") %>% 
  mutate(plot = as.character(plot),
         plotid = paste0(plot, subplot),
         plant = as.character(plant),
         notes = other_notes %>% str_match("tagnf|nf|dead|chewed|eaten")  %>% factor %>% 
           fct_explicit_na("blank") %>% recode(eaten="chewed"),
         flowering = flowering %>% tolower %>% factor %>% fct_explicit_na("blank"),
         n_rosettes = as.integer(as.character(rosettes)),
         rosettes = recode(as.character(rosettes), ID = "indist", `0` = "zero" , `NULL`="blank",.default="one_or_more", ),
         status = recode(paste(notes, rosettes, flowering), !!!translate_nrf_19, .default="recheck") %>% fct_relevel(status_priority),
         r1_longest	= ifelse(rosettes=="zero", NA, r1_longest),
         r1_leaves = ifelse(rosettes=="zero", NA, r1_leaves),
         id = na_if(id, "-")) %>% 
  drop_na(id) %>%
  left_join(treatments, by="plotid") %>%
  mutate_if(is.character, as.factor)

cen20 <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="2020") %>% 
  select(-ends_with("_19")) %>% 
  mutate(plot = as.character(plot),
         plotid = paste0(plot, subplot),
         plant = as.character(plant),
         notes = notes %>% fct_explicit_na("blank"),
         flowering = flowering %>% as.character %>% fct_explicit_na("blank"),
         n_rosettes = as.integer(as.character(rosettes)),
         rosettes = recode(as.character(rosettes), ID = "indist", `0` = "zero" , `NULL`="blank",.default="one_or_more", ),
         status = recode(paste(notes, rosettes, flowering), !!!translate_nrf_20, .default="recheck") %>% fct_relevel(status_priority),
         id = na_if(id, "-"),
         plotid = paste0(plot, subplot)) %>% 
  drop_na(id) %>% 
  left_join(treatments, by="plotid") %>%
  mutate_if(is.character, as.factor)

cen <- full_join(cen20, rename_all(cen19, paste0, "_19"), by=c("id"="id_19")) 

# phenology ---------------------------------------------------------------

ph18 <- read_sheet(filter(datasheets, name=="2020 Maxfield Phenology"), sheet="2018") %>% 
  mutate(plot = as.character(plot), 
         subplot = toupper(subplot),
         plotid = paste0(plot, subplot),
         plantid = paste0(plotid, plant),
         julian = factor(yday(date)),
         open = rowSums(select(., starts_with("open")), na.rm=T),
         buds = rowSums(select(., starts_with("buds")), na.rm=T)) %>% 
  left_join(treatments, by=c("plotid","plot","subplot")) %>%
  mutate_if(is.character, as.factor) 

ph19 <- read_sheet(filter(datasheets, name=="2020 Maxfield Phenology"), sheet="2019") %>% 
  mutate(plot = as.character(plot), 
         plotid = paste0(plot, subplot),
         plantid = paste0(plotid, plant),
         plant = as.character(plant),
         julian = factor(yday(date))) %>% 
  left_join(treatments, by=c("plotid","plot","subplot")) %>%
  mutate_if(is.character, as.factor)

ph20 <- read_sheet(filter(datasheets, name=="2020 Maxfield Phenology"), sheet="2020") %>% 
  mutate(plot = as.character(plot), 
         plotid = paste0(plot, subplot),
         plantid = paste0(plotid, plant),
         plant = as.character(plant),
         julian = factor(yday(date)),
         open = rowSums(select(., starts_with("open")), na.rm=T),
         buds = rowSums(select(., starts_with("buds")), na.rm=T),
         eggs = rowSums(select(., starts_with("eggs")), na.rm=T)) %>% 
  left_join(treatments, by=c("plotid","plot","subplot")) %>%
  mutate_if(is.character, as.factor) 

ph20 <- ph20 %>% 
  complete(nesting(plantid, plotid, plot, subplot, snow, water),nesting(julian, date), fill=list(open=0,buds=0)) %>% #add zeros to weeks the plant was not counted
  mutate(flowering = open + buds > 0,
         has_egg = eggs > 0)

# leaftraits --------------------------------------------------------------

lt_sheets <- filter(datasheets, name=="2020 Maxfield Leaf Traits")
lt <- bind_rows(lapply(sheet_names(lt_sheets), function(x) read_sheet(lt_sheets, sheet=x))) %>% 
  mutate(plot = as.character(plot),
         plotid = paste0(plot, subplot),
         year = year(date_collected),
         round = factor(ifelse(yday(date_collected) > 200,1,2)),
         trichome_density = trichomes / leaf_area_cm2,
         sla = leaf_area_cm2 / dry_weight_g,
         sla_wet = leaf_area_cm2 / wet_weight_g,
         water_content = (wet_weight_g-dry_weight_g)/wet_weight_g) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor)

lt.subplotround <- lt %>% 
  group_by(year, round, plot, subplot, plotid, water, water4, snow) %>% 
  summarize_if(is.numeric, mean, na.rm=T)

#average by plant, then subplot
lt.means <- lt %>% mutate_at(c("plotid","plant"), as.character) %>% 
  group_by(plotid, plant) %>% summarize_if(is.numeric, mean, na.rm=T) %>% 
  group_by(plotid) %>% summarize_if(is.numeric, mean, na.rm=T) %>% 
  left_join(treatments, by="plotid")

# floraltraits ------------------------------------------------------------------

mt <- read_sheet(filter(datasheets, name=="2020 Maxfield Floral Traits"), sheet="Morphology") %>% 
  bind_rows(read_sheet(filter(datasheets, name=="2020 Maxfield Floral Traits"), sheet="Morphology1819")) %>% 
  mutate(plant = as.character(plant), plantid = paste0(plotid, plant),
         year = factor(ifelse(is.na(year), year(date), year))) %>% 
  left_join(treatments, by="plotid") %>%
  mutate_if(is.character, as.factor)

#Nectar
nt <- read_sheet(filter(datasheets, name=="2020 Maxfield Floral Traits"), sheet="Nectar") %>% 
  bind_rows(read_sheet(filter(datasheets, name=="2020 Maxfield Floral Traits"), sheet="Nectar1819") %>% mutate(plant=as.list(plant))) %>% 
  mutate(plant = as.character(plant),  plantid = paste0(plotid, plant), 
         year = factor(ifelse(is.na(year), year(date), year))) %>% 
  left_join(treatments, by="plotid") %>%
  mutate_if(is.character, as.factor)

#Merge nectar and morphology, 
#add avg soil moisture and snow melt date, and
#pick traits for analysis
mt <- bind_rows(mt, nt)  %>% 
  left_join(sm.subplotyear) %>% 
  left_join(meltdates)

#average by plant and year, then by subplot
mt.subplot <- mt %>% mutate_at(c("plotid","plant"), as.character) %>% 
  group_by(year, water, snow, plot, plotid, plant) %>% summarize_if(is.numeric, mean, na.rm=T) %>% 
  group_by(year, water, snow, plot, plotid) %>% summarize_if(is.numeric, mean, na.rm=T)

#average by plant and year
mt.plantyr <- mt %>% mutate_at(c("plotid","plant"), as.character) %>% 
  group_by(year, water, snow, plot, plotid, plant) %>% summarize_if(is.numeric, mean, na.rm=T) %>% ungroup

# floralvolatiles ---------------------------------------------------------------

vt <- read_sheet(filter(datasheets, name=="2020 Maxfield Floral Volatiles"), sheet="2020") %>% 
  mutate(plotid = paste0(plot, subplot),
         plant = as.character(plant),
         plantid = paste0(plotid, plant)) %>% 
  left_join(treatments, by="plotid") %>%
  mutate_if(is.character, as.factor)

# seeds -------------------------------------------------------------------

sds18 <- read_sheet(filter(datasheets, name=="2020 Maxfield Seeds"), sheet="2018") %>% 
  mutate(plot = as.character(plot), 
         subplot = toupper(subplot),
         plotid = paste0(plot, subplot),
         plantid = paste0(plotid, plant)) %>% 
  left_join(treatments, by=c("plotid","plot","subplot")) %>%
  mutate_if(is.character, as.factor) 

# export ------------------------------------------------------------------

remove(wx)
save.image("data/maxfield_data.rda")

alldata <- list("treatments"=treatments,
                "meltdates"=meltdates,
                "soil_moisture"=sm,
                "census_2018"=cen18,
                "census_2019"=cen19,
                "census_2020"=cen20,
                "phenology_2018"=ph18,
                "phenology_2019"=ph19,
                "phenology_2020"=ph20,
                "leaf_traits"=lt,
                "floral_traits"=mt,
                "floral_volatiles"=vt,
                "seeds_2018"=sds18)

purrr::walk(names(alldata), ~write_tsv(alldata[[.]], paste0("data/",., ".tsv")))
