# Load in all the Maxfield Meadow Ipomopsis data
# Includes treatments, loggers, census, phenology, floral and leaf traits, and seeds

# TODO write metadata  to /metadata for the final data files
# TODO add the following datasets:
#       + floral volatiles 2018, 2019
#       + leaf physiology measures (LICOR) 2018, 2019, 2020

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

# hobo --------------------------------------------------------------------

# The experimental treatments. "water4" keeps all four treatments, "water" lumps control and mock rainout together.
treatments <- read_sheet(filter(datasheets, name=="2020 Maxfield Soil Moisture & Snow"), sheet="treatments") %>%
  mutate(
    water4 = factor(recode(water, 
                           "control"="Control", "mock rainout"="Mock rainout", "addition"="Addition", "rainout"="Reduction")),
    water = factor(recode(water4, "Mock rainout"="Control")),
    snow = factor(recode(snow, "early"="Early","normal"="Normal"))) %>% 
  separate(plotid, c("plot","subplot"), sep=1, remove=F)

treatments_map <- treatments
treatments <- select(treatments, !ends_with(c("x","y")))

# Tarps were placed on plots 2, 4 and 5 in the spring to make the snow melt faster.
# HOBO pendant temperature/light loggers (model UA-002-xx?) were placed at the corners of each plot.
# HOBO manual: https://www.onsetcomp.com/files/manual_pdfs/9556-M%20UA-002%20Manual.pdf
# When the snow melts, the temperature rises above freezing and the light intensifies.
# The day of snowmelt can be calculated from the following thresholds.
# DRC also provided a list of the melt_date - this uses 
# HOBO light thresholds of 10000 lux, except for 2 plots in 2020.

melt_threshold_sun  <- 10000 # light units (probably lux = lumen / m2)
melt_threshold_warm <- 1    # degrees Celsius difference (positive or negative) from 0C (inside snow)

hobo <- drive_download(filter(datasheets, name=="maxfield_hobo_data.csv"), overwrite = T)$local_path %>% read_csv() %>% 
  mutate(time=with_tz(time, "Etc/GMT+6"), #the HOBOs do not account for DST, and were likely all set up during MST
         sun = intensity > melt_threshold_sun, 
         warm = abs(temp) > melt_threshold_warm) %>% 
  mutate_at(c("year","plot"), as.factor)

meltdates <- full_join(
  hobo %>% filter(month(time)%in%4:6, sun==TRUE) %>% 
    group_by(year, plot) %>% summarize_at("time", min) %>% 
    mutate(sun_date=yday(time)) %>%  rename(sun_time=time),
  hobo %>% filter(month(time)%in%4:6, warm==TRUE) %>% 
    group_by(year, plot) %>% summarize_at("time", min) %>% 
    mutate(warm_date=yday(time)) %>%  rename(warm_time=time)) %>% 
  # In 2018 the HOBO logger for plot 4 did not record any data - impute data from plot 5's melt_date
  bind_rows(data.frame(year=factor("2018"), plot=factor("4"), sun_date=114, warm_date=114)) %>% 
  left_join(treatments %>% select(plot, snow) %>% distinct) %>% 
  left_join(read_sheet(filter(datasheets, name=="2020 Maxfield Soil Moisture & Snow"), sheet="snowmelt") %>% 
              mutate(year=factor(year), plot=as.factor(plot)) %>% select(-notes),
            snow_new = factor(recode(snow_new, "early"="Early","normal"="Normal"))) %>% 
  mutate(plot=factor(plot),
         melt_time = parse_date_time(paste(year,melt_date,12,0), "y j H M", tz="Etc/GMT+6")) %>%
  # In 2019 the avalanche caused plot 2 to melt after the normal snowmelt plots, so this is recoded as normal
  # This still records plot 5 as an "early" (1 day earlier)
  mutate(snow = factor(ifelse(year=="2019" & plot=="2", "Normal", as.character(snow))))

# Calculate an offset from the mean sun_date in the normal plots
meltdates <- meltdates %>% left_join(meltdates %>% group_by(year,snow) %>% summarize_at("sun_date", mean) %>% 
            filter(snow=="Normal") %>% select(-snow) %>% rename(sun_date_normal_mean = sun_date)) %>% 
  mutate(melt_offset = sun_date - sun_date_normal_mean) 

# Update the treatments with the actual meltdates and updated early/normal codes in each year
treatments <- treatments %>% select(-snow) %>% left_join(meltdates)

# Color palettes for water and snow treatments
water4_pal <- setNames(brewer.pal(9,name="Set1")[c(2,9,7,1)], levels(treatments$water4))
water_pal <- setNames(brewer.pal(9,name="Set1")[c(2,9,1)], levels(treatments$water))
snow_pal <- setNames(brewer.pal(3, name="Dark2")[c(2,1)], levels(treatments$snow))
year_pal <- setNames(brewer.pal(8, name="Set2")[c(2,3,6)], levels(treatments$year))

# weather -----------------------------------------------------------------

# Weather data from billy barr's RMBL station in Gothic, CO, USA provided by the Western Regional Climate Center
# Original (not FPA) data from 2017-2020 in metric units
# Downloaded as xls (csv download broken) from https://wrcc.dri.edu/cgi-bin/rawMAIN.pl?corbil
# Columns renamed from WRCC names according to /metadata/weather_wrcc_metadata.csv

wx <- drive_download(filter(datasheets, name=="billy_rmbl_wrcc_weather_2017_2020.csv"), overwrite = T)$local_path %>% 
  read_tsv(col_types="Dtnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn")
daily_precip <- wx %>% group_by(date) %>% summarize_at("precip_mm", sum) %>% mutate(year=factor(year(date)), julian=yday(date))
daily_temp <- wx %>% group_by(date) %>% summarize_at("av_temp_2m_C", mean, na.rm=T) %>% mutate(year=factor(year(date)), julian=yday(date))

snowcloth <- read_sheet(filter(datasheets, name=="2020 Maxfield Soil Moisture & Snow"), sheet="snowcloth") %>% 
  mutate(plots = as.character(plots), date = force_tz(date, "America/Denver")+hours(12), year=factor(year)) %>% 
  bind_rows(meltdates %>% group_by(year) %>% filter(snow=="Normal") %>% summarize(date=mean(sun_time)) %>% 
              mutate(mean_snow_depth_cm=0, notes="Mean time of snowmelt (sun_time) in normal plots"))

groundcover <- read_sheet(filter(datasheets, name=="2020 Maxfield Soil Moisture & Snow"), sheet="groundcover") %>% 
  mutate(across(starts_with("first"), list(day=yday)), year=year(first_0_cm))

# soil --------------------------------------------------------------------

# Soil moisture readings were taken with a Campbell Scientific HydroSense II https://www.campbellsci.com/hs2
# For each subplot the corners and middle of the plot were tested with 12 cm (or 20 cm?) rods
# The units are volumetric water content expressed as a percentage
# Data from 2018-06-06 and 2018-06-07 did not include all plots and are excluded

sm_sheets <- filter(datasheets, name=="2020 Maxfield Soil Moisture & Snow")
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

breaks_20 <- seq(160, 240, by=20)
sm.subplot20d <- sm.subplot %>% 
  mutate(mo20d = cut(yday(date), breaks_20, labels=paste(breaks_20[1:4]+1, breaks_20[2:5], sep=" - "))) %>% 
  group_by(year, mo20d, plot, subplot, plotid, water, snow) %>% 
  summarize(VWC = mean(VWC, na.rm=T), .groups="drop")

sm.subplotmonth <- sm.subplot %>% mutate(mo = month(date)) %>% group_by(year, mo, plot, subplot, plotid, water, snow) %>% 
  summarize(VWC = mean(VWC, na.rm=T), .groups="drop")

sm.subplotyear <- sm.subplot %>% group_by(year, plot, subplot, plotid, water, snow) %>% 
  summarize(VWC = mean(VWC, na.rm=T), .groups="drop")

#sm.subplotearly <- sm.subplot %>% filter(yday(date)<=183 & yday(date)>=175) %>% group_by(year, plot, subplot, plotid, water, snow) %>% 
#  summarize(VWC = mean(VWC, na.rm=T), .groups="drop")

# census ------------------------------------------------------------------

# Census of all plants at the beginning of each season to mark if they are dead, vegetative, or flowering
# For vegetative plants counted number of leaves on each rosette and longest leaf

# After loading the censuses, check the allowed states of notes/rosettes/flowering for each year and 
# [translate](https://docs.google.com/spreadsheets/d/17qSYLvl4aO3TNIHwz1YW_XdfjX_grSqoy38m4mKYXe8/edit#gid=816285957) 
# each combination into a simplified plant "status". If there is a problem, status is flagged "recheck". 

allowed_nrf <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="status")
status_priority <- c("flowering", "vegetative", "recheck",  "dead_nf", "tagnf") # determines order for deduplicating
translate_nrf_18 <- with(filter(allowed_nrf, year==2018),  setNames(status, nrf))
translate_nrf_19 <- with(filter(allowed_nrf, year==2019),  setNames(status, nrf))
translate_nrf_20 <- with(filter(allowed_nrf, year==2020),  setNames(status, nrf))
status_ok <- c("flowering", "vegetative", "dead_nf")
status_pal <- setNames(brewer.pal(3, name="Set1")[c(1,3,2)], status_ok)

cen18 <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="2018") %>% 
  mutate_at(vars(ends_with(c("longest","leaves"))), ~as.integer(as.character(.))) %>% 
  mutate(plot = as.character(plot),
         plotid = paste0(plot, subplot),
         plant = as.character(plant),
         plantid = paste0(plotid,plant),
         year=factor(year(date)), round=1,
         notes =  fct_explicit_na(notes, "blank"),
         flowering = fct_explicit_na(flowering, "blank"),
         n_rosettes = as.integer(as.character(rosettes)),
         rosettes = recode(as.character(rosettes), `NULL`="blank",.default="one_or_more"),
         status = recode(paste(notes, rosettes, flowering), !!!translate_nrf_18, .default="recheck") %>% fct_relevel(status_priority),) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor)

cen18r2 <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="2018round2") %>% 
  mutate_at(vars(ends_with(c("longest","leaves"))), ~as.integer(as.character(.))) %>% 
  mutate(plant = as.character(plant),
         plantid = paste0(plotid,plant),
         plot = factor(str_sub(plotid,1,1)), subplot = factor(str_sub(plotid,2,2)),
         year="2018",round=2,date=cen18$date[1]+months(2), #TODO find date of second 2018 census
         flowering = recode(notes, "flowering"="y",.missing="n",.default="n"),
         notes =  notes %>% fct_explicit_na("blank") %>% recode("flowering"="blank"),
         n_rosettes = as.integer(as.character(rosettes)),
         rosettes = recode(as.character(rosettes), .missing="blank",.default="one_or_more"),
         status = recode(factor(paste(notes, rosettes, flowering)), !!!translate_nrf_18, .default="recheck") %>% 
           fct_expand(status_priority) %>% fct_relevel(status_priority)) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor)

cen19 <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="2019foranalysis") %>% 
  mutate(plot = as.character(plot),
         plotid = paste0(plot, subplot),
         plant = as.character(plant),
         plantid = paste0(plotid,plant),
         year=factor(year(date)),
         notes = other_notes %>% str_match("tagnf|nf|dead|chewed|eaten")  %>% factor %>% 
           fct_explicit_na("blank") %>% recode(eaten="chewed"),
         flowering = flowering %>% tolower %>% factor %>% fct_explicit_na("blank"),
         n_rosettes = as.integer(as.character(rosettes)),
         rosettes = recode(as.character(rosettes), ID = "indist", `0` = "zero" , `NULL`="blank",.default="one_or_more", ),
         status = recode(paste(notes, rosettes, flowering), !!!translate_nrf_19, .default="recheck") %>% fct_relevel(status_priority),
         r1_longest	= ifelse(rosettes=="zero", NA, r1_longest),
         r1_leaves = ifelse(rosettes=="zero", NA, r1_leaves)) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor)

cen20 <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="2020") %>% 
  select(-ends_with("_19")) %>% 
  mutate(plot = as.character(plot),
         plotid = paste0(plot, subplot),
         plant = as.character(plant),
         plantid = paste0(plotid,plant),
         year=factor(year(date)),
         notes = notes %>% fct_explicit_na("blank"),
         flowering = flowering %>% as.character %>% fct_explicit_na("blank"),
         n_rosettes = as.integer(as.character(rosettes)),
         rosettes = recode(as.character(rosettes), ID = "indist", `0` = "zero" , `NULL`="blank",.default="one_or_more", ),
         status = recode(paste(notes, rosettes, flowering), !!!translate_nrf_20, .default="recheck") %>% fct_relevel(status_priority),
         plotid = paste0(plot, subplot)) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor)

# duplicates --------------------------------------------------------------
# When there is more than one record for a plant tag in the census, which duplicate should be kept? 
#   Currently resolved by the following rules:
#   
# 1. sort by status (keep flowering > vegetative > recheck > dead_nf > tagnf)
# 2. sort by date within status (if statuses are the same, keep the observation with the newest date)
# 3. keep the first row and discard the other duplicates
# 
# This is automatic, but it may be better to resolve duplicates manually in some cases. The duplicates are written out to [spreadsheets](https://docs.google.com/spreadsheets/d/1vuQXBF8fCOjw-mYFLBk7Pm_iRUzSAZ8QmBCSJtt6IhA/edit#gid=367829026)

# Write out the duplicates
censuses <- c("2018","2018r2","2019","2020")
joiners <- c("plot","subplot","plotid","plant","plantid")
cen <- list(cen18, cen18r2, cen19, cen20) %>% 
#  skip writing for speed
#  walk2(.y=censuses, .f = ~ .x %>% group_by(plantid) %>% filter(n()>1) %>% arrange(plantid, status, desc(date)) %>% 
#          write_sheet(ss=filter(datasheets, name=="Maxfield Results"), sheet=paste0(.y,"duplicates"))) %>% 
  map(~ .x %>% arrange(plantid, status, desc(date)) %>% distinct(plantid, .keep_all = TRUE)) %>% # distinct() keeps the first row after sorting, discards the other rows
  map2(paste0("_",str_remove(censuses, "20")), 
       ~ .x %>% rename_with(paste0, !any_of(joiners), .y)) %>% 
  reduce(full_join)
                                                                                                                           
cen.status <- cen %>% select(any_of(joiners)|starts_with("status")) %>% 
  mutate(across(starts_with("status"), 
                . %>% fct_explicit_na("no_record") %>% fct_expand(status_priority) %>% fct_relevel(status_priority)),
         plantid = fct_reorder(factor(plantid), .desc=T, 
                          paste(as.integer(status_18),as.integer(status_19),as.integer(status_20),as.integer(status_18r2))),
         transition_1819 = paste(status_18, status_19, sep=">"),
         transition_1920 = paste(status_19, status_20, sep=">"),
         transition_181920 = paste(status_18, status_19, status_20, sep=">"))

cen.status.long <- cen.status %>% select(!starts_with("transition")) %>% 
  pivot_longer(starts_with("status"),names_to="census",values_to="status") %>% 
  mutate(census = paste0("20",str_remove(census, "status_"))) %>% 
  mutate(across(.fns=factor))

cen.size.long <- cen %>% select(plantid|starts_with("n_rosettes")) %>% 
  pivot_longer(starts_with("n_rosettes"),names_to="census",values_to="n_rosettes") %>%
  mutate(census = paste0("20",str_remove(census, "n_rosettes_")))

# phenology ---------------------------------------------------------------

ph18.raw <- read_sheet(filter(datasheets, name=="2020 Maxfield Phenology"), sheet="2018") %>% 
  mutate(plot = as.character(plot), 
         subplot = toupper(subplot),
         plotid = paste0(plot, subplot),
         plantid = paste0(plotid, plant),
         plant=as.character(plant),
         julian = factor(yday(date)),
         year="2018",
         open = rowSums(select(., starts_with("open")), na.rm=T),
         buds = rowSums(select(., starts_with("buds")), na.rm=T),
         eggs=eggs_total) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor) 

ph18 <- ph18.raw %>% 
  complete(nesting(plantid, plotid, plot, subplot, snow, water),nesting(julian, date), fill=list(open=0,buds=0)) %>% #add zeros to weeks the plant was not counted
  mutate(flowering = open + buds > 0,
         has_egg = eggs > 0)

ph19.raw <- read_sheet(filter(datasheets, name=="2020 Maxfield Phenology"), sheet="2019") %>% 
  mutate(plot = as.character(plot), 
         plotid = paste0(plot, subplot),
         plantid = paste0(plotid, plant),
         plant = as.character(plant),
         julian = factor(yday(date)),
         year="2019",
         open = rowSums(select(., starts_with("open")), na.rm=T),
         buds = rowSums(select(., starts_with("buds")), na.rm=T),
         eggs = rowSums(select(., starts_with("eggs")), na.rm=T)) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor)

ph19 <- ph19.raw %>% 
  complete(nesting(plantid, plotid, plot, subplot, snow, water),nesting(julian, date), fill=list(open=0,buds=0)) %>% #add zeros to weeks the plant was not counted
  mutate(flowering = open + buds > 0,
         has_egg = eggs > 0)

ph20.raw <- read_sheet(filter(datasheets, name=="2020 Maxfield Phenology"), sheet="2020") %>% 
  mutate(plot = as.character(plot), 
         plotid = paste0(plot, subplot),
         plantid = paste0(plotid, plant),
         plant = as.character(plant),
         julian = factor(yday(date)),
         year="2020",
         open = rowSums(select(., starts_with("open")), na.rm=T),
         buds = rowSums(select(., starts_with("buds")), na.rm=T),
         eggs = rowSums(select(., starts_with("eggs")), na.rm=T)) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor) 

ph20 <- ph20.raw %>% 
  complete(nesting(plantid, plotid, plot, subplot, snow, water),nesting(julian, date), fill=list(open=0,buds=0)) %>% #add zeros to weeks the plant was not counted
  mutate(flowering = open + buds > 0,
         has_egg = eggs > 0)

ph <- bind_rows(list("2018" = ph18, "2019" = ph19, "2020" = ph20), .id="year") %>% 
  mutate(julian=as.integer(as.character(julian)))

ph.plantid <- ph %>% left_join(sm.subplotyear) %>% group_by(year, water, water4, snow, snow_new, plot, plotid, plant, plantid) %>% summarize_at(c("open","height_cm", "VWC"), mean, na.rm=T) 

# leaftraits --------------------------------------------------------------

lt_sheets <- filter(datasheets, name=="2020 Maxfield Leaf Traits")
lt <- bind_rows(lapply(sheet_names(lt_sheets), function(x) read_sheet(lt_sheets, sheet=x))) %>% 
  mutate(plot = as.character(plot),
         plotid = paste0(plot, subplot),
         plant=as.character(plant),
         year = factor(year(date_collected)),
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

# floraltraits ------------------------------------------------------------------

mt <- read_sheet(filter(datasheets, name=="2020 Maxfield Floral Traits"), sheet="Morphology") %>% 
  bind_rows(read_sheet(filter(datasheets, name=="2020 Maxfield Floral Traits"), sheet="Morphology1819")) %>% 
  mutate(plant = as.character(plant), plantid = paste0(plotid, plant), 
         plot = factor(str_sub(plotid,1,1)), subplot = factor(str_sub(plotid,2,2)),
         year = factor(ifelse(is.na(year), year(date), year))) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor)

#Nectar
nt <- read_sheet(filter(datasheets, name=="2020 Maxfield Floral Traits"), sheet="Nectar") %>% 
  bind_rows(read_sheet(filter(datasheets, name=="2020 Maxfield Floral Traits"), sheet="Nectar1819") %>% mutate(plant=as.list(plant))) %>% 
  mutate(plant = as.character(plant),  plantid = paste0(plotid, plant), 
         plot = factor(str_sub(plotid,1,1)), subplot = factor(str_sub(plotid,2,2)),
         year = factor(ifelse(is.na(year), year(date), year)),
         nectar_24_h_ul = nectar_48_h_mm * 5 /(2 * 32), #5-uL microcapillary tube 32 mm in length, 2 days
         # units for nectar_conc are degrees Brix = 1 g sucrose / 100 g solution (percentage by mass)
         nectar_density = 1.852e-5 * nectar_conc^2 + 3.665e-3 * nectar_conc + 1, # mg/uL of sucrose solution at 20 C 
         # Density depends on percentage sucrose by mass. polynomial fit to table at:
         # https://www.mt.com/us/en/home/supportive_content/concentration-tables-ana/Sucrose_de_e.html
         nectar_sugar_24_h_mg = nectar_24_h_ul * nectar_density * nectar_conc ) %>% 
  select(-c("nectar_48_h_mm","nectar_density")) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor)

#Merge morphology and nectar
mt <- bind_rows(mt, nt)  %>% 
  left_join(sm.subplotyear) #add avg soil moisture and snow melt date

#average by plant and year
mt.plantyr <- mt %>% mutate_at(c("plotid","plant"), as.character) %>% 
  group_by(year, water, water4, snow, snow_new, plot, plotid, plant) %>% summarize_if(is.numeric, mean, na.rm=T) %>% ungroup %>% 
  mutate_if(is.numeric, ~replace(., is.nan(.), NA))

#average by plant and year, then by subplot
mt.subplot <- mt.plantyr %>%  
  group_by(year, water, water4, snow, snow_new, plot, plotid) %>% summarize_if(is.numeric, mean, na.rm=T)%>% ungroup %>% 
  mutate_if(is.numeric, ~replace(., is.nan(.), NA))

# floralvolatiles ---------------------------------------------------------------

vt <- read_sheet(filter(datasheets, name=="2020 Maxfield Floral Volatiles"), sheet="2020") %>% 
  mutate(plotid = paste0(plot, subplot), plant = as.character(plant),
         plot = as.character(plot),
         plantid = paste0(plotid, plant),
         year=factor(year(date))) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor)

# seeds -------------------------------------------------------------------

sds18 <- read_sheet(filter(datasheets, name=="2020 Maxfield Seeds"), sheet="2018", skip=2) %>% 
  mutate(plot = as.character(plot), 
         subplot = toupper(subplot),
         plotid = paste0(plot, subplot),
         plantid = paste0(plotid, plant),
         plant = as.character(plant),
         year = "2018",
         seeds_mature = seeds - seeds_early,
         fruits_mature = fruits - fruits_early_countable,
         aborts = aborted_nofruit + fruits_aborted,
         flowers_buds = flowers_buds_collected_early + flowers_buds_collected_last) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor) 

counts_19_20 <- c("seeds","fruits","fruits_split","aborts","fruits_fly_no_seeds","fruits_fly_with_seeds","seeds_fly","fruits_caterpillar","fruits_early_uncountable","flowers_buds","flowers_buds_collected_early","flowers_buds_collected_last")

sds19 <- read_sheet(filter(datasheets, name=="2020 Maxfield Seeds"), sheet="2019", skip=1) %>% 
  mutate(plot = as.character(plot), 
         plotid = paste0(plot, subplot),
         plantid = paste0(plotid, plant),
         plant = as.character(plant),
         year = "2019",
         flowers_buds = flowers_buds_collected_early + flowers_buds_collected_last) %>%
  mutate(across(all_of(counts_19_20), replace_na, 0)) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor) 

sds20.date <- read_sheet(filter(datasheets, name=="2020 Maxfield Seeds"), sheet="2020", skip=1) %>% 
  mutate(plot = as.character(plot), 
         plotid = paste0(plot, subplot),
         plantid = paste0(plotid, plant),
         plant = as.character(plant),
         julian = factor(yday(date)),
         year = "2020") %>% 
  mutate(flowers_buds_collected_last = ifelse(yday(date) >= 238, flowers_buds, 0),
         flowers_buds_collected_early = ifelse(yday(date) < 238, flowers_buds, 0)) %>% 
  mutate(across(counts_19_20, replace_na, 0)) %>%
  left_join(treatments) %>%
  mutate_if(is.character, as.factor) 

sds20 <- sds20.date %>% 
  group_by(across(c("plant", "plantid", colnames(treatments)))) %>% 
  summarize(across(all_of(counts_19_20), sum, na.rm=T), dates=n(), .groups="drop")

sds20.date.zeroed <- sds20.date %>% 
  mutate(julian=recode(julian, "215"="216")) %>% 
  filter(date %in% (sds20.date %>% drop_na(date) %>% count(date) %>% filter(n>1) %>% pull(date))) %>% 
  complete(nesting(plantid, plotid, plot, subplot, snow, water), nesting(julian, date), fill=list(seeds=0, fruits=0, fruits_split=0, fruits_fly_no_seeds=0, fruits_fly_with_seeds=0, seeds_fly=0, fruits_caterpillar=0)) %>% 
  mutate(julian=as.integer(as.character(julian)),
         seeds_per_fruit = seeds/fruits,
         fruits_with_seeds = fruits + fruits_split + fruits_fly_with_seeds,
         fruits_nonaborted = fruits_with_seeds + fruits_fly_no_seeds + fruits_caterpillar + fruits_split,
         prop_infested	= (fruits_fly_no_seeds + fruits_fly_with_seeds + fruits_caterpillar) / fruits_nonaborted,
         fruiting = fruits_nonaborted > 0)

# megatally ---------------------------------------------------------------

megatally <- list(
  cen.status,
  cen18 %>% group_by(plotid, plant) %>% tally(name="census18"),
  cen18r2 %>% group_by(plotid, plant) %>% tally(name="census18r2"),
  cen19 %>% group_by(plotid, plant) %>% tally(name="census19"),
  cen20 %>% group_by(plotid, plant) %>% tally(name="census20"),
  lt %>% filter(year=="2018") %>% group_by(plotid, plant) %>% tally(name="leaf18"),
  lt %>% filter(year=="2019") %>% group_by(plotid, plant) %>% tally(name="leaf19"),
  lt %>% filter(year=="2020") %>% group_by(plotid, plant) %>% tally(name="leaf20"),
  mt %>% filter(year=="2018") %>% group_by(plotid, plant) %>% tally(name="floral18"),
  mt %>% filter(year=="2019") %>% group_by(plotid, plant) %>% tally(name="floral19"),
  mt %>% filter(year=="2020") %>% group_by(plotid, plant) %>% tally(name="floral20"),
  vt %>% filter(year=="2020") %>% group_by(plotid, plant) %>% tally(name="volatiles20"),
  ph18.raw %>% group_by(plotid, plant) %>% tally(name="pheno18"),
  ph19.raw %>% group_by(plotid, plant) %>% tally(name="pheno19"),
  ph20.raw %>% group_by(plotid, plant) %>% tally(name="pheno20"),
  sds18 %>% group_by(plotid, plant) %>% tally(name="seeds18"),
  sds19 %>% group_by(plotid, plant) %>% tally(name="seeds19"),
  sds20 %>% group_by(plotid, plant) %>% tally(name="seeds20")
) %>% 
  reduce(full_join) %>% arrange(plotid, plant) %>% 
  mutate(across(where(is.integer), replace_na, 0),
         plantid = factor(paste0(plotid,plant)),
         plot = factor(str_sub(plotid,1,1)), subplot = factor(str_sub(plotid,2,2))) %>% 
  rowwise() %>%  mutate(total_obs = sum(c_across(census18:seeds20), na.rm=T)) %>% ungroup

megatally %>% 
  write_csv("data/megatally.csv")
  #write_sheet(ss=filter(datasheets, name=="Maxfield Results"), sheet="megatally")

megatally.status <- megatally %>% select(!starts_with("census")) %>% 
  mutate(across(starts_with(c("volatiles","seeds")), recode, `0`="no_record", .default="flowering"),
         across(starts_with(c("pheno","floral")), recode, `0`="no_record", `1` = "recheck", .default="flowering"),
         across(starts_with("leaf"), recode, `0`="no_record", .default="vegetative"),
         across(starts_with(c("pheno","floral","volatiles","seeds","leaf")), fct_expand, status_priority),
         across(where(is.factor), fct_explicit_na, "no_record")) %>% 
  arrange(status_18, pheno18, seeds18, status_19, pheno19, seeds19, status_20, pheno20, seeds20, leaf18, leaf19, leaf20, status_18r2, total_obs) %>% 
  mutate(plantid = fct_reorder(plantid, row_number(), .desc=T))

megatally.status.long <- megatally.status %>% select(!starts_with(c("transition","total_obs"))) %>% 
  pivot_longer(status_18:seeds20,names_to="dataset",values_to="status") %>% 
  mutate(across(.fns=factor),
         dataset = fct_relevel(dataset, c("status_18","status_18r2","status_19","status_20","pheno18","seeds18","floral18","pheno19","seeds19","floral19","pheno20","seeds20","floral20","volatiles20","leaf18","leaf19","leaf20")))


sds <- bind_rows(sds18, sds19, sds20) %>%
  left_join(sm.subplotyear) %>% 
  mutate(plantid = paste0(plotid,plant)) %>% 
  left_join(megatally %>% select(plantid|starts_with("floral")) %>% 
              pivot_longer(starts_with("floral"),names_to="year",values_to="flowers_morphology_nectar") %>% 
              mutate(year = paste0("20",str_remove(year,"floral")))) %>% 
  left_join(megatally %>% select(plantid|starts_with("volatiles")) %>% 
              pivot_longer(starts_with("volatiles"),names_to="year",values_to="flowers_volatiles") %>% 
              mutate(year = paste0("20",str_remove(year,"volatiles")))) %>% 
  mutate(
    flowers_volatiles = replace_na(flowers_volatiles, 0), #TODO get volatiles for 2018 & 2019
    flowers_destroyed = flowers_morphology_nectar + flowers_volatiles, #TODO + flowers_color
    seeds_per_fruit	= seeds / fruits,
    fruits_aborted	= aborts + flowers_buds_collected_last,
    seeds_est	= seeds + seeds_fly + (flowers_buds_collected_early + flowers_destroyed + fruits_early_uncountable) * (seeds/(fruits + fruits_aborted + fruits_fly_with_seeds + fruits_fly_no_seeds + fruits_caterpillar)) + fruits_split * seeds_per_fruit, 
    fruits_with_seeds	= fruits + fruits_split + fruits_fly_with_seeds, 
    fruits_nonaborted	= fruits_with_seeds + fruits_fly_no_seeds + fruits_caterpillar + fruits_split,
    flowers_est	= fruits_nonaborted + fruits_aborted + flowers_buds + flowers_destroyed,
    prop_infested	= (fruits_fly_no_seeds + fruits_fly_with_seeds + fruits_caterpillar) / fruits_nonaborted,
    prop_aborted	= fruits_aborted / (fruits_nonaborted + fruits_aborted),
    prop_nonaborted	= fruits_nonaborted / (fruits_nonaborted + fruits_aborted),
    seeds_per_flower = seeds_est / flowers_est)

######## Merge morphology, nectar, phenology, and seeds ####
mnps <- bind_rows(mt, ph, sds)

#average by plant and year
mnps.plantyr <- mnps %>% mutate_at(c("plotid","plant"), as.character) %>% 
  group_by(year, water, water4, snow, snow_new, plot, plotid, plant) %>% summarize_if(is.numeric, mean, na.rm=T) %>% ungroup %>% 
  mutate_if(is.numeric, ~replace(., is.nan(.), NA))

#average by plant and year, then by subplot
mnps.subplot <- mnps.plantyr %>%  
  group_by(year, water, water4, snow, snow_new, plot, plotid) %>% summarize_if(is.numeric, mean, na.rm=T)%>% ungroup %>% 
  mutate_if(is.numeric, ~replace(., is.nan(.), NA))

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
                "seeds_2018"=sds18,
                "seeds_2019"=sds19,
                "seeds_2020"=sds20)

purrr::walk(names(alldata), ~write_tsv(alldata[[.]], paste0("data/",., ".tsv")))
