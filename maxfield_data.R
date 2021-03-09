# Load in all the Maxfield Meadow Ipomopsis data
# Includes treatments, loggers, census, phenology, floral and leaf traits, and seeds

# TODO write metadata  to /metadata for the final data files
# TODO add the following datasets:
#       + floral volatiles 2018, 2019

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

melt_threshold_sun  <- 10000 # light units (probably lux = lumen / m2)
melt_threshold_warm <- 1    # degrees Celsius difference (positive or negative) from 0C (inside snow)

hobo <- read_csv("data/maxfield_hobo_data.csv") %>% 
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
  # In 2018 the HOBO logger for plot 4 did not record any data - impute data from plot 5's sun_date
  bind_rows(data.frame(year=factor("2018"), plot=factor("4"), sun_date=114, warm_date=114)) %>% 
  left_join(treatments %>% select(plot, snow) %>% distinct) %>% 
  # In 2019 the avalanche caused plot 2 to melt after the normal snowmelt plots, so this is recoded as normal
  # This still records plot 5 as an "early" (1 day earlier)
  mutate(snow = factor(ifelse(year=="2019" & plot=="2", "Normal", as.character(snow))),
         plot = factor(plot))

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

#Get average snowmelt dates in plots
normal_meltdate <- meltdates %>% group_by(year) %>% filter(snow=="Normal") %>% summarize(date=mean(sun_time)) %>%
  mutate(day = yday(date), mean_snow_depth_cm=0, notes="Mean time of snowmelt (sun_time) in normal plots")

snowcloth <- read_sheet(filter(datasheets, name=="2020 Maxfield Soil Moisture & Snow"), sheet="snowcloth") %>% 
  mutate(plots = as.character(plots), date = force_tz(date, "America/Denver")+hours(12), year=factor(year)) %>% 
  bind_rows(normal_meltdate)

groundcover <- read_sheet(filter(datasheets, name=="2020 Maxfield Soil Moisture & Snow"), sheet="groundcover") %>% 
  mutate(across(starts_with("first"), list(day=yday)), year=year(first_0_cm))

waterdates <- read_sheet(filter(datasheets, name=="2020 Maxfield Soil Moisture & Snow"), sheet="water_dates")
# weather -----------------------------------------------------------------

# Load output from maxfield_soil_moisture.R that fit soil moisture ~ precipitation model
load("data/daily_precip_est.rda")

# Weather data from billy barr's RMBL station in Gothic, CO, USA provided by the Western Regional Climate Center
# Original (not FPA - unavailable) data from 2009-2021 in metric units
# Downloaded as xls (which results in a tsv) from https://wrcc.dri.edu/cgi-bin/rawMAIN.pl?corbil
# Columns renamed from WRCC names according to /metadata/WRCC/weather_wrcc_metadata.csv
# Flags on each column indicate 0 = raw data, (M)issing data, (E)stimated data, (N)on-measuring time interval, per https://wrcc.dri.edu/cgi-bin/wea_listex.pl?flags

wrcc_metadata <- read_tsv("metadata/WRCC/weather_wrcc_metadata.tsv") %>% mutate(index = row_number(), flag=F)
wrcc_metadata <- bind_rows(wrcc_metadata, wrcc_metadata %>% select(name_metric, col_type, index) %>% filter(! name_metric %in% c("date","time")) %>% 
  mutate(name_metric=paste0("flag_",name_metric), col_type="f",flag=T)) %>% arrange(index)

tenmin_CORBIL <- read_tsv("data/wrcc_corbil.tsv.gz", skip = 4, 
               col_types=paste0(wrcc_metadata$col_type, collapse=""), col_names=wrcc_metadata$name_metric) %>% 
  filter(date >= as.Date("2015-07-01"))#when precip starts to have reasonable values, temp starts 2015-10-01
hourly_CORBIL <- tenmin_CORBIL %>% mutate(hour=hour(time)) %>% group_by(date, hour) %>% summarize(precip_mm = sum(precip_mm), av_temp_2m_C=mean(av_temp_2m_C)) #take the mean of the precip rate per hour that is reported every 10 min
daily_CORBIL <- hourly_CORBIL %>% group_by(date) %>% summarize(precip_mm = sum(precip_mm), av_temp_2m_C=mean(av_temp_2m_C)) %>% 
  mutate(year=factor(year(date)), julian=yday(date)) #sum the hourly precip rates over each 24 hrs

# Other data (no precip) for Kettle Ponds S of Gothic, CO and Judd Falls E of Gothic, CO

# Weather data from Gold Link station in Mt. Crested Butte. Michelle Newcomer and David Brian Rogers. 2020. Gap-filled meteorological data (2011-2020) and modeled potential evapotranspiration data from the KCOMTCRE2 WeatherUnderground weather station, from the East River Watershed, Colorado. ESS-DIVE: Deep Insight for Earth Science Data. doi:10.15485/1734790, version: ess-dive-ce2fd798f9d81f5-20201210T030125380.

daily_KCOMTCRE2 <- read_csv("data/2011_2020_KComCret_1day_ETModeled.csv")

# Weather data from EPA CASTNET GTH161 station, Gothic Research Meadow https://www3.epa.gov/castnet/site_pages/GTH161.html
hourly_GTH161 <- read_csv("data/CASTNET-GTH161-Meteorological - Hourly.csv.gz", 
                          col_types=cols(DATE_TIME=col_datetime("%m/%d/%Y %H:%M:%S"), QA_CODE=col_character())) %>% 
  filter(year(DATE_TIME) < 2012)
daily_GTH161 <- hourly_GTH161 %>% mutate(date = date(DATE_TIME)) %>% group_by(date) %>% 
  summarize(TEMPERATURE=mean(TEMPERATURE, na.rm=T), PRECIPITATION = sum(PRECIPITATION, na.rm=T))

library(rnoaa)
#station_data <- ghcnd_stations()# long download
#maxfield_coords <- data.frame(id = "RMBL", latitude = 38.9495, longitude = -106.9908)
#meteo_nearby_stations(lat_lon_df = maxfield_coords, station_data = station_data,
#                      radius = 10, var = c("PRCP", "TAVG"))#stations within 10 km
# billy barr's cabin in Gothic, CO https://www.ncdc.noaa.gov/cdo-web/datasets/GHCND/stations/GHCND:US1COGN0018/detail
#All units are in tenths of mm, converted to mm here.
 
daily_billy <- meteo_pull_monitors("US1COGN0018", keep_flags = T) %>% 
  mutate(across(contains("flag"), as.factor)) %>% 
  mutate(across(ghcnd_vars <- c("prcp","snow","snwd","wesd","wesf"), ~ as.integer(.x)/10, .names="{.col}_mm"), .keep="unused") %>% 
  mutate(date = date - days(1)) # comparision to other 3 stations shows the precip is recorded the next day

#library(GSODR)
#gsodr.inventory <- get_inventory() %>% #50 km away from Maxfield Meadow - only get Crested Butte and Aspen airports
#  filter(STNID %in% nearest_stations(38.9495, -106.9908,50))

#Combine weather station data
stations <- c(GTH161="EPA_RsrchMdw",KCOMTCRE2="ESSDIVE_GoldLink",CORBIL="WRCC_billy",billy="NOAA_billy")#station="source_location"
first_snow_2020 <- as.POSIXct(c("2020-10-25", "2020-12-31"), tz="UTC") #TODO check this first permanent snow of 2020 (not yet on billy's website)
daily_all <- daily_billy %>% full_join(daily_GTH161) %>% full_join(daily_KCOMTCRE2) %>% full_join(daily_CORBIL) %>% 
  rename(EPA_RsrchMdw=PRECIPITATION, ESSDIVE_GoldLink=Precip, WRCC_billy=precip_mm, NOAA_billy=prcp_mm) %>% 
  mutate(yr = year(date), mo = month(date, label=F),
         ground_covered = factor(c("smmr","wntr"))[1+date %in% do.call(c, map2(c(groundcover$first_snow, first_snow_2020[1]), c(groundcover$first_0_cm, first_snow_2020[2]), ~ seq(date(.x), date(.y), by="day")))]) %>% 
  mutate(WRCC_billy = ifelse(WRCC_billy>50, NA, WRCC_billy)) #cut errors >50mm/day from daily_CORBIL

#correlate daily EPA_RsrchMdw and NOAA_billy to predict what EPA would look like for 2012-2020
EPA_NOAA_daily_model <-  lm(EPA_RsrchMdw ~ NOAA_billy:ground_covered+0, data=daily_all)
#EPA_NOAA_monthly_model <-lm(EPA_RsrchMdw ~ NOAA_billy*season, data=monthly_all %>% mutate(season=factor(c("wntr","smmr"))[1+mo %in% 6:9]))

daily_all <- daily_all %>% 
  mutate(EPA_predicted = is.na(EPA_RsrchMdw) & !is.na(NOAA_billy),
         EPA_NOAA_filled = ifelse(is.na(EPA_RsrchMdw), predict(EPA_NOAA_daily_model, newdata=daily_all), EPA_RsrchMdw),)
stations <- c("EPA_NOAA_filled", stations)
monthly_all <- daily_all %>% group_by(yr, mo) %>% 
  summarize(across(unname(stations), ~ ifelse(sum(is.na(.x))>5, NA, sum(.x, na.rm=T))), .groups="drop") %>% arrange(yr, mo) 

#Add summer precip estimates to treatments from NOAA-filled EPA Research Meadow dataset
# snowmelt date in each plot - start of watering (100% of precip)
# start of watering - last morphology/nectar measurement (100% for controls, 50% for water reduction, or 100% + 1.75 mm/day for water addition)
precip_total <- function(year_start, day_start, day_end) {
  daily_all %>% filter(year(date)==year_start, yday(date) > day_start, yday(date) < day_end) %>% 
    summarize(tot=sum(EPA_NOAA_filled, na.rm=T)) %>% pull(tot)
}
mt.lastday <- data.frame(year=factor(2018:2020), last_day=c(212,219,210))
treatments <- treatments %>% 
  inner_join(waterdates %>% select(-date) %>% pivot_wider(names_from=precip_treatments, values_from=day) %>% 
               mutate(year=factor(year)) %>% rename(water_begin=started, water_end=ended)) %>% 
  inner_join(mt.lastday) %>% # Maxfield Results - timings (mt)
  #inner_join(mt %>% group_by(year) %>% summarize(last_day=yday(max(date, na.rm=T)))) %>% #can't call mt before its made
  #mnps %>% filter(is.na(seeds)) would include phenology measurements
  mutate(precip_prewater_mm = pmap_dbl(list(year, sun_date, water_begin), precip_total),
         precip_postwater_mm = pmap_dbl(list(year, water_begin, last_day), precip_total) * ifelse(water=="Reduction",.5,1) + 
           ifelse(water=="Addition", (14 / 4 / 2) * (last_day - water_begin),0),
         precip_est_mm = precip_prewater_mm + precip_postwater_mm)

#Add summer precip for longer climate record, from snowmelt date to average last morphology/nectar day
groundcover <- groundcover %>% mutate(precip_est_mm = pmap_dbl(list(year, first_0_cm_day, round(mean(mt.lastday$last_day))), precip_total), precip_est_mm = ifelse(year<1990, NA, precip_est_mm))

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
sm.subplot <- sm %>% group_by(year, date, plot, subplot, plotid, water, water4, snow) %>% 
  summarize(VWC = mean(VWC, na.rm=T), .groups="drop") %>% mutate(day=yday(date))

breaks_20 <- seq(160, 240, by=20)
sm.subplot20d <- sm.subplot %>% 
  mutate(mo20d = cut(yday(date), breaks_20, labels=paste(breaks_20[1:4]+1, breaks_20[2:5], sep=" - "))) %>% 
  group_by(year, mo20d, plot, subplot, plotid, water, water4, snow) %>% 
  summarize(VWC = mean(VWC, na.rm=T), .groups="drop")

sm.subplotyear <- sm.subplot %>% group_by(year, plot, subplot, plotid, water, water4, snow) %>% 
  summarize(VWC = mean(VWC, na.rm=T), .groups="drop")

# fixes -------------------------------------------------------------------

fixes <- read_sheet(filter(datasheets, name=="Maxfield Results"), sheet="megatally", col_types="icccccccccccccccicciciiiiiicciiiiiicciiiiiii") %>% 
  separate(corrections, into=c("transition_fixed", "commands"), sep="/ ", fill="right", extra="merge") %>% 
  separate(commands, into=c("command", "separate2"), sep="/ separate", fill="right") %>%
  separate(command, into=c("verb","datasets"), sep=" ", extra="merge") %>% 
  separate(datasets, into=c("datasets","plantid_new"), sep=" \"", fill="right") %>%
  mutate(across(c("transition_fixed", "verb","datasets","plantid_new","separate2"), ~na_if(trimws(.x, whitespace="[ \"]"),""))) %>% 
  pivot_wider(names_from=verb, values_from=datasets) %>% select(!c("NA")) %>%  
  unite("separate", separate, separate2, sep=" ", na.rm=T) %>% 
  mutate(across(c("separate","modify","discard"), replace_na, "")) %>% 
  mutate(transition_fixed = ifelse(is.na(transition_fixed), transition_181920, transition_fixed)) %>% 
  separate(transition_fixed, into=paste("status",18:20,"fixed", sep="_"), sep=" > ", remove=F)

fixes.fixed <- fixes %>% mutate(plantid = ifelse(is.na(plantid_new), plantid, plantid_new))
mergefixes <- fixes.fixed %>% group_by(plantid) %>% tally() %>% filter(n> 1) %>% pull(plantid)
#View(fixes.fixed %>% filter(plantid %in% mergefixes))
#View(fixes %>% select(!starts_with(c("notes","other","census","floral","pheno","leaf","volatiles","seeds"))))

# take a dataset for one year (last two digits of name) and make the changes in fixes
# make sure to run this function before merging with dfs like treatments that use the plotid
# combining the censuses will only work if they are already padded with the insertions to no_records:
# 1. for no_record lines not needing mergefixes, add line to census with new status
# 2. for lines with mergefixes, just change the status and the ID in the "modify" census
fix_dataset <- function(dat, name) { #name must end with two-digit year 20XX
  yr <- str_sub(name, -2) # _year added to the separated plantids - groups all records in a year!
  
  dat_fixed <- dat %>% mutate(plantid = as.character(plantid)) %>% 
    left_join(fixes %>% select(plantid, discard, separate, modify, plantid_new, ends_with("fixed")), by="plantid") %>% 
    filter(!str_detect(discard,  name)) %>%  
    mutate(modified = str_detect(separate, name) | str_detect(modify, name),
           plantid = ifelse(str_detect(separate, name), paste(plantid, yr, sep="_"), plantid),
           plantid = ifelse(str_detect(modify,   name), plantid_new, plantid))
  if(str_detect(name, "census") & name != "census18r2") { # don't have info to modify census18r2
    status_yr       <- paste("status",yr, sep="_")
    status_yr_fixed <- paste("status",yr,"fixed", sep="_") 
    dat_fixed <- dat_fixed %>% mutate(status = !!sym(status_yr_fixed), .keep="unused") %>% 
    mutate(status = ifelse(is.na(status), #fill in the status from the census columns if missing
                           recode(paste(notes, rosettes, flowering), 
                                  !!!translate_nrf[[paste0("20",yr)]], .default="recheck") %>% 
                              fct_expand(status_priority) %>% fct_relevel(status_priority) %>% as.character, 
                           status)) %>% 
    #add on rows for where plant is missing from census & it hasn't just been created by a modify command &  
    # there is something to fill it in status_yr_fixed
      bind_rows(fixes %>% mutate(modified=TRUE) %>% 
                  filter(!!sym(status_yr) == "no_record" & !(plantid %in% dat_fixed$plantid) &
                           !is.na(!!sym(status_yr_fixed)) & !!sym(status_yr_fixed) != "no_record") %>% 
                  select(plantid, discard, separate, modify, plantid_new, ends_with("fixed")) %>% 
                  rename(status=sym(status_yr_fixed))) %>% fill(year)
  }
  if(yr!="20" & name !="census18r2") { #add survival and flowering. can't say "flowering" - a census column name
    status_nextyr_fixed <- paste("status",as.integer(yr)+1,"fixed", sep="_") 
    dat_fixed <- dat_fixed %>% mutate(
      flowered = recode(!!sym(status_nextyr_fixed), "flowering"=1,"vegetative"=0,"dead_nf"=0,.default=as.double(NA)),
      survived =  recode(!!sym(status_nextyr_fixed), "flowering"=1,"vegetative"=1,"dead_nf"=0,.default=as.double(NA)))
  }
  dat_fixed %>% mutate(plot =    str_sub(plantid,1,1), 
                      subplot = str_sub(plantid,2,2),
                      plotid =  str_sub(plantid,1,2),
                      plant =   str_sub(plantid,3)) %>% return()
}

# census ------------------------------------------------------------------

# Census of all plants at the beginning of each season to mark if they are dead, vegetative, or flowering
# For vegetative plants counted number of leaves on each rosette and longest leaf

# After loading the censuses, check the allowed states of notes/rosettes/flowering for each year and 
# [translate](https://docs.google.com/spreadsheets/d/17qSYLvl4aO3TNIHwz1YW_XdfjX_grSqoy38m4mKYXe8/edit#gid=816285957) 
# each combination into a simplified plant "status". If there is a problem, status is flagged "recheck". 

allowed_nrf <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="status")
status_priority <- c("flowering", "vegetative", "recheck",  "dead_nf", "tagnf") # determines order for deduplicating
translate_nrf <- map(set_names(2018:2020), ~ with(filter(allowed_nrf, year==.x), setNames(status, nrf)))
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
         rosettes = recode(as.character(rosettes), `NULL`="blank",.default="one_or_more")) %>% 
  fix_dataset("census18") %>% left_join(treatments) %>% 
  mutate_if(is.character, as.factor)

cen18r2 <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="2018round2") %>% 
  mutate_at(vars(ends_with(c("longest","leaves"))), ~as.integer(as.character(.))) %>% 
  mutate(plant = as.character(plant),
         plantid = paste0(plotid,plant),
         plot = factor(str_sub(plotid,1,1)), subplot = factor(str_sub(plotid,2,2)),
         year="2018",round=2,date=cen18$date[1]+months(2), #TODO find date of second 2018 census in August
         flowering = recode(notes, "flowering"="y",.missing="n",.default="n"),
         notes =  notes %>% fct_explicit_na("blank") %>% recode("flowering"="blank"),
         n_rosettes = as.integer(as.character(rosettes)),
         rosettes = recode(as.character(rosettes), .missing="blank",.default="one_or_more"),
         status = recode(factor(paste(notes, rosettes, flowering)), !!!translate_nrf[["2018"]], .default="recheck") %>% 
           fct_expand(status_priority) %>% fct_relevel(status_priority)) %>% 
  fix_dataset("census18r2") %>% left_join(treatments) %>% 
  mutate_if(is.character, as.factor)

cen19 <- read_sheet(filter(datasheets, name=="2020 Maxfield Rosettes"), sheet="2019foranalysis") %>% drop_na(plot) %>% 
  mutate(plot = as.character(plot),
         plotid = paste0(plot, subplot),
         plant = as.character(plant),
         plantid = paste0(plotid,plant),
         year=factor(year(date)),
         notes = other_notes %>% str_match("tagnf|nf|dead|chewed|eaten|ts")  %>% factor %>% 
           fct_explicit_na("blank") %>% recode(eaten="chewed"),
         flowering = flowering %>% tolower %>% factor %>% fct_explicit_na("blank"),
         n_rosettes = as.integer(as.character(rosettes)),
         rosettes = recode(as.character(rosettes), ID = "indist", `0` = "zero" , `NULL`="blank",.default="one_or_more", ),
         r1_longest	= ifelse(rosettes=="zero", NA, r1_longest),
         r1_leaves = ifelse(rosettes=="zero", NA, r1_leaves)) %>% 
  fix_dataset("census19") %>% left_join(treatments) %>% 
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
         plotid = paste0(plot, subplot)) %>% 
  fix_dataset("census20") %>% left_join(treatments) %>% 
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
  walk2(.y=censuses, .f = ~ .x %>% group_by(plantid) %>% filter(n()>1) %>% arrange(plantid, status, desc(date)) %>% 
          write_sheet(ss=filter(datasheets, name=="Maxfield Results"), sheet=paste0(.y,"duplicates"))) %>% 
  map(~ .x %>% arrange(plantid, status, desc(date)) %>% distinct(plantid, .keep_all = TRUE)) %>% # distinct() keeps the first row after sorting, discards the other rows
  map2(paste0("_",str_remove(censuses, "20")), 
       ~ .x %>% rename_with(paste0, !any_of(joiners), .y)) %>% 
  reduce(full_join)

#if there is no_record, roll over dead_nf and assume flowering plants died
# this should be already taken care of by new transitions, but just in case 
assumedead <- c("dead_nf","flowering")
cen.status <- cen %>% select(any_of(joiners) | (starts_with(c("status","check")) & !contains("fixed")) | contains("notes")) %>% 
  mutate(across(starts_with("status"), 
                ~ .x %>% fct_explicit_na("no_record") %>% fct_expand(status_priority) %>% fct_relevel(status_priority)),
         status_19 = ifelse(status_19 == "no_record" & (status_18 %in% assumedead | status_18r2 %in% assumedead), "dead_nf", as.character(status_19)), 
         status_20 = ifelse(status_20 == "no_record" & status_19 %in% assumedead, "dead_nf", as.character(status_20)), 
         across(starts_with("status"), 
                ~ .x %>% fct_explicit_na("no_record") %>% fct_expand(status_priority) %>% fct_relevel(status_priority)),
         plantid = fct_reorder(factor(plantid), .desc=T, 
                          paste(as.integer(status_18),as.integer(status_19),as.integer(status_20),as.integer(status_18r2))),
         transition_1819 = paste(status_18, status_19, sep=" > "),
         transition_1920 = paste(status_19, status_20, sep=" > "),
         transition_181920 = paste(status_18, status_19, status_20, sep=" > "),
         across(starts_with("status"), .names="flowering_{.col}", 
                ~ recode(.x, "flowering"=1,"vegetative"=0,"dead_nf"=0,.default=as.double(NA))),
         across(starts_with("status"), .names="alive_{.col}", 
                ~ recode(.x, "flowering"=1,"vegetative"=1,"dead_nf"=0,.default=as.double(NA))))
census.status %>% write_tsv(file="data/census_status.tsv") %>% 
  write_sheet(ss=filter(datasheets, name=="Maxfield Results"), sheet="census_status")

cen.status.long <- cen.status %>% select(!starts_with("transition|flowering|alive")) %>% 
  pivot_longer(starts_with("status"),names_to="census",values_to="status") %>% 
  mutate(census = paste0("20",str_remove(census, "status_"))) %>% 
  mutate(across(.fns=factor)) %>% 
  mutate(survived = as.integer(case_when(census == "2018" ~ alive_status_19,
                              census == "2019" ~ alive_status_20))-1,
         flowered = as.integer(case_when(census == "2018" ~ flowering_status_19,
                              census == "2019" ~ flowering_status_20))-1) 

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
  fix_dataset("pheno18") %>% left_join(treatments) %>% 
  mutate_if(is.character, as.factor) 

ph18 <- ph18.raw %>% 
  complete(nesting(plantid, plotid, plot, subplot, snow, water, water4),nesting(julian, date), fill=list(open=0,buds=0)) %>% #add zeros to weeks the plant was not counted
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
  fix_dataset("pheno19") %>% left_join(treatments) %>% 
  mutate_if(is.character, as.factor)

ph19 <- ph19.raw %>% 
  complete(nesting(plantid, plotid, plot, subplot, snow, water, water4),nesting(julian, date), fill=list(open=0,buds=0)) %>% #add zeros to weeks the plant was not counted
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
  fix_dataset("pheno20") %>% left_join(treatments) %>% 
  mutate_if(is.character, as.factor) 

ph20 <- ph20.raw %>% 
  complete(nesting(plantid, plotid, plot, subplot, snow, water, water4),nesting(julian, date), fill=list(open=0,buds=0)) %>% #add zeros to weeks the plant was not counted
  mutate(flowering = open + buds > 0,
         has_egg = eggs > 0)

ph <- bind_rows(list("2018" = ph18, "2019" = ph19, "2020" = ph20), .id="year") %>% 
  mutate(julian=as.integer(as.character(julian))) %>% 
  left_join(sm.subplotyear) 

ph.plantid <- ph %>% group_by(year, water, water4, snow, plot, plotid, plant, plantid) %>% summarize_at(c("open","height_cm", "VWC"), mean, na.rm=T) 

# leaftraits --------------------------------------------------------------

lt_sheets <- filter(datasheets, name=="2020 Maxfield Leaf Traits")
lt <- map_dfr(sheet_names(lt_sheets), ~ read_sheet(lt_sheets, sheet=.x)) %>% 
  mutate(plot = as.character(plot),
         plotid = paste0(plot, subplot),
         plant=as.character(plant),
         plantid=paste0(plotid,plant),
         year = factor(year(date_collected)),
         round = factor(ifelse(yday(date_collected) < 200,1,2)),
         year.round = factor(paste(year, round, sep=".")),
         trichome_density = trichomes / leaf_area_cm2,
         sla = leaf_area_cm2 / dry_weight_g,
         sla_wet = leaf_area_cm2 / wet_weight_g,
         water_content = (wet_weight_g-dry_weight_g)/wet_weight_g) %>% 
  rename(date = date_collected) %>% 
  group_by(year) %>% group_split() %>% map2_dfr(paste0("leaf",18:20), fix_dataset) %>% 
  left_join(treatments) %>% left_join(sm.subplotyear) %>% 
  mutate_if(is.character, as.factor) 

lt.plantyr <- lt %>% 
  group_by(year, round, year.round, plot, subplot, plotid, plant, plantid, water, water4, snow) %>% 
  summarize_if(is.numeric, mean, na.rm=T)

lt.subplotround <- lt.plantyr %>% 
  group_by(year, round, year.round, plot, subplot, plotid, water, water4, snow) %>% 
  summarize_if(is.numeric, mean, na.rm=T)

# licor -------------------------------------------------------------------
# Parser for LICOR files adapted from https://www.ericrscott.com/post/li-cor-wrangling/

licor.path <- "../licor/tsv/"
licor.files <- list.files(licor.path)
licor.dates <- coalesce(mdy(licor.files), ymd(licor.files))

map_chr(paste0(licor.path, licor.files), read_file) %>% 
  str_split("\\$STARTOFDATA\\$") %>% # separate headers from actual data
  map(`[`, 1) %>% #extract just the second element, the actual data
  flatten_chr() %>% paste(collapse="\n\n\n\n") %>%
  write_file(file="data/licor_headers.xml") #not really xml

licor <- map_chr(paste0(licor.path, licor.files), read_file) %>% 
  str_split("\\$STARTOFDATA\\$") %>% # separate headers from actual data
  map(`[`, 2) %>% #extract just the second element, the actual data
  flatten_chr() %>% #converts to a vector
  set_names(licor.dates) %>% #add the dates as a column
  map_dfr( ~ read_tsv(.x, col_types = cols(.default = col_double(), Obs = col_character(), HHMMSS = col_time())), 
           skip = 1, .id="date") %>% #skip the in/out row
  mutate(remark = ifelse(is.na(as.integer(Obs)), Obs, "NA"), #copy those remarks to the remark column
         Obs = ifelse(!is.na(as.integer(Obs)), Obs, NA), #remove remarks from Obs column
         sampleID = str_extract(remark, "[1-6][ABCDabcd]_*.*$"), .after="date") %>% #find Maxfield ID tags - in format 5A123 or 5a_123
  filter(!xor(remark == "NA" , is.na(sampleID))) %>%  #get rid of other remarks  
  fill(sampleID) %>% #fill down the sample ID column
  filter(complete.cases(.)) %>% #get rid of the rest of the remark rows
  select(-remark) %>% #get rid of the remark column
  mutate(plantid = sampleID %>% str_remove("_") %>% str_to_upper %>% str_replace("([A-D])0+","\\1") %>% str_sub(1,5), # standardize plantid. Get rid of underscore, convert to uppercase, delete leading zeros in plant, cap at 5 chars to remove "ACTUAL"
         untagged = str_detect(plantid,"UNT"),
         plantid = ifelse(untagged, paste(plantid,date,Obs, sep="_"), plantid), #keep untagged plants separate by giving them indices of the date and Obs number
         plot = str_sub(plantid,1,1), subplot = str_sub(plantid,2,2), plotid = str_sub(plantid,1,2), plant = str_sub(plantid,3),
         date=as.Date(date), year=factor(year(date)), Obs=as.integer(Obs), .after="date")

#output raw licor data to add VWC and leaf area
licor %>% 
  select(year, date, sampleID, Obs, HHMMSS, Photo, Cond, Ci) %>% group_by(year, date, sampleID) %>% 
  summarize(first_Obs = first(Obs), all_Obs = paste(Obs, collapse=","), n_Obs = n(), 
            HHMMSS = paste(HHMMSS, collapse=","), across(where(is.numeric), mean)) %>% 
  arrange(date, first_Obs) %>%  write_tsv("../licor_sampleID.tsv")

# group by plantid and tally for each year
#note that this averages across dates if a plant was remeasured, including if loose-formatted SampleID is different
licor.plantid <- licor %>% group_by(untagged, plantid, plot, subplot, plotid, plant, year) %>% add_tally() %>% 
  summarize(across(where(is.numeric), mean, na.rm=T), .groups="keep") 

# tally for just the tagged plants in wide format
licor.tally <- licor.plantid %>% select(group_cols(),n) %>% 
  filter(!untagged) %>% ungroup %>% select(-untagged) %>% 
  mutate(year=str_sub(year,3)) %>% arrange(year) %>% 
  pivot_wider(names_from=year, names_prefix="licor", values_from=n)

read_csv("data/megatally.csv") %>% mutate(plot=as.character(plot)) %>% #read in the original megatally
  full_join(licor.tally) %>% write_tsv("data/megatally_licor.tsv", na="") #add licor columns

#All measurements leaf area and VWC available from scans - filled into licor_sampleID.tsv
pt.leafarea <- read_sheet(filter(datasheets, name=="2020 Maxfield Physiology"), sheet="area_VWC", 
                           col_types="dcccccddcdTcddcdddddd")

# read in the corrected WUE survival analysis, with dates and VWC added back in
pt1819 <- read_sheet(filter(datasheets, name=="2020 Maxfield Physiology"), sheet="WUE_1819_cleanedDRC") %>%
  filter(!is.na(photosynthesis)) %>% #some of the rows are duplicated with no licor data
  mutate(plant=as.character(plant), plantid = paste0(plotid,plant), plot = str_sub(plotid,1,1), subplot=str_sub(plotid,2,2), 
         year=factor(year))  

#This sheet has incorrect data for 2018 and 2019 - just use 2020
pt20 <- read_sheet(filter(datasheets, name=="2020 Maxfield Physiology"), sheet="WUE_181920") %>%
  mutate(plot=as.character(plot), year=factor(year), plant=as.character(plant),  subplot=str_sub(plotid,2,2)) %>%  
  filter(year==2020, !is.na(photosynthesis)) #the blank licor rows show the subplots that have missing data

pt.raw <- bind_rows(pt1819, pt20) %>% # lump 2019-08-15 entry with previous day
  mutate(round = as.character(date), round=recode(round, "2019-08-15"="2019-08-14"), round=factor(round))

#add VWC based on subplot mean taken closest in time
pt.sm <- pt.raw %>%   
  left_join(sm.subplot, by = c("year","plot","subplot","plotid"), suffix=c("",".sm")) %>%
  mutate(date_diff = abs(date - date.sm) + (date - date.sm)/100) %>%  #later date breaks ties
  group_by(plotid, date) %>%  filter(date_diff == min(date_diff)) %>% ungroup

pt <- pt.raw %>% left_join(pt.sm %>% select(plantid, date, duplicate, date.sm, date_diff, VWC.sm)) %>%
  rename(VWC.plant = VWC) %>% 
  group_by(year) %>% group_split() %>% map2_dfr(paste0("licor",18:20), fix_dataset) %>% 
  left_join(treatments) %>% left_join(sm.subplotyear) %>% 
  rename(VWC.sm.year = VWC)

# TODO does not lump leaves from the same plant measured in multiple rounds
pt.plantyr <- pt %>% group_by(year, round, water, water4, snow, plot, plotid, plant, plantid) %>% 
  summarize_if(is.numeric, mean, na.rm=T) %>% ungroup %>% 
  mutate_if(is.numeric, ~replace(., is.nan(.), NA))

pt.subplot <- pt.plantyr %>% group_by(year, water, water4, snow, plot, plotid) %>% 
  summarize_if(is.numeric, mean, na.rm=T) %>% ungroup %>% 
  mutate_if(is.numeric, ~replace(., is.nan(.), NA))

# floraltraits ------------------------------------------------------------------

mt <- read_sheet(filter(datasheets, name=="2020 Maxfield Floral Traits"), sheet="Morphology") %>% 
  bind_rows(read_sheet(filter(datasheets, name=="2020 Maxfield Floral Traits"), sheet="Morphology1819")) %>% 
  mutate(plant = as.character(plant), plantid = paste0(plotid, plant), 
         plot = factor(str_sub(plotid,1,1)), subplot = factor(str_sub(plotid,2,2)),
         year = factor(ifelse(is.na(year), year(date), year))) %>% 
  group_by(year) %>% group_split() %>% map2_dfr(paste0("floral",18:20), fix_dataset) %>% 
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
         nectar_sugar_24_h_mg = nectar_24_h_ul * nectar_density * nectar_conc/100 ) %>% 
  select(-c("nectar_48_h_mm","nectar_density")) %>% 
  group_by(year) %>% group_split() %>% map2_dfr(paste0("floral",18:20), fix_dataset) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor)

#Merge morphology and nectar
mt <- bind_rows(mt, nt)  %>% 
  left_join(sm.subplotyear) #add avg soil moisture and snow melt date

#average by plant and year
mt.plantyr <- mt %>% mutate_at(c("plotid","plant"), as.character) %>% 
  group_by(year, water, water4, snow, plot, plotid, plant) %>% summarize_if(is.numeric, mean, na.rm=T) %>% ungroup %>% 
  mutate_if(is.numeric, ~replace(., is.nan(.), NA))

#average by plant and year, then by subplot
mt.subplot <- mt.plantyr %>%  
  group_by(year, water, water4, snow, plot, plotid) %>% summarize_if(is.numeric, mean, na.rm=T)%>% ungroup %>% 
  mutate_if(is.numeric, ~replace(., is.nan(.), NA))

# floralvolatiles ---------------------------------------------------------------

vt <- read_sheet(filter(datasheets, name=="2020 Maxfield Floral Volatiles"), sheet="2020") %>% 
  mutate(plotid = paste0(plot, subplot), plant = as.character(plant),
         plot = as.character(plot),
         plantid = paste0(plotid, plant),
         year=factor(year(date))) %>% 
  group_by(year) %>% group_split() %>% map2_dfr(paste0("volatiles",18:20), fix_dataset) %>% 
  left_join(treatments) %>%
  mutate_if(is.character, as.factor)

# seeds -------------------------------------------------------------------

sds18 <- read_sheet(filter(datasheets, name=="2020 Maxfield Seeds"), sheet="2018", skip=2) %>% 
  filter(!is.na(plot)) %>%  # discard envelope with no plot
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
  fix_dataset("seeds18") %>% left_join(treatments) %>%
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
  fix_dataset("seeds19") %>% left_join(treatments) %>%
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
  mutate(across(all_of(counts_19_20), replace_na, 0)) %>%
  fix_dataset("seeds20") %>% left_join(treatments) %>%
  mutate_if(is.character, as.factor) 

sds20 <- sds20.date %>% 
  group_by(across(c("plant", "plantid", colnames(treatments)))) %>% 
  summarize(across(all_of(counts_19_20), sum, na.rm=T), dates=n(), .groups="drop")

sds20.date.zeroed <- sds20.date %>% 
  mutate(julian=recode(julian, "215"="216")) %>% 
  filter(date %in% (sds20.date %>% drop_na(date) %>% count(date) %>% filter(n>1) %>% pull(date))) %>% 
  complete(nesting(plantid, plotid, plot, subplot, snow, water, water4), nesting(julian, date), fill=list(seeds=0, fruits=0, fruits_split=0, fruits_fly_no_seeds=0, fruits_fly_with_seeds=0, seeds_fly=0, fruits_caterpillar=0)) %>% 
  mutate(julian=as.integer(as.character(julian)),
         seeds_per_fruit = seeds/fruits,
         fruits_with_seeds = fruits + fruits_split + fruits_fly_with_seeds,
         fruits_nonaborted = fruits_with_seeds + fruits_fly_no_seeds + fruits_caterpillar + fruits_split,
         prop_infested	= (fruits_fly_no_seeds + fruits_fly_with_seeds + fruits_caterpillar) / fruits_nonaborted,
         fruiting = fruits_nonaborted > 0)

# megatally ---------------------------------------------------------------

#this will now be the fixed up megatally with all the solutions incorporated
megatally <- list(
  cen.status,
  cen18 %>% group_by(plotid, plant) %>% tally(name="census18"),
  cen18r2 %>% group_by(plotid, plant) %>% tally(name="census18r2"),
  cen19 %>% group_by(plotid, plant) %>% tally(name="census19"),
  cen20 %>% group_by(plotid, plant) %>% tally(name="census20"),
  lt %>% filter(year=="2018") %>% group_by(plotid, plant) %>% tally(name="leaf18"),
  lt %>% filter(year=="2019") %>% group_by(plotid, plant) %>% tally(name="leaf19"),
  lt %>% filter(year=="2020") %>% group_by(plotid, plant) %>% tally(name="leaf20"),
  licor %>% filter(year=="2018") %>% group_by(plotid, plant) %>% tally(name="licor18"),
  licor %>% filter(year=="2019") %>% group_by(plotid, plant) %>% tally(name="licor19"),
  licor %>% filter(year=="2020") %>% group_by(plotid, plant) %>% tally(name="licor20"),
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
  mutate(index=row_number()) %>% 
  mutate(across(census18:seeds20, na_if, 0)) %>% 
  select(index,plot,subplot,plotid,plant,plantid,status_18,status_18r2,status_19,status_20,transition_1819,transition_1920,transition_181920,check_19_20_20,census18,notes_18,other_notes_18,census18r2,other_notes_18r2,pheno18,floral18,seeds18,leaf18,licor18,census19,notes_19,other_notes_19,pheno19,floral19,seeds19,leaf19,licor19,census20,notes_20,other_notes_20,pheno20,floral20,volatiles20,seeds20,leaf20,licor20,total_obs) %>% 
  write_tsv("data/megatally_fixed.tsv", na="") %>% 
  write_sheet(ss=filter(datasheets, name=="Maxfield Results"), sheet="megatally_fixed")

megatally.status <- megatally %>% select(!starts_with(c("census"))) %>% 
  mutate(across(starts_with(c("volatiles","seeds")), recode, `0`="no_record", .default="flowering"),
         across(starts_with(c("pheno","floral")), recode, `0`="no_record", `1` = "recheck", .default="flowering"),
         across(starts_with("leaf"), recode, `0`="no_record", .default="vegetative"),
         across(starts_with("licor"), recode, `0`="no_record", .default="vegetative"),
         across(starts_with(c("pheno","floral","volatiles","seeds","leaf")), fct_expand, status_priority),
         across(where(is.factor), fct_explicit_na, "no_record")) %>% 
  arrange(status_18, pheno18, seeds18, status_19, pheno19, seeds19, status_20, pheno20, seeds20, leaf18, leaf19, leaf20, status_18r2, total_obs) %>% 
  mutate(plantid = fct_reorder(plantid, row_number(), .desc=T))

megatally.status.long <- megatally.status %>% select(!starts_with(c("transition","total_obs","flowering","alive")) & !contains(c("notes","check"))) %>% 
  pivot_longer(status_18:seeds20,names_to="dataset",values_to="status") %>% 
  mutate(across(.fns=factor),
         status = fct_relevel(status, status_priority),
         dataset = fct_relevel(dataset, c("status_18","status_18r2","status_19","status_20","pheno18","seeds18","floral18","pheno19","seeds19","floral19","pheno20","seeds20","floral20","volatiles20","leaf18","licor18","leaf19","licor19","leaf20","licor20")))


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
    flowers_est	= fruits_nonaborted + aborts + flowers_buds + flowers_destroyed,# aborts used to be fruits_aborted, which includes flowers_buds_collected_last, but these are already in flowers_buds
    prop_infested	= (fruits_fly_no_seeds + fruits_fly_with_seeds + fruits_caterpillar) / fruits_nonaborted,
    prop_aborted	= fruits_aborted / (fruits_nonaborted + fruits_aborted),
    prop_nonaborted	= fruits_nonaborted / (fruits_nonaborted + fruits_aborted),
    seeds_per_flower = seeds_est / flowers_est)

######## Merge morphology, nectar, phenology, and seeds ####
mnps <- bind_rows(mt, ph, sds)

#average by plant and year
mnps.plantyr <- mnps %>% mutate_at(c("plotid","plant"), as.character) %>% 
  group_by(year, water, water4, snow, plot, plotid, plant) %>% summarize_if(is.numeric, mean, na.rm=T) %>% ungroup %>% 
  mutate_if(is.numeric, ~replace(., is.nan(.), NA))

#average by plant and year, then by subplot
mnps.subplot <- mnps.plantyr %>%  
  group_by(year, water, water4, snow, plot, plotid) %>% summarize_if(is.numeric, mean, na.rm=T)%>% ungroup %>% 
  mutate_if(is.numeric, ~replace(., is.nan(.), NA)) %>% drop_na(water4) #TODO figure out what entries in ph are causing NAs
write_sheet(mnps.subplot %>% select(!starts_with(c("open_","buds_","nr_","eggs_","fit."))), 
                                    ss=filter(datasheets, name=="Maxfield Results"), sheet="subplot_means")

# timings -----------------------------------------------------------------

timings <- bind_rows(
  snowcloth=snowcloth %>% mutate(plots = recode(plots, "5"="2,5")) %>% # assume that plot 2 was uncovered with plot 5 in 2019
    select(year, day, cloth, plots) %>% drop_na(cloth) %>% group_by(year,plots) %>% pivot_wider(names_from=cloth, values_from=day) %>% rename(begin=added, end=removed),
  meltdates=meltdates %>% group_by(year,snow) %>% summarize_at("sun_date", mean) %>% mutate(sun_date=round(sun_date)) %>% pivot_wider(names_from=snow, values_from=sun_date) %>% rename(begin=Early, end=Normal),
  #summer_precip_timing = summer_precip %>% group_by(year) %>% summarize(begin=min(sun_date), end=max(last_day)),
  waterdates = waterdates %>% select(-date) %>% pivot_wider(names_from=precip_treatments, values_from=day) %>% mutate(year=factor(year)) %>% rename(begin=started, end=ended),
  sm = sm %>% group_by(year) %>% summarize(begin = min(yday(date)), end = max(yday(date))),
  mt = mt %>% drop_na(corolla_length) %>% group_by(year) %>% summarize(begin = min(yday(date), na.rm=T), end = max(yday(date), na.rm=T)),
  nt = nt %>% group_by(year) %>% summarize(begin = min(yday(date), na.rm=T), end = max(yday(date), na.rm=T)),
  ph = ph %>% drop_na(height_cm) %>% group_by(year) %>% summarize(begin = min(julian), end = max(julian)),
  sds= sds20.date %>% group_by(year) %>% summarize(begin = min(yday(date), na.rm=T), end = max(yday(date), na.rm=T)), 
  .id="variable") %>% ungroup

timings %>% mutate(range=paste(begin,end,sep=" - "), .keep="unused") %>% 
  pivot_wider(names_from=year, values_from=range, names_sort=T) %>% 
  write_sheet(ss=filter(datasheets, name=="Maxfield Results"), sheet="timings")

# traitnames --------------------------------------------------------------

traits <- c("corolla_length", "style_length", "corolla_width", "sepal_width", "nectar_24_h_ul", "nectar_conc","nectar_sugar_24_h_mg","height_cm","open","seeds_per_fruit", "seeds_est", "fruits_nonaborted",  "flowers_est", "prop_infested", "prop_aborted", "seeds_per_flower") #exclude anthers
seedtraits <- traits[10:16] 
leaftraits <- c("sla","trichome_density","water_content")
phystraits <- c("photosynthesis", "conductance", "WUE")

traitnames <- set_names(c(
  "Corolla length (mm)", "Style length (mm)", "Corolla width (mm)", "Sepal width (mm)", 
  "Nectar production (\U00B5L/day)", "Nectar conc. (% by mass)", "Nectar sucrose (mg/day)", 
  "Inflorescence height (cm)", 
  "Open flowers", "Seeds per fruit", "Estimated total seeds", "Nonaborted fruits", "Estimated total flowers", 
  "Prop. nonaborted fruits infested","Prop. fruits aborted","Estimated seeds per flower",
  "Specific leaf area (cm\U00B2 g\U207B\U00B9)", "Trichome density (cm\U207B\U00B2)", "Water content", 
  "Photosynthetic rate (\U00B5mol CO\U2082 m\U207B\U00B2 s\U207B\U00B9)", 
  "Stomatal conductance (mol H\U2082O m\U207B\U00B2 s\U207B\U00B9)", 
  "Intrinsic water-use efficiency (\U00B5mol CO\U2082 mol\U207B\U00B9 H\U2082O)"), c(traits, leaftraits, phystraits))

# export ------------------------------------------------------------------

remove(hourly_GTH161, tenmin_CORBIL)
#load("data/maxfield_data.rda")
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
                "licor"=licor,
                "phys_traits"=pt,
                "floral_traits"=mt,
                "floral_volatiles"=vt,
                "seeds_2018"=sds18,
                "seeds_2019"=sds19,
                "seeds_2020"=sds20)

purrr::walk(names(alldata), ~write_tsv(alldata[[.]], paste0("data/",., ".tsv")))
