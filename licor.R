# Parse LICOR files for Maxfield Meadow Ipomopsis experiment
# adapted from https://www.ericrscott.com/post/li-cor-wrangling/

library(tidyverse)
library(lubridate)

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

write_tsv(licor, file="data/licor.tsv")

# group by plantid and tally for each year
#note that this averages across dates if a plant was remeasured, including if loose-formatted SampleID is differenct
licor.plantid <- licor %>% group_by(untagged, plantid, plot, subplot, plotid, plant, year) %>% add_tally() %>% 
  summarize(across(where(is.numeric), mean, na.rm=T), .groups="keep") 

# tally for just the tagged plants in wide format
licor.tally <- licor.plantid %>% select(group_cols(),n) %>% 
  filter(!untagged) %>% ungroup %>% select(-untagged) %>% 
  mutate(year=str_sub(year,3)) %>% arrange(year) %>% 
  pivot_wider(names_from=year, names_prefix="licor", values_from=n)

read_csv("data/megatally.csv") %>% mutate(plot=as.character(plot)) %>% 
  full_join(licor.tally) %>% write_tsv("data/megatally_licor.tsv", na="")

# plot all the LICOR variables by date
licor_date_pal <- sample(rainbow(nlevels(factor(licor$date))))

licor %>% pivot_longer(Photo:StableF) %>% mutate(index=row_number()) %>% filter(!name %in% c("Area","BLC_1","BLCond","StmRat")) %>% 
  ggplot(aes(x=index, y= value, color=factor(date), shape=year)) + facet_wrap(vars(name), scales="free_y", ncol=4) + geom_point() + scale_color_manual(values=licor_date_pal, guide=F) + theme_dark() + labs(x="Measurement", y="", shape="Year")

licor %>% ggplot(aes(x=Cond, y= Photo, color=factor(date), linetype=year, shape=year)) + geom_point() + geom_smooth(method="lm", se=F) + scale_linetype_manual(values=c(3,1,4)) + scale_color_manual(values=licor_date_pal, guide=F) + theme_dark()

licor %>% ggplot(aes(y=Photo/Cond, x=factor(date), color=year)) + geom_boxplot(size=1.5) + 
  scale_y_continuous(limits=c(0,350)) + scale_color_manual(values=year_pal) + scale_x_discrete(guide = guide_axis(angle = 90)) + theme_minimal()
licor %>% left_join(treatments)%>% ggplot(aes(y=Photo/Cond, x=factor(date), color=water)) + 
  geom_boxplot(position=position_dodge2(preserve="single")) + 
  scale_y_continuous(limits=c(0,350)) + scale_color_manual(values=water_pal) + scale_x_discrete(guide = guide_axis(angle = 90)) + theme_minimal()
