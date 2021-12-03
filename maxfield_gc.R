library(tidyverse)
library(reshape2)
library(lubridate)
library(vegan)
library(ggvegan)

setwd("~/MyDocs/MEGA/UCI/Schiedea/Analysis/scent/rmbl/Maxfield")
source("../read_shimadzu.R")

# Read chromatograms ------------------------------------------------------

#reading in the Shimadzu search output is slow, skip it
# setwd("./data")
# datafiles <- list.files(pattern=".txt")
# #safemadzu <- safely(read.shimadzu)
# #shimadzu.data <- map(datafiles, safemadzu)
# #map(shimadzu.data, "error") %>% map_lgl(~!is.null(.)) %>% sum #tally import errors
# maxf.data <-  datafiles %>% set_names() %>% map(read.shimadzu) %>% bind_rows(.id="batch")
# setwd("..")
# #save(maxf.data, file="maxfield210827.Rdata")

load("maxfield210827.Rdata")
maxf.data <- maxf.data %>% filter(batch != "ipo2015.2016.txt")
sort(unique(maxf.data$batch))
#TODO need to filter these a lot better here or later with metadata

maxf.all <- dcast(maxf.data, Filename~Name, sum, value.var="Area")
rownames(maxf.all) <- maxf.all[,1]
maxf.all[,1] <- NULL
maxf.cut <- maxf.all[,colSums(maxf.all)>5e8]#arbitrary cutoff

# k-means  --------------------------------------------------------------------
k <- 40
set.seed(1)
km <- kmeans(decostand(maxf.cut, method="log"), k, nstart=3)
#save(km, file="km30.Rdata")

maxf.km <- tibble(FileName=rownames(maxf.all)) %>% 
  mutate(rowSum = rowSums(maxf.all),
         Project = str_extract(FileName, "Blank|[aA]ir|Aeven|Conditioning|Corydalis|PH|SE|SP|Sterile|Yeast|Veg|Ambient|GNA|SC|Lupinus") %>% replace_na("sample"),
         Type = fct_collapse(Project, blank="Blank", air=c("air","Air"), other_level = "sample"),
         nameBlank = Type=="blank",
         runYear = str_extract(FileName, "2018|2019|2020|2021") %>% replace_na("2018") %>% factor,
         Cluster = km$cluster) %>% # Figure out which k-means clusters are the blanks
  mutate(kBlank = Cluster %in% (count(., nameBlank, Cluster) %>% filter(nameBlank, n>2) %>% pull(Cluster)),
         Mixup = nameBlank != kBlank)

with(maxf.km, table(kBlank, nameBlank))

#  Filter Maxfield samples from inventory -----------------------------------------------------------
# Moved this "Maxfield runs" section from markes_sequence.R, code above that section outputs this Rdata:
load("../Inventory/markes_sequence.rda")

#Try to get all 2021 data to annotate, not just maxfield samples
#811-8132021 are for picking out 2021 Lupinus argenteus samples
#There are some 2021 GNA, VF runs in Test instrument runs too, but adding that folder would break pre-2021 verdict order
Maxfield <- str_detect(sequ.summary$FullName, "Maxfield|07262018|June28_07292018|Ipomopsis 2020|Corydalis|8112021|8122021|8132021")
maxfield.batchids <- sequ.summary %>% filter(Maxfield) %>% select(id) %>% unique() %>% na.omit() 
maxfgc <- sequ.summary %>% filter(id %in% maxfield.batchids$id | Maxfield) #get entire batch if it had a sample that matches

maxfgc <- maxfgc %>% left_join(maxf.km %>% select(FileName, nameBlank, Mixup, kBlank, Cluster)) %>% 
  mutate(verdict="", sample="", index=row_number()) %>% 
  left_join(maxf.data %>% rename(FileName=Filename) %>% group_by(FileName,batch) %>% tally(name="n_peaks")) %>% #get batch file names from Shimadzu output
  select(c("index", "sequence.start", "batch", "Desorb.Start.Time", "CreationTime", "eithertime", "status", 
           "Tube", "markes_n", "GC_n", "either_n", "markes_GC", "create_desorb", "desorb.Start.diff", 
           "Mixup", "nameBlank", "kBlank", "Cluster", "n_peaks", "verdict", "FileName", "sample", "user", "FullName", "id"))
write_csv(maxfgc, "maxfield_all210827_updated.csv")

#Output only Maxfield samples to split filenames into parts
data_inventory <- "1X8oo7qZlo1p6MVl_CBeBe6CUTHEAcd-FWQzfHud3Qws" #"RMBL GC-MS Data Inventory"
read_sheet(data_inventory, sheet="maxfield_all210827annot", na="NA") %>% 
  filter(user=="D Campbell : Maxfield") %>% 
  mutate(sample = na_if(sample,"") %>% coalesce(FileName) %>% str_remove(".qgd")) %>% 
  select(index, batch, Desorb.Start.Time, verdict, sample) %>% 
  separate(sample, into=paste0("file",1:9), remove=F) #%>% 
#range_write(ss=data_inventory, sheet="maxfield_meta",range="A:N")

#get hand-split filename metadata
#TODO work more on this metadata, some just have vial numbers
maxfmeta <- read_sheet(data_inventory, sheet="maxfield_meta", guess_max=2000, col_types="c") %>% 
  select(index, type:vial) %>% as.data.frame %>% 
  drop_na(type) %>% #ditch the gunk for now
  distinct(index, .keep_all = T) %>% # TODO investigate these dupes with two Shimadzu batches
  left_join(maxfgc %>% distinct(index, .keep_all = T) %>% mutate(index=as.character(index))) %>%
  mutate(plotid=ifelse(type=="floral", str_sub(plantid, 1,2),NA),
         plant= ifelse(type=="floral", str_sub(plantid, 3),NA),
         sampledate=ymd(ifelse(year(sequence.start)==2018, paste0("2018",sampledate), sampledate)))

#subset data to only Maxfield samples with an annotated type (floral/ambient) from maxfmeta
maxf <- maxf.all[maxfmeta$FileName,]
rownames(maxf) <- rownames(maxfmeta) <- maxfmeta$index
save(maxf, maxfmeta, file="~/MyDocs/MEGA/UCI/RMBL 2020/maxfield/data/volatiles/maxfield_volatiles.rda")

# NMDS of blanks and all samples --------------------------------------------------------------------

nmds.maxf <- metaMDS(sqrt(maxf.cut), dist="bray", autotransform = FALSE, try=1, trymax=1)
nmds.points <- fortify(nmds.maxf) %>% as_tibble() %>% 
  filter(Score=="sites") %>% left_join(maxf.km, by=c("Label"="FileName"))

ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, color=log(rowSum))) + 
  geom_point() + scale_color_viridis_c() + theme_dark()

ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, color=Cluster, shape=Type)) + geom_point() +
  scale_color_gradientn(colors=turbo(k)) + theme_dark()

ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, color=nameBlank, shape=Type)) + geom_point()

ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, label=Cluster, color=Type)) + geom_text(size=3)

ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, color=runYear, alpha=nameBlank)) + geom_point() +
  scale_alpha_manual(values=c(0.2,1))

# CAP of blanks vs all samples ---------------------------------------------------------------------

maxf.cap <- capscale(maxf.cut ~ kBlank * runYear, distance="bray", metaMDSdist = F, data=maxf.km)
maxf.cap.points <-  fortify(maxf.cap) %>% as_tibble() %>% 
  filter(Score=="sites") %>% left_join(maxf.km, by=c("Label"="FileName"))
ggplot(maxf.cap.points, aes(x=CAP1, y=CAP2, color=runYear, alpha=nameBlank)) + geom_point() +
  scale_alpha_manual(values=c(0.2,1))
