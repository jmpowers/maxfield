library(tidyverse)
library(reshape2)
library(vegan)

# Read chromatograms ------------------------------------------------------

#reading in the Shimadzu search output is slow, skip it
# setwd("~/MyDocs/MEGA/UCI/Schiedea/Analysis/scent/rmbl/Maxfield")
# source("../read_shimadzu.R")
# setwd("./data")
# datafiles <- list.files(pattern=".txt")
# #safemadzu <- safely(read.shimadzu)
# #shimadzu.data <- map(datafiles, safemadzu)
# #map(shimadzu.data, "error") %>% map_lgl(~!is.null(.)) %>% sum #tally import errors
# maxf.data <-  datafiles %>% set_names() %>% map(read.shimadzu) %>% bind_rows(.id="batch")
# setwd("..")
# #save(maxf.data, file="maxfield210827.Rdata") #copied this to maxfield folder:

load("data/volatiles/maxfield210827.Rdata")
maxf.data <- maxf.data %>% filter(!batch %in% c("ipo2015.2016.txt", #old files
                                  "heather_maxfield190813_new.txt")) #subset of heather_maxfield190813_new001.txt
# Duplicate filenames in Shimadzu output:
# 2 unique chomatograms for Blank1_762021_01.qgd
#These two will get added together in maxf.all
maxf.data %>% count(Filename,batch) %>% group_by(Filename) %>% filter(n()>1)

maxf.all <- dcast(maxf.data, Filename~Name, sum, value.var="Area")
rownames(maxf.all) <- maxf.all[,1]
maxf.all[,1] <- NULL
maxf.cut <- maxf.all[,colSums(maxf.all)>5e8]#arbitrary cutoff

# k-means  --------------------------------------------------------------------
k <- 40
set.seed(1)
km <- kmeans(decostand(maxf.cut, method="log"), k, nstart=3)

maxf.km <- tibble(FileName=rownames(maxf.all)) %>% 
  mutate(rowSum = rowSums(maxf.all),
         Project = str_extract(FileName, "Blank|[aA]ir|Aeven|Conditioning|Corydalis|PH|SE|SP|Sterile|Yeast|Veg|Ambient|GNA|SC|Lupinus") %>% replace_na("sample"),
         Type = fct_collapse(Project, blank="Blank", air=c("air","Air"), other_level = "sample"),
         nameBlank = Type=="blank",
         runYear = str_extract(FileName, "2018|2019|2020|2021") %>% replace_na("2018") %>% factor,
         Cluster = km$cluster) %>% # Figure out which k-means clusters are the blanks
  mutate(kBlank = Cluster %in% (count(., nameBlank, Cluster) %>% filter(nameBlank, n>2) %>% pull(Cluster)),
         Mixup = nameBlank != kBlank) #mixup if filename=Blank but doesn't cluster with blanks, or vice versa

with(maxf.km, table(kBlank, nameBlank))

#  Filter Maxfield samples from inventory -----------------------------------------------------------
# Moved this "Maxfield runs" section from markes_sequence.R, code above that section outputs this Rdata:
#fixed: the new "../Inventory/markes_sequence.rda" is likely causing the discrepancies with 210827_updated - 2022 data processed with different fuzzy join parameters. original +0 offset, +-16 min tolerance, now +18 +-15
#trying the old markes_sequence.rda (2021-12-02), also see maxfield/markes_sequence.R (same date)
#TODO make sure this old version can be generated with that old code
load("data/volatiles/markes_sequence.rda") #loads sequ, sequ.file, sequ.summary, qgdfiles

#Try to get all 2021 data to annotate, not just maxfield samples
#811-8132021 are for picking out 2021 Lupinus argenteus samples
#There are some 2021 GNA, VF runs in Test instrument runs too, but adding that folder would break pre-2021 verdict order
Maxfield <- str_detect(sequ.summary$FullName, "Maxfield|07262018|June28_07292018|Ipomopsis 2020|Corydalis|8112021|8122021|8132021")
maxfield.batchids <- sequ.summary %>% filter(Maxfield) %>% select(id) %>% unique() %>% na.omit() 
maxfgc.sequ <- sequ.summary %>% filter(id %in% maxfield.batchids$id | Maxfield) #get entire batch if it had a sample that matches

# TODO Duplicate filenames in Markes/GC sequence: 
#   2x2 unique chromatograms for both Blank1_070621_1.qgd, Corydalis_200727_12_redo_762021_08.qgd (OK) - gives join warning message below
#   8x2 mismatch fuzzy matchings with one file, two close desorb times (can remove later with verdict=mismatch?)
maxfgc.sequ %>% drop_na(FileName) %>% group_by(FileName) %>% filter(n()>1) %>% arrange(FileName)

#join Markes/GC sequence, k-means, and batch file names
#many-tomany warning is from duplicate filenames in maxf.data (2 blanks), and duplicate filenames in Markes/GC sequence (mismatches)
maxfgc <- maxfgc.sequ %>% 
  left_join(maxf.km %>% select(FileName, nameBlank, Mixup, kBlank, Cluster)) %>% 
  mutate(verdict="", sample="", index=row_number()) %>% 
  left_join(maxf.data %>% rename(FileName=Filename) %>% group_by(FileName,batch) %>% tally(name="n_peaks")) %>% #get batch file names from Shimadzu output
  select(c("index", "sequence.start", "batch", "Desorb.Start.Time", "CreationTime", "eithertime", "status", 
           "Tube", "markes_n", "GC_n", "either_n", "markes_GC", "create_desorb", "desorb.Start.diff", 
           "Mixup", "nameBlank", "kBlank", "Cluster", "n_peaks", "verdict", "FileName", "sample", "user", "FullName", "id")) %>% 
  write_csv("data/volatiles/maxfield_all210827_lessdupes.csv")#used to output to the rmbl/Maxfield folder. Original from project folder backed up in RMBL 2020/bak

#Output only Maxfield samples to split filenames into parts
# library(googlesheets4)
# gs4_auth(email = T)
# data_inventory <- "1X8oo7qZlo1p6MVl_CBeBe6CUTHEAcd-FWQzfHud3Qws" #"RMBL GC-MS Data Inventory"
# read_sheet(data_inventory, sheet="maxfield_all210827annot", na="NA") %>% 
#   filter(user=="D Campbell : Maxfield") %>% 
#   mutate(sample = na_if(sample,"") %>% coalesce(FileName) %>% str_remove(".qgd")) %>% #renames skips
#   select(index, batch, Desorb.Start.Time, verdict, sample) %>% 
#   separate(sample, into=paste0("file",1:9), remove=F) #%>% 
# #range_write(ss=data_inventory, sheet="maxfield_meta",range="A:N") #output to sheet - do not uncomment!

#get hand-split filename metadata
#join with maxfgc is very fragile (by index) so can't ever change maxfgc order
#some of these only have vial numbers, but it doesn't matter since they were sampled in 2021 OTC experiment
#index is from maxfield_all210827annot, sample accounts for skips
#FIXED maxfgc and maxfield_meta gsheet index columns do not line up so this join is broken - 
# looks like it broke with the maxfield_all210827_updated.csv (versus maxfield_all210827.csv)
#read_sheet(data_inventory, sheet="maxfield_meta", guess_max=2000, col_types="c") %>% 
maxfmeta <- read_tsv("data/volatiles/RMBL GC-MS Data Inventory - maxfield_meta.tsv") %>% 
  filter(batch!="heather_maxfield190813_new.txt") %>% #drop duplicate batch
  select(index, sample, type:vial) %>% as.data.frame %>% #added sample (renamed)
  drop_na(type) %>% # exclude files that don't have a manual "type" entry - blanks, leaks, skips, OTC experiment
  #distinct(index, .keep_all = T) %>% # FIXED dupes with two Shimadzu batches
  left_join(maxfgc %>% distinct(index, .keep_all = T) %>% select(-sample)) %>% #took off sample (empty) 
  mutate(plotid=ifelse(type=="floral", str_sub(plantid, 1,2),NA),
         plant= ifelse(type=="floral", str_sub(plantid, 3),NA),
         sampledate=ymd(ifelse(year(sequence.start)==2018, paste0("2018",sampledate), sampledate)),
         index=as.character(index)) %>% 
  write_csv("data/volatiles/maxfmeta_lessdupes.csv")

# #testing to see if duplicates can be removed without breaking index orders of maxfgc and maxfield_meta
# read_tsv("data/volatiles/RMBL GC-MS Data Inventory - maxfield_meta.tsv") %>% filter(index %in% index[duplicated(index)]) %>% with(table(sample, batch))
# read_tsv("data/volatiles/RMBL GC-MS Data Inventory - maxfield_meta.tsv") %>% filter(batch=="heather_maxfield190813_new.txt")
# maxfgc %>% filter(index %in% index[duplicated(index)])
# maxfgc %>% filter(batch=="heather_maxfield190813_new.txt")

# #testing differences between these versions, decided that oldseq will work
# origpath <- "~/MyDocs/MEGA/UCI/Schiedea/Analysis/scent/rmbl/Maxfield/"
# maxfgcnew <- read_csv(paste0(origpath, "maxfield_all210827_updated.csv"))# different from the other two
# maxfgcmid <- read_csv("data/volatiles/maxfield_all210827_updated.csv")#almost identical to old
# maxfgcold <- read_csv(paste0(origpath, "maxfield_all210827.csv"))
# plot(maxfgcmid$Cluster,maxfgcold$Cluster, pch=19, col=alpha("black",0.1))#kmeans clusters are different though
# plot(match(maxfgcnew$Desorb.Start.Time, maxfgcold$Desorb.Start.Time))
# plot(match(maxfgcold$Desorb.Start.Time, maxfgcnew$Desorb.Start.Time))
# plot(maxfgcnew$Desorb.Start.Time)
# plot(maxfgcold$Desorb.Start.Time)
# setdiff(maxfgcold$FileName, maxfgcnew$FileName) #some old filenames deleted - 2018 and 2019
# setdiff(maxfgcnew$FileName, maxfgcold$FileName) #no new filenames
# maxfgcold$deleted <- maxfgcold$FileName %in% setdiff(maxfgcold$FileName, maxfgcnew$FileName) 
# plot(maxfgcold$deleted)
# oldnew <- full_join(maxfgcold, 
#                     maxfgcnew %>% rename(newindex=index) %>% 
#                       select(-CreationTime, -eithertime, -create_desorb, -Mixup, -kBlank, -Cluster))
# sum(is.na(oldnew$newindex))/nrow(oldnew)
# plot(oldnew$index, oldnew$newindex)
# maxfgcoldseq <- read_csv("data/volatiles/maxfield_all210827_oldseq.csv")#newly created from maxfield/data/volatiles/markes_sequence.rda (2021-12-02), identical to mid except for quotes around "NA NA"

#subset wide data to only Maxfield samples with an annotated type (floral/ambient) from maxfmeta
maxf <- maxf.all[maxfmeta$FileName,]
rownames(maxf) <- rownames(maxfmeta) <- maxfmeta$index
save(maxf, maxfmeta, file="data/volatiles/maxfield_volatiles_lessdupes.rda")

# # NMDS of blanks and all samples --------------------------------------------------------------------
# ggplot(maxfgc, aes(x=n_peaks, fill=paste(nameBlank,kBlank))) + geom_histogram() + facet_wrap(vars(year(eithertime)))
# 
# library(vegan)
# library(ggvegan)
# nmds.maxf <- metaMDS(sqrt(maxf.cut), dist="bray", autotransform = FALSE, try=1, trymax=1)
# nmds.points <- fortify(nmds.maxf) %>% as_tibble() %>% 
#   filter(Score=="sites") %>% left_join(maxf.km, by=c("Label"="FileName"))
# 
# ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, color=log(rowSum))) + 
#   geom_point() + scale_color_viridis_c() + theme_dark()
# 
# ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, color=Cluster, shape=Type)) + geom_point() +
#   scale_color_gradientn(colors=turbo(k)) + theme_dark()
# 
# ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, color=nameBlank, shape=Type)) + geom_point()
# 
# ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, label=Cluster, color=Type)) + geom_text(size=3)
# 
# ggplot(nmds.points, aes(x=NMDS1, y=NMDS2, color=runYear, alpha=nameBlank)) + geom_point() +
#   scale_alpha_manual(values=c(0.2,1))
# 
# # CAP of blanks vs all samples ---------------------------------------------------------------------
# 
# maxf.cap <- capscale(maxf.cut ~ kBlank * runYear, distance="bray", metaMDSdist = F, data=maxf.km)
# maxf.cap.points <-  fortify(maxf.cap) %>% as_tibble() %>% 
#   filter(Score=="sites") %>% left_join(maxf.km, by=c("Label"="FileName"))
# ggplot(maxf.cap.points, aes(x=CAP1, y=CAP2, color=runYear, alpha=nameBlank)) + geom_point() +
#   scale_alpha_manual(values=c(0.2,1))
