# Preprocess GCMS chromatograms with xcms and metaMS
# John Powers 2022-04-21

# First section based on xcms (xcms_vignette_GC.R) vignette
library(xcms)
library(tidyverse)

#Parallel processing
register(bpstart(MulticoreParam(3)))

# Maxfield meadow Ipomopsis floral samples, air controls, blanks, and alkane ladders
#TODO get standards to run through this pipeline

## Get the full path to the CDF files
#cdfdir <- "~/MyDocs/MEGA/UCI/Schiedea/GCMS/RMBL GC CDF/maxfield_CDF"
cdfdir <- "/media/john/Bonzai/Backups_Documents/GCMSsolution_200109/Data/Project1_210827/Campbell RMBL Data/Maxfield_CDF"
cdfs <- list.files(cdfdir, pattern = ".CDF", full.names = TRUE, ignore.case = TRUE)

## Create a phenodata data.frame
pd <- tibble(sample_name = str_remove(basename(cdfs), fixed(".CDF"))) %>% 
  mutate(sample_group = ifelse(str_detect(cdfs, "Blank|blank"), "blank", "floral"),
         sample_group = ifelse(str_detect(cdfs, "Air|air|AIR"), "ambient", sample_group),
         sample_group = ifelse(str_detect(cdfs, "ALK|Alk"), "alkanes", sample_group),
         year = str_extract(sample_name, "2021|2019|2018"))

raw_data <- readMSData(files = cdfs, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk")
#save(raw_data, file=paste0(cdfdir, "/raw_data.rda")) 

## Peak detection
peak_params <- MatchedFilterParam(binSize=0.5, fwhm=2, snthresh=2)
xdata <- findChromPeaks(raw_data, param = peak_params)

## Group peaks across samples
pdp <- PeakDensityParam(sampleGroups = xdata$sample_group,
                        minFraction = 0.1, bw = 1)
#TODO try other peak grouping methods, PeakDensity is based on RT only
xdata <- groupChromPeaks(xdata, param = pdp)
#save(xdata, file=paste0(cdfdir, "/xdata_peaks.rda"))

## Gap filling of missing peaks
#run out of memory with more than 2 cores here
xdata <- fillChromPeaks(xdata, param = ChromPeakAreaParam(), BPPARAM = MulticoreParam(2))
#save(xdata, file=paste0(cdfdir, "/xdata_filled.rda"))


# Second section based on metaMS (metaMS_runGC_NIST.R) vignette

library(metaMS)
data(FEMsettings) #default settings for GC-MS in TSQXLS.GC


#GCset <- xdata %>% xcms::split(f=fromFile(xdata)) 
#GCset <- lapply(GCset, as, "xcmsSet")
#doesn't work, use CAMERA directly

#xann <- xsAnnotate(xset)
#xgrouped <- CAMERA::groupFWHM(xann, perfwhm = 1)
#xgrouped.msp <- metaMS::to.msp(xgrouped, file=NULL, settings = metaSetting(TSQXLS.GC, "DBconstruction"))
# running CAMERA combined dataset doesn't give right format for to.msp

## Convert XCMSnExp to list of old xcmsSets (one per file) for metaMS compatibility
xset <- as(xdata, "xcmsSet")
xset.list <- split(xset, factor(sampnames(xset), levels = sampnames(xset)))

## Cluster single m/z peaks with similar retention time into pseudospectra
allSamples <- lapply(xset.list, runCAMERA, chrom = "GC", settings = metaSetting(TSQXLS.GC, "CAMERA")) 

## Convert to MSP text format - takes ~10 min
allSamples.msp <- lapply(allSamples, to.msp, file = NULL, settings = metaSetting(TSQXLS.GC, "DBconstruction"))
#save(allSamples.msp, file=paste0(cdfdir, "/allSamples_msp.rda"))

## Match to dummy library of 3 standards, likely with wrong RT
#data(threeStdsDB) #TODO spin our own internal library based on standards run on same method
#DB.treated <- treat.DB(DB)
#allSam.matches <- 
#  matchSamples2DB(allSamples.msp, DB = DB.treated, 
#                  settings = metaSetting(TSQXLS.GC, "match2DB"), 
#                  quick = FALSE)

## Match samples to samples to find unknowns
allSamples.msp.scaled <- lapply(allSamples.msp, treat.DB, isMSP = FALSE)
save(allSamples.msp.scaled, file=paste0(cdfdir, "/allSamples_msp_scaled.rda"))

btnSamp.settings <- metaSetting(TSQXLS.GC, "betweenSamples")
btnSamp.settings$min.class.fraction <- 0.1 #at least 10% of samples
btnSamp.settings$min.class.size <- 60 #at least 60 samples
btnSamp.settings$simthresh <- 0.9 #default 0.95 match
btnSamp.settings$rtdiff <- 5/60 #5 s, default 3 s

#under the hood - matches all samples to each other
library(microbenchmark)
microbenchmark(metaMS:::match.unannot.patterns(allSamples.msp[[1]], allSamples.msp[[2]], 
                                               settings=btnSamp.settings), times=5) #664 ns

# First cut to the data - retention time < 18
bind_rows(map(allSamples.msp[[1]], as_tibble), .id="ps_id") %>% 
  ggplot(aes(x=rt, y=mz, color=maxo^(1/4)))+geom_point(size=0.2)+
  scale_color_viridis_c()+theme_dark()+geom_vline(xintercept=18)

# Second cut to the data - sum of pseudospectrum < 2e5
bind_rows(map(allSamples.msp[[1]], as_tibble), .id="ps_id") %>% 
  group_by(ps_id) %>% summarize(maxo=sum(maxo)) %>% 
  ggplot(aes(x=maxo))+ geom_histogram(binwidth=1e5) + geom_vline(xintercept=2e5)

cutfilter <- map(allSamples.msp, ~map_lgl(.x, function(ps) {
  ps[,"rt"][1] <= 18 & sum(ps[,"maxo"]) >= 0}))#1e5
hist(map_dbl(cutfilter, length), main="Unfiltered pseudospectra in each sample", xlim=c(0,500))
hist(map_dbl(cutfilter, sum), main="Filtered pseudospectra in each sample", xlim=c(0,500))

#Cut and benchmark again
allSamples.msp.cut <-        map2(allSamples.msp,        cutfilter, function(x,y) x[y])
allSamples.msp.scaled.cut <- lapply(allSamples.msp.cut, treat.DB, isMSP = FALSE)

microbenchmark(metaMS:::match.unannot.patterns(allSamples.msp.cut[[1]], allSamples.msp.cut[[2]], 
                                               settings=btnSamp.settings), times=5) # 86 ns after RT

# profvis says most of the time is spent in match.unnannot.patterns, and total time scales quadratically
runS2S <- function(n) {
  matchSamples2Samples(allSamples.msp.scaled.cut[1:n], 
                       allSamples.msp.cut[1:n], 
                       annotations = NULL,
                       settings = btnSamp.settings)
  }

timings <- tibble(n=c(2, 3, 4, 5, 10)) %>% 
  mutate(time = map_dbl(n, ~ microbenchmark({runS2S(n=.x)}, times=1)$time))
ggplot(timings, aes(x=n, y=time)) + geom_point() + geom_smooth(method="lm", formula=y ~ poly(x, 2, raw=TRUE))
time.poly <- lm(time ~ poly(n, 2), data=timings)
(predtimings <- tibble(n=c(2,10,100,670)) %>% 
  mutate(time_s = predict(time.poly, newdata=.)/1e9))

#very slow (est 56 hrs, never ended) on original allSamples.msp
#one core, scales quadratically with samples
#now with cut dataset, estimate 3 hrs
nt <- 50 #number of samples for testing
allSam.matches <- 
  matchSamples2Samples(allSamples.msp.scaled.cut[1:nt], 
                       allSamples.msp.cut[1:nt], 
                       annotations = NULL,#allSam.matches$annotations, 
                       settings = btnSamp.settings)

## Make dataframe of unknowns
features.df <- getFeatureInfo(stdDB = DB, allMatches = allSam.matches, 
                              sampleList = allSamples.msp.cut[1:nt])

## Get list of pseudospectra of unknowns
PseudoSpectra <- constructExpPseudoSpectra(allMatches = allSam.matches, 
                                           standardsDB = DB)

## Get intensities relative to pseudospectrumcode 
ann.df <- getAnnotationMat(exp.msp = allSamples.msp[1:nt], pspectra = PseudoSpectra, 
                           allMatches = allSam.matches)

## Convert to peak areas/heights (probably height?)
ann.df2 <- sweep(ann.df, 1, sapply(PseudoSpectra, 
                                   function(x) max(x$pspectrum[, 2])), FUN = "*")

## Search NIST database for matches to pseudospectra of unknowns
NISTpath <- "/home/john/.wine/drive_c/Program Files (x86)/NISTMS/MSSEARCH/"
NISTpath.short <- "C:\\PROG~5P2\\NISTMS\\MSSEARCH\\"
MSPfile <- "pseudospectra.msp"
metaMS::write.msp(PseudoSpectra, file=MSPfile) 
file.copy(from=MSPfile, to=NISTpath, overwrite=TRUE)
cat(paste0(NISTpath.short, "FILSPEC.FIL"), file=paste0(NISTpath, "AUTOIMP.MSD"))
cat(paste0(NISTpath.short, MSPfile, " Overwrite"), file=paste0(NISTpath, "FILSPEC.FIL"))#deleted every time
system(paste0("wine \'", NISTpath, "nistms$.exe\' /par=2"))
while(!file.exists(paste0(NISTpath,"SRCREADY.TXT")))  Sys.sleep(3)
srch <- readLines(paste0(NISTpath, "SRCRESLT.TXT"))

parse_hit <- function(hit) {
  hit %>% str_replace("Hit [0-9]*  ", "Hit") %>% 
    str_replace(";<", "; Formula: <") %>% 
    str_replace(".$", "") %>% 
    str_replace_all("<|>", "") %>% 
    str_replace_all("\\':", "\\'|") %>%
    str_replace_all(": ", ":") %>% 
    str_replace_all("; ", "\n") %>% 
    read_delim(delim = ":", col_names=c("name","value"), col_types="cc") %>% 
    mutate(HitN = str_match(hit, "Hit ([0-9]*)")[1,2])
}

get_hits <- function(srch){
  headers <- str_which(srch, fixed("Unknown: "))
  chemnames <- map_chr(srch[headers], ~ str_match(.x, "Unknown: (.*?)   ")[1,2])
  cil <-   map_dbl(srch[headers], ~ as.integer(str_match(.x, "Compound in Library Factor = ([-]?[0-9]*)$")[1,2]))
  result <- as.list(headers) # construct with same size
  for(i in 1:length(headers)) {
    start <- headers[i] + 1
    end <- c(headers, length(srch)+1)[i+1] - 1 # add end line
    if(start < end) result[[i]] <- map_dfr(srch[start:end], parse_hit)
    else            result[[i]] <- tibble(name="Hit", value=NA, HitN="1") #*No Hits found*
    result[[i]]$Name <- chemnames[i]
    result[[i]]$CIL <- cil[i]
  }
  bind_rows(result) %>% pivot_wider() %>% type_convert() %>% 
    mutate(CAS = na_if(CAS, "0-00-0"), 
           RI = na_if(RI, 0))
}

nist.results <- get_hits(srch)

## Combine peak intensities, unknowns metadata, and NIST search results
ann.df3 <- bind_cols(as.data.frame(ann.df2), left_join(features.df, filter(nist.results, HitN==1)))

## Make table of peak intensities with samples in rows and compounds in columns
final.df <- ann.df3 %>% 
  select(Name, Hit, all_of(1:nt)) %>% #tools::file_path_sans_ext(basename(cdffiles))
  mutate(Hit = coalesce(Hit, Name)) %>% 
  group_by(Hit) %>%  summarize(across(where(is.numeric), sum)) %>% 
  column_to_rownames("Hit") %>% 
  t() %>% as.data.frame()
rownames(final.df) <- tools::file_path_sans_ext(basename(cdfs[1:nt]))


#final.df.cut <- final.df[,-c(143,144, 105)]#problem character encodings
final.df.cut <- final.df[,colSums(final.df>0)> nt/5] 
library(pheatmap)
pheatmap(t(as.matrix(final.df.cut))^(1/4), 
         annotation_col=pd[1:nt,] %>% column_to_rownames("sample_name"))
mean(final.df$`(1R)-2,6,6-Trimethylbicyclo[3.1.1]hept-2-ene`>0)
