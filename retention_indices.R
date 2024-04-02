#Get Kovats retention indices of compounds in quant integrations

library(tidyverse)
quanttable <- read_tsv("data/volatiles/quant_RT_ratios - quanttable57.tsv") %>% 
  mutate(name = str_replace_all(name, ".alfa.","alpha")) #alfa-copaene

library(PubChemR)
cids <- get_cids(quanttable$name)
cids.combined <- cids %>% group_by(Identifier) %>% summarize(CID = paste(CID, collapse=", "))
cids.first <- cids %>% mutate(CID=as.integer(CID)) %>% arrange(Identifier, CID) %>% group_by(Identifier) %>% slice(1)

library(jsonlite)
get_median_RI <- function(cid) {
  pugview.url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/", cid, "/JSON/?heading=Kovats+Retention+Index")
  urlerror <- function(cond) { print(paste("Can't open connection for CID", cid)); NA }
  ri.json <- tryCatch(fromJSON(pugview.url), warning = urlerror, error=urlerror)
  if(is.na(ri.json)) return(NA)
  ri.df <- ri.json[["Record"]][["Section"]][["Section"]][[1]][["Section"]][[1]][["Information"]][[1]] 
  snp.row <- which(ri.df$Name == "Standard non-polar")
  if(length(snp.row) == 0) {
    snp.row <- which(ri.df$Name == "Semi-standard non-polar")
    if(length(snp.row) == 0) {
      print(paste("Can't find non-polar RI for CID", cid))
      return(NA)
    }
    print(paste("Used semi-standard non-polar RI for CID", cid))
  }
  ri.all <- unlist(ri.df$Value$Number[snp.row])
  return(median(ri.all))
}

quanttable.RI <- quanttable %>% left_join(cids.first, by=c("name"="Identifier")) %>% 
  mutate(RI_pubchem = map_dbl(CID, get_median_RI))

quanttable.RI %>% select(CID, RI_pubchem)

#quanttable.RI %>% mutate(maxfield = shortname %in% ipogood) %>% select(maxfield)

#from 2018 maxfield quant data with 51 volatiles
quant.RT.2018 <- read_tsv("data/volatiles/quant_RT.tsv") %>% filter(year==2018) %>% mutate(RI_2018 = get_kovats(RT.median, "2018-09-03")) %>% 
  rename(name = Name, RT_2018 = RT.median)
quanttable %>% left_join(quant.RT.2018) %>% select(name, RI_2018)
