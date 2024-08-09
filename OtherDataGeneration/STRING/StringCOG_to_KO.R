source('KEGGInterface.R')

load('StringToKEGG.RData')
allcogs <- 'COG.links.detailed.v12.0.txt'

#cogtable <- read.table(allcogs)
all_cogs_string <- tapply(StringToKEGG[,5], StringToKEGG$orthogroup, list)

string_KOs <- vector('character', length(all_cogs_string))
for(i in seq_along(string_KOs)){
  cat(i, '/', length(string_KOs), '\r')
  if(grepl("^K[0-9]{5}$", string_KOs[i]))
    next
  res <- NA_character_
  j <- 1L
  geneset <- all_cogs_string[[i]]
  while(is.na(res)){
    res <- ko_for_gene(geneset[j])
    if(is.null(res) || !grepl("^K[0-9]{5}$", res)) res <- NA_character_
    j <- j + 1L
  }
  string_KOs[i] <- res
}
