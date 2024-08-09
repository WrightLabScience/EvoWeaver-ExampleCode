library(DECIPHER)

Args <- commandArgs(trailingOnly = TRUE)
Process <- as.integer(Args[1L])+1L

load('StringCogsInList.RData')
source(KEGGInterface.R)

i <- Process
geneset <- all_cogs_string[[i]]
res <- NA_character_
for(gene in geneset){
  res <- ko_for_gene(gene)
  if(is.null(res) || !grepl("^K[0-9]{5}$", res))
    res <- NA_character_
  else break
}
RETURNDATA <- list(Process=Process, KO=res)
save(RETURNDATA, file=paste0('LookupResult_', Process, '.RData'))
