library(xlsx)
curdir <- getwd()
datadir <- file.path(curdir, "Data","OtherData","RunCladeOScope")
Taxa <- read.xlsx(file.path(dirname(datadir), 'KEGGOrganismBreakdown.xlsx'), sheetIndex=1L)
Taxa <- Taxa[,c(2,4:7)]

USE_ONLY_EUKARYA <- TRUE

AllResultsFileName <- file.path(curdir, "Data", "SupplementalData", "CORUM", "ExternalAlgorithmResultsCORUM.RData")
load(AllResultFileName)
bitscores <- AllResults$bitscores

all_taxa <- colnames(bitscores)
Taxa <- Taxa[Taxa[,1] %in% all_taxa,]

if(USE_ONLY_EUKARYA){
  Taxa <- Taxa[Taxa[,2]=="Eukaryotes",]
  bitscores <- bitscores[,Taxa[,1]]
}

idMapping <- seq_len(nrow(Taxa))
names(idMapping) <- Taxa[,1]

TaxaJsonFile <- file.path(datadir,
                          ifelse(USE_ONLY_EUKARYA,
                                 "clades_euktaxid.json",
                                 "clades_taxid.json"))
cat("{\n", file=TaxaJsonFile)

for(i in seq(5,2)){
  all_taxa <- unique(Taxa[,i])
  for(tax in all_taxa){
    if(tax == "Unknown") next
    ids <- Taxa[Taxa[,i]==tax,1]
    idstring1 <- paste(ids, collapse='",\n\t\t\t"')
    idstring2 <- paste(idMapping[ids], collapse='",\n\t\t\t"')

    header <- paste0('\t"', tax, '": {\n\t\t"size": ', length(ids), ',\n')
    taxname <- paste0('\t\t"taxname": [\n\t\t\t"', idstring1, '"\n\t\t]\n\t},\n')
    idstring <- paste0('\t\t"taxid": [\n\t\t\t"', idstring2, '"\n\t\t],\n')
    cat(header, file=TaxaJsonFile, append=TRUE)
    cat(idstring, file=TaxaJsonFile, append=TRUE)
    cat(taxname, file=TaxaJsonFile, append=TRUE)
  }
}
cat("}", file=TaxaJsonFile, append=TRUE)

bitscores <- AllResults$bitscores
if(USE_ONLY_EUKARYA){
  bitscores <- bitscores[,Taxa[,1]]
}
bitscores[bitscores<20.4] <- 0
bitscores <- cbind(hsa=bitscores[,'hsa'], bitscores[,colnames(bitscores) != "hsa"])
normv <- apply(bitscores, 1, max)
normv[normv==0] <- 1
colnames(bitscores) <- as.character(idMapping[colnames(bitscores)])
bitscores <- bitscores / normv
#bitscores <- bitscores[,-1]
write.table(bitscores, sep='\t', quote=FALSE,
            file=file.path(datadir, ifelse(USE_ONLY_EUKARYA,
                                           "corum_euk_bitscores.tsv",
                                           "corum_bitscores.tsv")))

bitscores[bitscores!=0] <- log2(bitscores[bitscores!=0])
bitscores <- t(bitscores)
bitscores <- bitscores - colMeans(bitscores)
bitscoresSd <- apply(bitscores, 1L, sd)
bitscoresSd[bitscoresSd == 0] <- 1
bitscores <- bitscores / bitscoresSd
bitscores <- t(bitscores)
write.table(bitscores, sep='\t', quote=FALSE,
            file=file.path(datadir, ifelse(USE_ONLY_EUKARYA,
                                           "corum_euk_npp.tsv",
                                           "corum_npp.tsv")))

## run cladeoscope before running these lines
## note that the json file has to be manually modified to remove trailing comma
stop("Run cladeoscope now")

## these were run with the cloned repository, you'll have to run these manually and swap in the output file names
if(USE_ONLY_EUKARYA){
  outname <- "CladeOScopeEuk"
  #results <- readLines("CladeOScope/outputdata-euk.csv")
} else {
  outname <- "CladeOScopeAll"
  #results <- readLines("CladeOScope/outputdata-allscores.csv")
}
results <- strsplit(results, ',')[[1]]
results <- as.numeric(results)
cladeoscores <- matrix(0, nrow=nrow(bitscores), ncol=nrow(bitscores))
cladeoscores[lower.tri(cladeoscores)] <- results
cladeoscores <- t(cladeoscores)
cladeoscores[lower.tri(cladeoscores)] <- results
rownames(cladeoscores) <- colnames(cladeoscores) <- rownames(bitscores)

cladepairscores <- numeric(nrow(AllResults$Dataset))
for(i in seq_along(cladepairscores)){
  r <- AllResults$Dataset[i,]
  cladepairscores[i] <- cladeoscores[r[[1]], r[[2]]]
}

source(file.path(curdir, 'HelperScripts', 'PredictionCheck.R'))

res <- vcheckans(cladepairscores, AllResults$Dataset$isTP)
res$AUROC
AllResults[[outname]] <- list(Result=res, Scores=cladepairscores)
#save(AllResults, file=AllResultFileName)

plot(res$FPR, res$TPR, xaxs='i', yaxs='i', type='l', col='red')
abline(a=0,b=1)
