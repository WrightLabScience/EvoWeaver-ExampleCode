load(file.path(datadir, "SupplementalData", "CORUM", "TrainTestSet.RData"))
load(file.path(datadir, "Modules", "ModulesSpeciesTree.RData"))

FullFinalDataset <- FinalDataset
nuctransport <- read.delim(file.path(datadir, "SupplementalData", "CORUM", 'nuclear_transport_complexes.txt'))
cmplx_nuctransport <- nuctransport$ComplexName
gene_names <-  unique(unlist(strsplit(nuctransport$subunits.Gene.name., ';')))
FinalDataset <- FullFinalDataset[FullFinalDataset$ComplexName %in% cmplx_nuctransport,]
tmp <- FinalDataset[,1:2]
FinalDataset[,1] <- pmin(tmp[,1], tmp[,2])
FinalDataset[,2] <- pmax(tmp[,1], tmp[,2])
FinalDataset <- FinalDataset[!duplicated(FinalDataset[,1:2]),]
AllUsed <- unique(unlist(FinalDataset[,1:2]))
gene_names <- lapply(gene_names, \(x) x[x%in%AllUsed])
set.seed(123L)
NegPairs <- FinalDataset
for(i in seq_len(nrow(NegPairs))){
  cat(i, '/', nrow(NegPairs), '\r')
  FOUND <- FALSE
  while(!FOUND){
    p <- sample(3, 2)
    v <- c(sample(gene_names[[p[1]]], 1), sample(gene_names[[p[2]]], 1))
    v <- sort(v)
    curseen <- rbind(FinalDataset[,1:2], NegPairs[seq_len(i-1),1:2])

    if(any((v[1]==curseen[,1] & v[2] == curseen[,2]))){
      next
    }
    FOUND <- TRUE
    NegPairs[i,1:2] <- v
  }
}
NegPairs[,3] <- "None"
NegPairs[,4] <- FALSE
FinalDataset <- rbind(FinalDataset, NegPairs)
sort(table(unlist(FinalDataset[,1:2])))

genes <- unique(unlist(strsplit(nuctransport$subunits.Gene.name., ';')))
FinalDataset <- FullFinalDataset[FullFinalDataset$Gene1 %in% genes | FullFinalDataset$Gene2%in% genes,]
#Balance the dataset
set.seed(389L)
to_samp <- min(table(FinalDataset$isTP))
FinalDataset <- FinalDataset[c(sample(which(FinalDataset$isTP), to_samp),
                               sample(which(!FinalDataset$isTP), to_samp)),]
save(FinalDataset, file=file.path(outdir, "CORUMBenchmarkPairs.RData"))
