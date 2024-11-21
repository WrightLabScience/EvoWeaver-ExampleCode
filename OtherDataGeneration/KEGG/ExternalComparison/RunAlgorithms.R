basedir <- getwd()
load(file.path(basedir, 'Data','SupplementalData','CORUM/NormalizedBLASTProfiles.RData'))
load(file.path(basedir,'Data','SupplementalData','CORUM/CORUMNuclearPredictions.RData'))
source(file.path(basedir, 'Data','HelperScripts','PredictionCheck.R'))

# ## pairs are in FinalDataset
# all_genes <- unique(c(FinalDataset[,1], FinalDataset[,2]))
#
# subProfiles <- HSAProfiles[all_genes]
#
# all_genomes <- lapply(subProfiles, \(x) x$Genome) |> unlist() |> unique()
#
# bitscores <- matrix(0, nrow=length(subProfiles), ncol=length(all_genomes))
# colnames(bitscores) <- all_genomes
# rownames(bitscores) <- all_genes
#
#
# for(i in seq_along(subProfiles)){
#   v <- subProfiles[[i]]
#   if(is.null(v)) next
#
#   bitscores[i,v$Genome] <- v$Bitscore
# }

## human self-hit bitscores are modified in GetHSASelfHit.R
AllResultFileName <- file.path(basedir, "Data","SupplementalData","CORUM","ExternalAlgorithmResultsCORUM.RData")
load(AllResultFileName, v=T)
bitscores <- AllResults$bitscores

###########
# SVD-phy #
###########

svdBits <- bitscores
# SVD-phy thresholds to 60
svdBits[svdBits<=60] <- 0

# each protein is normalized by its best hit
svdBits <- svdBits / apply(svdBits, 1, max)

# then SVD on the matrix
res <- svd(svdBits)
resSVD <- res$u

## retain the first 30% of columns
resSVD <- resSVD[,seq_len(ncol(resSVD)*0.3)]
rownames(resSVD) <- rownames(svdBits)

## normalize each row to unit vectors
resSVD <- resSVD / apply(resSVD, 1, \(x) sqrt(sum(x**2)))

EucDist <- \(x, y){
  sqrt(sum((x - y)**2))
}

SVDPhyScores <- numeric(nrow(FinalDataset))
for(i in seq_along(SVDPhyScores)){
  g <- FinalDataset[i,]
  SVDPhyScores[i] <- EucDist(resSVD[g[[1]],], resSVD[g[[2]],])
}

## lower distance is better
SVDPhyScores <- 1-SVDPhyScores
SVDResult <- vcheckans(SVDPhyScores, FinalDataset$isTP)
SVDResult$AUROC

#######
# PPP #
#######

## PPP does some stuff that's challenging to translate into ROCs
## reference: https://peerj.com/articles/3712/
## their score is formed by the following:
##  - threshold at bitscore of 70
##  - normalize rows as log2(p / max(p))
##  - normalize columns as (p-mu) / sd(p)
PPPbits <- bitscores
PPPbits[PPPbits <= 70] <- 0
PPPbits <- PPPbits / apply(PPPbits, 1L, max)
PPPbits[PPPbits > 0] <- log2(PPPbits[PPPbits > 0])
PPPbits <- t(PPPbits)
pppsds <- apply(PPPbits, 1L, sd)
pppsds[pppsds==0] <- 1
pppmeans <- apply(PPPbits, 1L, mean)
PPPbits <- (PPPbits - pppmeans) / pppsds
PPPbits <- t(PPPbits)

## Prediction uses the pearson correlation of the vectors
## "positives" are the top N ranked proteins
## therefore using the pearson correlation should be fine, since the ROC
## would naturally sweep values of N
PPPScores <- numeric(nrow(FinalDataset))
for(i in seq_along(SVDPhyScores)){
  g <- FinalDataset[i,]
  PPPScores[i] <- cor(PPPbits[g[[1]],], PPPbits[g[[2]],])
}
PPPResult <- vcheckans(PPPScores, FinalDataset$isTP)
PPPResult$AUROC

#######################
# Binary Co-occurence #
#######################

## Using the PPP method to generate a binarized profile
BinaryProfile <- bitscores
BinaryProfile[BinaryProfile <= 70] <- 0
BinaryProfile[BinaryProfile > 70] <- 1
BinaryProfile[] <- as.integer(BinaryProfile)
HammingScores <- MutInfScores <- JaccardScores <- numeric(nrow(FinalDataset))
for(i in seq_along(SVDPhyScores)){
  g <- FinalDataset[i,]
  v1 <- BinaryProfile[g[[1]],]
  v2 <- BinaryProfile[g[[2]],]
  HammingScores[i] <- sum(v1 != v2)
  JaccardScores[i] <- sum(v1 & v2) / sum(v1 | v2)
  crosstab <- table(v1, v2)
  crosstab <- crosstab / sum(crosstab)

  crosstab2 <- crosstab / (rep(rowSums(crosstab), 2) * rep(colSums(crosstab), each=2))
  crosstab2[crosstab2 == 0] <- 1 # if the value is zero, make it 1 so the log makes it zero
  crosstab2[is.infinite(crosstab2)] <- 1 # same with infinite values
  crosstab2 <- log2(crosstab2)
  MutInfScores[i] <- sum(crosstab * crosstab2)
}
## this is a distance, should be a similarity
HammingScores <- 1-HammingScores
HammingResult <- vcheckans(HammingScores, FinalDataset$isTP)
MutInfResult <- vcheckans(MutInfScores, FinalDataset$isTP)
JaccardResult <- vcheckans(JaccardScores, FinalDataset$isTP)
HammingResult$AUROC
MutInfResult$AUROC
JaccardResult$AUROC

###############
# CladeOScope #
###############

## CladeOScope input expects the following:
## 1. matrix of bitscores, rows are genes and columns are organism IDs
## 2. JSON file of clades:
##    - TAXNAME: { size, taxid[], taxname}
##    - taxid is the column names in the .tsv above
##    - taxname is the organism names
##    - note that all clades specified in terms of leaves (one org can be in multiple)
##


AllResults <- list(
  Hamming=list(Result=HammingResult, Scores=HammingScores),
  Jaccard=list(Result=JaccardResult, Scores=JaccardScores),
  MutInf=list(Result=MutInfResult, Scores=MutInfScores),
  PPP=list(Result=PPPResult, Scores=PPPScores),
  SVDPhy=list(Result=SVDResult, Scores=SVDPhyScores),
  Dataset=FinalDataset,
  bitscores=bitscores
)
save(AllResults, file=file.path(basedir, "Data","SupplementalData",'CORUM','ExternalAlgorithmResultsCORUM.RData'))
