## This script was used to generate AllPairs data for Modules on distributed compute

## We will need the following:
##	- Pairings.RData: to get the job id
##	- EvoWeaver object
##	- inds object
##  - original species tree
##	- KEGG species tree
USE_PROK <- FALSE

load('Pairings.RData') # 0MB, Pairings
load('EukaryoteEWData.RData')

## datafiles from Data/Modules/
load('ModuleStatistics.RData')
load("ModulePredsAllPairs.RData")

## match the predictions used against the samples
pos_used <- match(rownames(RawScores), rownames(AllPairs))
Subpairings <- Pairings[pos_used,]

## first create subpw objects to make processing faster
modules_used <- unique(unlist(Subpairings))

#save(pw_prok, coloc_proks, SpecTreeProk, treeProk, file="ProkaryoteEWData.RData")
#save(pw_euk, coloc_euks, SpecTreeEuk, treeEuk, file="ProkaryoteEWData.RData")

algs_nocoloc <- c("PAJaccard", "PAOverlap", "GLDistance", "GLMI",
                  "TreeDistance", "RPMirrorTree", "RPContextTree",
                  "SequenceInfo", "GeneVector")

algs_coloc <- c("GeneDistance", "MoransI", "OrientationMI")

n_total_algs <- length(algs_coloc) + length(algs_nocoloc)

if(USE_PROK){
  p_present <- vapply(pw_prok, \(x) length(labels(x)), integer(1L))
  p_present <- p_present > 2
  subpw <- pw_prok[p_present]
  coloc_modules <- coloc_proks[p_present]
  present_mods <- names(p_present[p_present])
  SpecTree <- SpecTreeProk
  tree <- treeProk
  outfile <- "ModulesProkOnlyResult.RData"
} else {
  p_present <- vapply(pw_euk, \(x) length(labels(x)), integer(1L))
  p_present <- p_present > 2
  subpw <- pw_euk[p_present]
  coloc_modules <- coloc_euks[p_present]
  present_mods <- names(p_present[p_present])
  SpecTree <- SpecTreeEuk
  tree <- treeEuk
  outfile <- "ModulesEukOnlyResult.RData"
}
Subpairings <- AllPairs[rownames(RawScores), 1:2]
pairings_used <- (Subpairings$Mod1 %in% present_mods) & (Subpairings$Mod2 %in% present_mods)
Subpairings2 <- Subpairings[pairings_used,]

pw_nocoloc_orig <- EvoWeaver(subpw, MySpeciesTree=SpecTree, NoWarn=TRUE)
pw_coloc_orig <- EvoWeaver(coloc_modules, MySpeciesTree=SpecTree, NoWarn=TRUE)
pw_nocoloc_kegg <- EvoWeaver(subpw, MySpeciesTree=tree, NoWarn=TRUE)
pw_coloc_kegg <- EvoWeaver(coloc_modules, MySpeciesTree=tree, NoWarn=TRUE)

print("Original Tree, no coloc")
set.seed(635L)
p1 <- predict(pw_nocoloc_orig, Method=algs_nocoloc, Verbose=TRUE, Subset=Subpairings2)
print("Colocalization")
p2 <- predict(pw_coloc_orig, Method=algs_coloc, Verbose=TRUE, Subset=Subpairings2)
print("KEGG tree, no coloc")
set.seed(635L)
p3 <- predict(pw_nocoloc_kegg, Method=algs_nocoloc, Verbose=TRUE, Subset=Subpairings2)
print("KEGG tree, coloc")
p4 <- predict(pw_coloc_kegg, Method=algs_coloc, Verbose=TRUE, Subset=Subpairings2)

OrigEukData <- merge(p1,p2)
KEGGEukData <- merge(p3,p4)
colnames(OrigEukData)[7] <- colnames(KEGGEukData)[7] <- "TreeDistance"
ActualTP <- AllPairs[,c(1:2,21)]
res <- merge(Subpairings2, ActualTP)
OrigEukData$isTP <- res$NoKOOverlapCategory == 1
KEGGEukData$isTP <- res$NoKOOverlapCategory == 1

source("Data/HelperScripts/PredictionCheck.R")
Results <- list(Original=list(), KEGG=list())
for(i in seq_len(ncol(OrigEukData)-3L)){
  Results$Original[[colnames(OrigEukData)[i+2L]]] <- vcheckans(OrigEukData[,i+2L], OrigEukData$isTP)
  Results$KEGG[[colnames(KEGGEukData)[i+2L]]] <- vcheckans(KEGGEukData[,i+2L], KEGGEukData$isTP)
}

vapply(Results$Original, \(x) x$AUROC, numeric(1L))
vapply(Results$KEGG, \(x) x$AUROC, numeric(1L))
Results$OrigData <- list(Orig=OrigEukData, KEGG=KEGGEukData)
Results$Pairings <- Subpairings2
save(Results, file='Data/Modules/EukOnlyModuleStatistics.RData')

## try running an ensemble model on CORUM based on this
load('Data/SupplementalData/CORUM/CORUMNuclearPredictions.RData', v=T)
library(randomForest)
library(neuralnet)
set.seed(823L)
logit_model <- glm(isTP ~ ., data=OrigEukData[,-(1:2)], family='binomial')
res <- predict(logit_model, as.data.frame(RawResults))
vcheckans(res, FinalDataset$isTP)$AUROC

logit2 <- glm(isTP ~ ., data=data.frame(RawResults, isTP=FinalDataset$isTP), family='binomial')
res <- predict(logit2, OrigEukData)
vcheckans(res, OrigEukData$isTP)$AUROC

rf <- randomForest(isTP ~ ., data=data.frame(OrigEukData[,-c(1:2,15)], isTP=as.factor(OrigEukData$isTP)), maxnodes=10)
res <- predict(rf, as.data.frame(RawResults), type='prob')[,"TRUE"]
vcheckans(res, FinalDataset$isTP)$AUROC

rf2 <- randomForest(isTP ~ ., data=data.frame(RawResults, isTP=as.factor(FinalDataset$isTP)), maxnodes=10)
res <- predict(rf2, OrigEukData, type='prob')[,"TRUE"]
vcheckans(res, OrigEukData$isTP)$AUROC

NotEuk <- RawScores[!pairings_used,]
res <- predict(logit_model, data=NotEuk)
vcheckans(res, NotEuk$isTP)$AUROC

tmp <- vcheckans(rowSums(OrigEukData[,-c(1:2,15)]), OrigEukData$isTP)
plot(tmp$FPR, tmp$TPR, type='l', xaxs='i', yaxs='i')
abline(a=0,b=1,lty=3)
tmp <- vcheckans(rowSums(RawScores), RawScores$isTP)
plot(tmp$FPR, tmp$TPR, type='l', xaxs='i', yaxs='i')
abline(a=0,b=1,lty=3)
tmp <- vcheckans(rowSums(RawResults), FinalDataset$isTP)
plot(tmp$FPR, tmp$TPR, type='l', xaxs='i', yaxs='i')
abline(a=0,b=1,lty=3)
tmp <- vcheckans(rowSums(EWScores[,1:12]), EWScores$Category <= 3)
plot(tmp$FPR, tmp$TPR, type='l', xaxs='i', yaxs='i')
plot(tmp$FPR, tmp$TPR, type='l', xaxs='i', yaxs='i', xlim=c(0,0.02))
abline(a=0,b=1,lty=3)

tmp <- OrigEukData[,-c(1:2,15)]
tmp <- RawScores
tmp <- cbind(rowMeans(tmp[,1:4]), rowMeans(tmp[,5:7]), rowMeans(tmp[,9:10]), rowMeans(tmp[,11:12]))
vcheckans(rowSums(tmp), RawScores$isTP)$AUROC
