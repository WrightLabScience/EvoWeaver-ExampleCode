## This script is to check the sequence similarity between
## high scoring class 1 misclasses in the multiclass benchmark

seedval <- 1284L
NUM_SAMPLES <- 10000L
fname <- "ModuleSequenceSet.RData"
if(USE_BLAST){
  outname <- "ModulePairwiseFP_BlastPID.RData"
} else {
  outname <- "ModulePairwiseFP_PID.RData"
}
datadir <- "Data"

load(file.path(datadir, "SupplementalData", "Homology", fname))
load(file.path(datadir, "Modules", "ModulePredsAllPairs.RData"))

# sequence sets must have at least 3 sequences
aalist_nogap <- aalist_nogap[lengths(aalist_nogap) > 2]

PID_falsepositive <- matrix(0, nrow=3, ncol=NUM_SAMPLES)

SubPairs <- which(AllPairs$Mod1 %in% names(aalist_nogap) & AllPairs$Mod2 %in% names(aalist_nogap))
SubPairs <- AllPairs[SubPairs,]
pos_fp3 <- which(SubPairs$EnsemblePredictedClass == 1 & SubPairs$NoKOOverlapCategory == 3)
pos_fp4 <- which(SubPairs$EnsemblePredictedClass == 1 & SubPairs$NoKOOverlapCategory == 4)
pos_fp5 <- which(SubPairs$EnsemblePredictedClass == 1 & SubPairs$NoKOOverlapCategory == 5)
positions <- list(pos_fp3, pos_fp4, pos_fp5)

set.seed(seedval)
tf <- tempfile()
for(i in seq_len(NUM_SAMPLES)){
  for(j in seq_len(3L)){
    cat((i-1)*3+j, '/', NUM_SAMPLES*3, '\r')
    p <- sample(positions[[j]], 1)
    m1 <- SubPairs$Mod1[p]
    m2 <- SubPairs$Mod2[p]
    s1 <- sample(aalist_nogap[[m1]], 1)
    s2 <- sample(aalist_nogap[[m2]], 1)
    if(USE_BLAST){
      seqset <- c(s1, s2)
      names(seqset) <- c("Seq1", "Seq2")
      writeXStringSet(seqset, tf)
      x <- read.table(text=system(paste0("blastp -query ", tf, " -subject ", tf, " -outfmt 6"), intern=TRUE))
      x <- x[x[,1] != x[,2],]
      val <- ifelse(nrow(x)==0, 0, max(x[,12]))
    } else {
      p <- AlignProfiles(s1, s2, processors=NULL)
      val <- DistanceMatrix(p, method='shortest', includeTerminalGaps=TRUE,
                            processors=NULL, verbose=FALSE)[2]
    }
    PID_falsepositive[j,i] <- val
  }
}

## output of DistanceMatrix line is a distance, not a similarity
if(!USE_BLAST) PID_falsepositive <- 1-PID_falsepositive
rownames(PID_falsepositive) <- c("Actual3", "Actual4", "Actual5")
PID_falsepositive <- t(PID_falsepositive)
save(PID_falsepositive, file=file.path(datadir, "SupplementalData", "Homology", outname))
