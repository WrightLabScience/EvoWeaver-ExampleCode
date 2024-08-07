## random seeds, just using different ones for modules vs sequences vs CORUM
if(USE_COMPLEXES){
  seedval <- 871L
  fname <- "ComplexSequenceSet.RData"
  outfile <- "ComplexPairwise_PID.RData"
} else {
  seedval <- 919L
  fname <- "ModuleSequenceSet.RData"
  outfile <- "ModulePairwise_PID.RData"
}
if(USE_CORUM){
  seedval <- 421L
  fname <- "CorumSequenceSet.RData"
  outfile <- "CorumPairwise_PID.RData"
}

load(file.path(datadir, fname), v=T)
if(USE_CORUM){
  aalist_nogap <- CORUMSequenceSet
  rm(CORUMSequenceSet)
}
# sequence sets must have at least 3 sequences
aalist_nogap <- aalist_nogap[lengths(aalist_nogap) > 2]

set.seed(seedval)
PID_positive <- PID_negative <- numeric(NUM_SAMPLES)

## build positive samples
s <- seq_along(aalist_nogap)
for(i in seq_len(NUM_SAMPLES)){
  cat(i, '/', NUM_SAMPLES, '\r')
  mod_sample <- sample(s, 1)
  seqs_sample <- sample(length(aalist_nogap[[mod_sample]]), 2)
  p <- AlignProfiles(aalist_nogap[[mod_sample]][seqs_sample[1]],
                     aalist_nogap[[mod_sample]][seqs_sample[2]],
                     processors=NULL)
  PID_positive[i] <- DistanceMatrix(p, method='shortest', includeTerminalGaps=TRUE,
                      processors=NULL, verbose=FALSE)[2]
}
PID_positive <- 1-PID_positive
cat("\n")

## build negative samples
for(i in seq_len(NUM_SAMPLES)){
  cat(i, '/', NUM_SAMPLES, '\r')
  mod_sample <- sample(s, 2)
  seqs_sample1 <- sample(length(aalist_nogap[[mod_sample[1]]]), 1)
  seqs_sample2 <- sample(length(aalist_nogap[[mod_sample[2]]]), 1)
  p <- AlignProfiles(aalist_nogap[[mod_sample[1]]][seqs_sample1],
                     aalist_nogap[[mod_sample[2]]][seqs_sample2],
                     processors=NULL)
  PID_negative[i] <- DistanceMatrix(p, method='shortest', includeTerminalGaps=TRUE,
                                    processors=NULL, verbose=FALSE)[2]
}
PID_negative <- 1-PID_negative

cat("Positive pairs:\n")
v <- c(sum(PID_positive > 0.2), sum(PID_positive > 0.4))
cat('\t', v[1], ' (', v[1]/NUM_SAMPLES, '%) > 20% PID\n', sep='')
cat('\t', v[2], ' (', v[2]/NUM_SAMPLES, '%) > 40% PID\n', sep='')
cat("Negative pairs:\n")
v <- c(sum(PID_negative > 0.2), sum(PID_negative > 0.4))
cat('\t', v[1], ' (', v[1]/NUM_SAMPLES, '%) > 20% PID\n', sep='')
cat('\t', v[2], ' (', v[2]/NUM_SAMPLES, '%) > 40% PID\n', sep='')

save(PID_positive, PID_negative, file=file.path(outdir, outfile))
