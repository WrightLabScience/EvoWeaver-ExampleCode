## random seeds, just using different ones for modules vs sequences vs CORUM
if(USE_COMPLEXES){
  seedval <- 871L
  fname <- "ComplexSequenceSet.RData"
  outfile <- "ComplexPairwise_BlastPID.RData"
} else {
  seedval <- 919L
  fname <- "ModuleSequenceSet.RData"
  outfile <- "ModulePairwise_BlastPID.RData"
}
if(USE_CORUM){
  seedval <- 421L
  fname <- "CorumSequenceSet.RData"
  outfile <- "CorumPairwise_BlastPID.RData"
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
tf <- tempfile()
for(i in seq_len(NUM_SAMPLES)){
  cat(i, '/', NUM_SAMPLES, '\r')
  mod_sample <- sample(s, 1)
  seqs_sample <- sample(length(aalist_nogap[[mod_sample]]), 2)
  seqset <- aalist_nogap[[mod_sample]][seqs_sample]
  names(seqset) <- c("Seq1", "Seq2")
  writeXStringSet(seqset, tf)
  x <- read.table(text=system(paste0("blastp -query ", tf, " -subject ", tf, " -outfmt 6"), intern=TRUE))
  x <- x[x[,1] != x[,2],]
  # Seq1 vs. Seq2, third column is PID, 12th is bitscore
  # max bitscore among all matches
  PID_positive[i] <- ifelse(nrow(x)==0, 0, max(x[,12]))

  if(i%%2000 == 0)plot.ecdf(log(PID_positive[seq_len(i)]+1), verticals=TRUE, pch='.', main=i)
}
cat("\n")
plot.ecdf(log(PID_positive+1), verticals=TRUE, pch='.', main="positives")
## build negative samples
for(i in seq_len(NUM_SAMPLES)){
  cat(i, '/', NUM_SAMPLES, '\r')
  mod_sample <- sample(s, 2)
  seqs_sample1 <- sample(length(aalist_nogap[[mod_sample[1]]]), 1)
  seqs_sample2 <- sample(length(aalist_nogap[[mod_sample[2]]]), 1)
  seqset <- c(aalist_nogap[[mod_sample[1]]][seqs_sample1],
              aalist_nogap[[mod_sample[2]]][seqs_sample2])
  names(seqset) <- c("Seq1", "Seq2")
  writeXStringSet(seqset, tf)
  x <- read.table(text=system(paste0("blastp -query ", tf, " -subject ", tf, " -outfmt 6"), intern=TRUE))
  x <- x[x[,1] != x[,2],]
  # colnames(x) <- c("QueryID", "SubjectID", "Perc.Ident",
  #                  "Alignment.Length", "Mismatches", "Gap.Openings",
  #                 "Q.start", "Q.end", "S.start", "S.end", "E", "Bits")
  PID_negative[i] <- ifelse(nrow(x)==0, 0, max(x[,12]))
  if(i%%2000 == 0){
    plot.ecdf(log(PID_positive+1), verticals=TRUE, pch='.', main=i, col='green')
    plot.ecdf(log(PID_negative[seq_len(i)]+1), verticals=TRUE, pch='.', col='red', add=TRUE)
  }
}
#PID_negative <- 1-PID_negative

cat("Positive pairs:\n")
v <- c(sum(PID_positive > 0.2), sum(PID_positive > 0.4))
cat('\t', v[1], ' (', v[1]/NUM_SAMPLES, '%) > 20% PID\n', sep='')
cat('\t', v[2], ' (', v[2]/NUM_SAMPLES, '%) > 40% PID\n', sep='')
cat("Negative pairs:\n")
v <- c(sum(PID_negative > 0.2), sum(PID_negative > 0.4))
cat('\t', v[1], ' (', v[1]/NUM_SAMPLES, '%) > 20% PID\n', sep='')
cat('\t', v[2], ' (', v[2]/NUM_SAMPLES, '%) > 40% PID\n', sep='')

plot.ecdf(PID_positive, verticals=TRUE, pch='.', main=i, col='green')
plot.ecdf(PID_negative, verticals=TRUE, pch='.', col='red', add=TRUE)
save(PID_positive, PID_negative, file=file.path(outdir, outfile))
