random_distmat <- function(labs, minv=0.01, maxv=1){
  n <- length(labs)
  v <- runif(n*(n-1)/2, minv, maxv)
  class(v) <- 'dist'
  attr(v, 'Size') <- n
  attr(v, 'Labels') <- labs
  attr(v, 'Diag') <- FALSE
  attr(v, 'Upper') <- FALSE
  v
}

benchmark_algs_onerun <- function(norgs, ew_llength=2, algsubset=NULL){
  orgs <- as.character(seq_len(norgs))
  orgs1 <- paste(orgs, 1, sample(c(0,1),norgs,r=TRUE), sample(100,norgs,r=TRUE), sep='_')
  orgs2 <- paste(orgs, 1, sample(c(0,1),norgs,r=TRUE), sample(100,norgs,r=TRUE), sep='_')
  pa_vec1 <- sample(c(0,1), norgs, TRUE)
  pa_vec2 <- sample(c(0,1), norgs, TRUE)
  if(sum(pa_vec1) == 0)
    pa_vec1[sample(seq_len(norgs), 3)] <- 1L
  if(sum(pa_vec2) == 0)
    pa_vec2[sample(seq_len(norgs), 3)] <- 1L
  n_inter <- sum(pa_vec1 * pa_vec2)
  pa_vec1 <- orgs1[pa_vec1==1]
  pa_vec2 <- orgs2[pa_vec2==1]

  npa1 <- length(pa_vec1)
  npa2 <- length(pa_vec2)

  # trees
  spec_d <- random_distmat(orgs, 0.01, 3)
  if(is.null(algsubset) || algsubset > 4){
    d1 <- random_distmat(pa_vec1, 0.01, 3)
    d2 <- random_distmat(pa_vec2, 0.01, 3)
    d1 <- TreeLine(myDistMatrix = d1, method='NJ', processors=NULL, verbose=FALSE)
    d2 <- TreeLine(myDistMatrix = d2, method='NJ', processors=NULL, verbose=FALSE)
  } else {
    d1 <- pa_vec1
    d2 <- pa_vec2
  }

  spectree <- TreeLine(myDistMatrix = spec_d, method='NJ', processors=NULL, verbose=FALSE)

  l <- list(d1, d2)
  sub_l <- rep(c(1,2), length.out=ew_llength)
  ew <- EvoWeaver(l[sub_l], MySpeciesTree=spectree, NoWarn=TRUE)
  ncomp <- ew_llength * (ew_llength-1L) / 2L

  runtimes <- numeric(16)

  # Phylogenetic profiling
  if(is.null(algsubset)) algsubset <- 16L
  seq_lengths <- c(10L, 100L, 500L)
  algs <- character(16)
  algs[1:10] <- c("PAJaccard", "GLMI", "PAOverlap", "GLDistance",
            "RPMirrorTree", "RPContextTree", "TreeDistance",
            "GeneDistance", "MoransI", "OrientationMI")
  algs[c(11,13,15)] <- paste0("GeneVector", seq_lengths)
  algs[c(12,14,16)] <- paste0("SequenceInfo", seq_lengths)
  names(runtimes) <- algs
  for(i in seq_len(min(algsubset, 10))){
    runtimes[i] <- system.time(predict(ew, Method=algs[i], Verbose=FALSE))[3]
  }

  ctr <- 11L
  for(i in c(10, 100, 1000)){
    if(ctr > algsubset) break
    seq1 <- paste(sample(AA_STANDARD, i, replace=TRUE), collapse='')
    seq2 <- paste(sample(AA_STANDARD, i, replace=TRUE), collapse='')
    d1 <- dendrapply(d1, \(x){
      if(is.leaf(x)) attr(x,'state') <- seq1
      x
    })
    d2 <- dendrapply(d1, \(x){
      if(is.leaf(x)) attr(x,'state') <- seq2
      x
    })
    ew <- EvoWeaver(list(d1, d2), NoWarn=TRUE)
    runtimes[ctr] <- system.time(predict(ew, Method="GeneVector", Verbose=FALSE))[3]
    runtimes[ctr+1] <- system.time(predict(ew, Method="SequenceInfo", Verbose=FALSE))[3]
    ctr <- ctr+2L
  }

  r <- c(norgs, npa1, npa2, n_inter, ncomp, round(runtimes,5))
  names(r)[1:5] <- c("TotalN","N1", "N2", "NIntersect", "NPairs")
  r
}

N_PER_TRIAL <- 25L
MAX_NODES <- 10000
START_NODES <- 2000L
SKIP <- 1000L
PAIR_LENS <- seq(2L, 52L, 15L)
seq_to_test <- seq(START_NODES, MAX_NODES, by=SKIP)
benchmark_timings <- matrix(nrow=N_PER_TRIAL*length(seq_to_test)*length(PAIR_LENS), ncol=21)
ctr <- 1L
for(i in seq_to_test){
  #cat(i, "nodes\n")
  for(j in seq_len(N_PER_TRIAL)){
    for(k in PAIR_LENS){
      #cat(j, '.', k, '/', N_PER_TRIAL, '    \r', sep='')
      #if(!is.na(benchmark_timings[ctr,1])){
      #  ctr <- ctr + 1L
      #  next
      #}
      v <- benchmark_algs_onerun(i,k)
      benchmark_timings[ctr,] <- v
      if(ctr == 1){
        colnames(benchmark_timings) <- names(v)
      }
      ctr <- ctr + 1L
    }
  }
  cat('\n')
}

save(benchmark_timings, file=file.path(outdir, 'EvoWeaver_RuntimeBenchmark.RData'))
