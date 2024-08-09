load_seqs_for_groups <- function(konames){
  d <- 'All_KEGG_KOs'
  aas <- vector('list', length(konames))
  for(i in seq_along(konames)){
    print(i)
    s <- file.path(d, paste0("Results_", konames[i], '.RData'))
    #if(!file.exists(s)){
    #  cat("ERROR finding", s, '\n')
    #  return("ERROR")
    #}
    if(file.exists(s)){
      load(s)
      n <- names(aa)
      names(aa) <- gsub('([a-z]{3,4}):.*', '\\1', n)
      aa <- aa[order(names(aa))]
      aas[[i]] <- aa
    } else {
      aas[[i]] <- AAStringSet()
    }
  }

  names(aas) <- konames
  aas
}

combine_aalist <- function(lst){
  r <- AAStringSet()
  for(i in seq_along(lst)){
    r <- c(r, lst[[i]])
  }
  r
}

rowSumsOfDist <- function(dv){
  n <- attr(dv, 'Size')
  l <- attr(dv, 'Labels')
  outSums <- numeric(n)
  for(i in seq_len(n)){
    lesservals <- seq_len(i-1)
    lesservals <- n*(lesservals-1) - lesservals*(lesservals-1)/2 + i - lesservals
    greatervals <- c()
    if(i+1 <= n){
      greatervals <- seq(i+1, n)
      greatervals <- n*(i-1) - i*(i-1)/2 + greatervals - i
    }
    if(all(is.na(dv[c(lesservals, greatervals)]))){
      outSums[i] <- NA
    } else {
      outSums[i] <- sum(dv[c(lesservals, greatervals)], na.rm=TRUE)
    }
  }
  names(outSums) <- l
  return(outSums)
}

trim_paralogs <- function(xss, verbose=TRUE){
  n <- unique(names(xss))
  lst <- vector('list', length(n))
  names(lst) <- n
  if(verbose) cat('Partitioning sequences...\n')
  if(verbose) pb <- txtProgressBar(style=3, max=length(xss))
  for(i in seq_along(xss)){
    ni <- names(xss)[i]
    if(is.null(lst[[ni]])){
      lst[[ni]] <- xss[i]
    } else {
      lst[[ni]] <- c(lst[[ni]], xss[i])
    }
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) cat('\nFinding best sequences...\n')

  outSeqs <- NULL
  if(verbose) pb <- txtProgressBar(style=3, max=length(lst))
  for(i in seq_along(lst)){
    if(length(lst[[i]]) <= 2){
      choice <- 1L
    } else {
      seqs <- lst[[i]]
      st <- seqtype(seqs)
      seqs <- AlignSeqs(seqs, verbose=FALSE, processors=NULL)
      seqs <- MaskAlignment(seqs)
      seqs <- as(seqs, paste0(st, "StringSet"))
      d <- DistanceMatrix(seqs, type='dist', verbose = FALSE)
      d[is.na(d)] <- 2*max(d, na.rm=TRUE)
      sums <- rowSumsOfDist(d)
      choice <- which.min(sums)
    }

    if(is.null(outSeqs)){
      outSeqs <- lst[[i]][choice]
    } else {
      outSeqs <- c(outSeqs, lst[[i]][choice])
    }
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) cat('\n')

  outSeqs
}

build_tree_for_group <- function(seqs, verbose=TRUE){
  seqs <- trim_paralogs(seqs, verbose=verbose)
  if(verbose) cat("Aligning:\n")
  ali <- AlignSeqs(seqs, processors=NULL)
  if(verbose) cat("Distance Matrix:\n")
  d <- DistanceMatrix(ali, processors=NULL)
  if(verbose) cat("Tree:\n")
  tree <- TreeLine(myXStringSet=ali, myDistMatrix=d,
                   method='NJ', processors=NULL)
  tree
}
