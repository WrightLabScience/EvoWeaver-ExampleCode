source('KEGGInterface.R')
load('AllKEGGModules.RData')
load("KEGGModuleDefinitions.RData")

allmods <- unlist(AllModules, use.names = FALSE)
ModsList <- vector('list', length(allmods))
names(ModsList) <- allmods

for(mod in allmods){
  print(mod)
  r <- get_module(mod)
  r$ortho <- NULL
  ModsList[[mod]] <- r
}

test1 <- ModsList[[1]]$defn
strsplit(test1, ' ')

#####
# 4 cases to handle
#
# 1. base case: K00001 K00002
#    This can just be strsplit on space
#
# 2. single depth parentheses (K00001,K00002) K00003
#    Content of parentheses should be collapsed
#
# 3. double depth parentheses (K00001, (K00002, K00003)) K00004
#    Here I should probably condense then un-collapse

fix_module_defn <- function(moddef){
  moddef <- strsplit(moddef, ' ')[[1]]
  i <- 1L
  l <- length(moddef)
  while(i < l){
    cat(i, l, '\n')
    num_open <- nchar(gsub('[^(]', '', moddef[i]))
    num_close <- nchar(gsub('[^)]', '', moddef[i]))
    if(num_open != num_close){
      moddef[i] <- paste(moddef[i], moddef[i+1], sep=' ')
      moddef <- moddef[-(i+1)]
    } else {
      i <- i + 1
    }
    l <- length(moddef)
  }

  .parse_piece <- function(piece){
    lst <- vector('list', length(piece))
    for(i in seq_along(piece)){
      mdi <- piece[i]
      num_open <- nchar(gsub('[^(]', '', piece[i]))
      if(num_open <= 1){
        mdi <- gsub('[(,)]', ' ', mdi)
        mdi <- strsplit(mdi, ' ')[[1]]
        mdi <- mdi[nchar(mdi) > 0]
      }
      # else {
        #mdi <- substr(mdi, 2L, nchar(mdi)-1)
        #mdi <- strsplit(mdi, ' ')[[1]]
        #mdi <- .parse_piece(mdi)
      # }
      lst[[i]] <- mdi
    }

    return(lst)
  }

  .parse_piece(moddef)
}


SplitMods <- ModsList
for(i in seq_along(SplitMods)){
  SplitMods[[i]] <- fix_module_defn(SplitMods[[i]]$defn)
}

# Remove modules comprised of modules or ambiguous codes
SubsetMods <- SplitMods[sapply(SplitMods, \(x){
  return(!any(grepl('[(M]', unlist(x))))
})]


# split up components of complexes (maybe we should concatenate?)
for(i in seq_along(SubsetMods)){
  mv <- SubsetMods[[i]]
  for(j in seq_along(mv)){
    if(any(grepl('[+-]', mv[[j]]))){
      mv[[j]] <- unlist(strsplit(mv[[j]], '[+-]'))
    }
  }
  SubsetMods[[i]] <- mv
}

# remove "" entries
for(i in seq_along(SubsetMods)){
  mv <- SubsetMods[[i]]
  for(j in seq_along(mv)){
    mv[[j]] <- (mv[[j]])[nchar(mv[[j]]) > 0]
  }
  SubsetMods[[i]] <- mv
}

# remove character(0) entries and multi-stage modules
curlen <- length(SubsetMods)
to_remove <- logical(curlen)
for(i in seq_len(curlen)){
  x <- SubsetMods[[i]]
  if(any(tmppos <- grepl('\\n', x))){
    pos <- which(tmppos)
    outlst <- vector('list', length(pos)+1L)
    pos <- c(0L, pos, length(x))
    for(j in seq_along(outlst)){
      tmp <- x[(pos[j]+1):(pos[j+1L])]
      tmp <- tmp[sapply(tmp, \(y) length(y) > 0)]
      tmp <- lapply(tmp, \(y) gsub('\\n', '', y))
      outlst[[j]] <- tmp
    }

    n <- names(SubsetMods)[i]
    newn <- paste(n, seq_along(outlst), sep='_')
    n <- c(names(SubsetMods), newn)
    SubsetMods <- c(SubsetMods, outlst)
    names(SubsetMods) <- n
    to_remove[i] <- TRUE
  }
}

SubsetMods <- SubsetMods[!to_remove]
# remove "" entries again
for(i in seq_along(SubsetMods)){
  mv <- SubsetMods[[i]]
  for(j in seq_along(mv)){
    mv[[j]] <- (mv[[j]])[nchar(mv[[j]]) > 0]
  }
  SubsetMods[[i]] <- mv[sapply(mv, \(y) length(y) > 0)]
}
SubsetMods <- SubsetMods[lengths(SubsetMods) > 1]
save(SubsetMods, file='KEGGTrimmedModuleBlocks.RData')

# Fill in modules with sequence sets
# These functions come from ModuleComboExperiment.R
source("ModuleComboExperiment.R")
outdir <- 'moduleblocks'
n <- names(SubsetMods)
SubsetModSeqs <- SubsetMods
for(i in seq_along(SubsetModSeqs)){
  cat(i, '/', length(SubsetModSeqs), '\n')
  NAME <- n[i]
  SEQS <- SubsetModSeqs[[i]]
  ofile <- file.path(outdir, paste0('Module_', NAME, '.RData'))
  if(file.exists(ofile)) next
  for(j in seq_along(SEQS)){
    tmp <- (SEQS[[j]])[grepl('K[0-9]*', SEQS[[j]])]
    tmp <- load_seqs_for_groups(tmp)
    if(!(is.character(tmp) && tmp == 'error'))
      tmp <- combine_aalist(tmp)
    #tmp <- trim_paralogs(tmp)
    SEQS[[j]] <- tmp
  }
  save(SEQS,NAME, file=ofile)
}

#### Finding complexed proteins
load('KEGGModuleDefinitions.RData')
subsetval <- vapply(ModsList, \(x) grepl('[+]', x$defn), logical(1))
SubModsList <- ModsList[subsetval]
SubModsList <- lapply(SubModsList, \(x){
  d <- x$defn
  d <- strsplit(d, ' ')[[1]]
  d[grepl('[+]', d) & !grepl('[()]', d)]
})

AllComplexes <- unlist(SubModsList, use.names = TRUE)
mistakes <- grepl('[^()K0-9+-]', AllComplexes)
# Going to have to hand curate these
mistakeEntries <- AllComplexes[mistakes]
AllComplexes <- AllComplexes[!mistakes]

# remove \n
mistakeEntries <- gsub('\\n', '', mistakeEntries)
# remove commas
mistakeEntries <- unlist(strsplit(mistakeEntries, ','))

AllComplexes <- c(AllComplexes, mistakeEntries)
KEGGModuleComplexes <- AllComplexes
save(KEGGModuleComplexes, file='KEGGModuleComplexes.RData')

KOs <- Mods <- Label <- character(0)
idxLab <- c('Optional', 'Required')
for(i in seq_along(KEGGModuleComplexes)){
  modname <- names(KEGGModuleComplexes)[i]
  entry <- KEGGModuleComplexes[i]

  # - indicates optional, + indicates required
  # order does not matter

  entry <- gsub('([+-])', ' \\1', entry)
  entry <- strsplit(entry, ' ')[[1]]
  # first entry is required if unmarked (it's always unmarked)
  if(!grepl('[+-]', entry[1]))
    entry[1] <- paste0('+', entry[1])

  Label <- c(Label, idxLab[grepl('+', entry, fixed=TRUE)+1L])
  KOs <- c(KOs, substr(entry, 2, 7))
  Mods <- c(Mods, rep(modname, length(entry)))
}

ComplexMapping <- data.frame(KO=KOs, Module=Mods, Label=Label)

save(KEGGModuleComplexes,
     ComplexMapping,
     UniqueKOs,
     file='KEGGModuleComplexes.RData')

Description <- Symbol <- character(nrow(ComplexMapping))
for(i in seq_len(nrow(ComplexMapping))){
  ko <- ComplexMapping[i,1L]
  print(ko)
  r <- get_ko_description(ko, v=FALSE)
  Symbol[i] <- r[1]
  Description[i] <- r[2]
}

ComplexMapping$Symbol <- Symbol
ComplexMapping$Description <- Description
