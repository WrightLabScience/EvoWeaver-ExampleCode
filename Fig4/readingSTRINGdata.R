#### Reading in STRING Data ####

# These files can be downloaded from STRING, also available on Zenodo
txtpath <- 'COG.mappings.v12.0.txt'
txtpathlinks <- 'COG.links.detailed.v12.0.txt'

# This loads in the backmapping object to later map KOs to modules
load('ModulePredsAllPairs.RData')

allcogmap <- read.table(txtpath, sep='\t')
pos_kegg <- which(grepl('KEGG', allcogmap[,5L]))

keggpos <- allcogmap[pos_kegg,]
kegg_annots <- keggpos[,5]
kegg_annots <- vapply(keggpos[,5], \(x){
  x <- strsplit(x, ' ')[[1]]
  p <- which(grepl('^[a-z]{3,4}:[^;]+$', x))
  if(length(p) > 0) return(x[p[1]])
  return('')
}, character(1L))
names(kegg_annots) <- NULL

keggpos2 <- keggpos[which(kegg_annots!=''),]
keggpos2 <- cbind(keggpos2[,1:4], kegg_annots[kegg_annots!=''])
colnames(keggpos2) <- c('protein', 'start_pos', 'end_pos', 'orthogroup', 'kegg_name')

StringToKEGG <- keggpos2

#### Finding consensus annotations for STRING COGs ####

StringSub <- StringToKEGG[StringToKEGG$KOGroup!='',]
cogs <- tapply(StringSub$KOGroup, StringSub$orthogroup, list)

maj <- numeric(length(cogs))
for(i in seq_along(maj)){
  maj[i] <- max(table(cogs[[i]]) / length(cogs[[i]]))
}

# Trying a consensus annotation
cutoff <- 0.5
subCogs <- cogs[maj >= cutoff]
SubCogIds <- character(length(subCogs))
for(i in seq_along(SubCogIds)){
  tab <- table(subCogs[[i]]) / length(subCogs[[i]])
  SubCogIds[i] <- names(tab)[which.max(tab)]
}
names(SubCogIds) <- names(subCogs)

# 0.5 cutoff
Cog50 <- SubCogIds

# 1.0 cutoff made by setting the cutoff value to 1.0
# Cog100 <- SubCogIds

#### Map STRING COGs to KEGG Modules ####

# Map cogs to modules
Cogset <- Cog50
Modset <- vector('list', length(Cogset))
for(i in seq_along(Modset)){
  pos <- which(vapply(backmapping, \(x) Cogset[i] %in% x, logical(1L)))
  Modset[[i]] <- names(backmapping)[pos]
}
names(Modset) <- names(Cogset)
SubModset <- Modset[lengths(Modset) > 0]

MappedMods <- vector('list', length(SubModset))
for(i in seq_along(MappedMods)){
  ch <- character()
  for(j in seq_along(SubModset[[i]])){
    pos <- which(vapply(ModMap, \(x) SubModset[[i]][j] %in% x, logical(1L)))
    ch <- c(ch, names(ModMap)[pos])
  }
  MappedMods[[i]] <- ch
}

MappedMods <- lapply(MappedMods, unique)
names(MappedMods) <- names(SubModset)

#### Connecting STRING Modules to EvoWeaver Data ####
tab <- read.table(txtpathlinks, header=TRUE)

SubPairs <- AllPairs[!AllPairs$HasComplex & AllPairs$NoKOOverlapCategory > 0,]
r_ind <- rep(NA_integer_, nrow(tab))
for(i in seq_len(nrow(tab))){
  cat(i, '/', nrow(tab), '\r')
  cog1 <- tab$group1[i]
  cog2 <- tab$group2[i]
  if(!(cog1 %in% names(MappedMods) && cog2 %in% names(MappedMods)))
    next
  sps <- expand.grid(MappedMods[[cog1]], MappedMods[[cog2]])
  if(nrow(sps) == 0) next
  sps <- t(apply(sps, 1, sort, decreasing=FALSE))
  ind <- rep(NA, nrow(sps))
  for(j in seq_len(nrow(sps))){
    nj <- which(sps[j,1] == Pairings$Mod1 & sps[j,2] == Pairings$Mod2)
    if(length(nj) > 0){
      ind[j] <- nj[1]
    } else {
      ind[j] <- NA
    }
  }
  ind <- ind[!is.na(ind)]
  res <- NA
  if(length(ind) > 0){
    scores <- vapply(ind, \(x) max(allpredictions[x,]), numeric(1L))
    res <- ind[which.max(scores)]
  }
  r_ind[i] <- res
}

pos_good <- which(!is.na(r_ind))
SubscoresString <- tab[pos_good,]
IndexMapping <- r_ind[pos_good]

ActualCat <- SubPairs$NoKOOverlapCategory[IndexMapping]

PredictedCategory <- vapply(IndexMapping,
                            \(x) which.max(allpredictions[x,]), integer(1L))
FullSubscoresString <- cbind(SubscoresString, Pairings[IndexMapping,],
                             ActualCat, allpredictions[IndexMapping,],
                             PredictedCategory)
EWScores <- SubPairs[IndexMapping,]

# EWScores in StringPredictions50.RData has some additional metadata added to it
# Use the file rather than this script's output
# PlotStringEW.R relies on the "MinTreeBranch" column, which is simply the minimum
# of the number of leaves in each pair of trees

# save(FullSubscoresString, EWScores, file='StringPredictions50.RData')
