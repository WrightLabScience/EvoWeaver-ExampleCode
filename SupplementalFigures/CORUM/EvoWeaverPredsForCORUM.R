## TODO: clean up script
## This file is provided for reference -- file paths may not match up correctly
basedir <- getwd()
load(file.path(basedir, 'scripts','CORUM_BlastGeneration','orthology_groups','index.RData'), v=T)
all_used <- unique(unlist(FinalDataset[,1:2]))
dendlist <- vector('list', length(all_used))
denddir <- file.path(basedir, 'scripts','CORUM_BlastGeneration','CORUM_BlastDendrograms')
for(i in seq_along(all_used)){
  cat(i, "/", length(all_used), '\r')
  ref <- which(all_used[i] == n)
  fpath <- file.path(denddir, paste0("ResultBlastDend_", ref, ".RData"))
  load(fpath)
  if(all_used[i] != RETURNDATA$gene) stop('mismatch')
  dendlist[[i]] <- RETURNDATA$tree
}
names(dendlist) <- all_used
dend_broken <- vapply(dendlist, \(x) length(x) > 1 && attr(x, 'members') > 3, logical(1L))
pos_to_remove <- names(which(!dend_broken))
FinalDataset <- FinalDataset[!(FinalDataset$Gene1 %in% pos_to_remove) & !(FinalDataset$Gene2 %in% pos_to_remove),]
table(FinalDataset$isTP)
dendlist <- dendlist[dend_broken]
vapply(dendlist, \(x) attr(x, 'members'), integer(1L))
# complexcounts <- table(FullFinalDataset$ComplexName)
# complexcounts <- complexcounts[complexcounts > 5]
# complexcounts <- complexcounts[names(complexcounts) != "None"]
# subcomp <- sample(complexcounts, 10)
# subcomp <- names(subcomp)
# FinalDataset <- FullFinalDataset[FullFinalDataset[,3] %in% subcomp,]
# FinalDataset <- rbind(FinalDataset, FullFinalDataset[sample(which(!FullFinalDataset$isTP), nrow(FinalDataset)),])
algs <- c("PAJaccard", "PAOverlap", "GLMI", "GLDistance")
res <- matrix(NA_real_, nrow=nrow(FinalDataset), ncol=length(algs))
colnames(res) <- algs
algs <- c("RPMirrorTree", "RPContextTree", "TreeDistance", "GeneVector", "SequenceInfo")
ew <- EvoWeaver(dendlist, MySpeciesTree=SpecTree)
tmp_res <- vector('list', length(algs))
for(i in seq_along(algs)){
  print(algs[i])
  tmp_res[[i]] <- predict(ew, Method=algs[i], Subset=FinalDataset[,1:2])
}
names(tmp_res) <- algs

load(file.path(basedir, "scripts","CORUM_BlastGeneration","OrthogroupsWithIndices.RData"))
ew <- EvoWeaver(all_genes, MySpeciesTree=SpecTree)
algs <- c("GeneDistance", "MoransI", "OrientationMI")
for(i in seq_len(nrow(FinalDataset))){
  cat(i, '/', nrow(FinalDataset), '\r')
  #ewmini <- EvoWeaver(all_genes[c(FinalDataset[i,1], FinalDataset[i,2])], MySpeciesTree=Spec)
  p <- predict(ew[c(FinalDataset[i,1], FinalDataset[i,2])], Method=algs, Verbose=FALSE, MySpeciesTree=SpecTree)
  tmp_resm[i,1:2] <- vapply(p, \(x) x[[2]], numeric(1L))
}

tmp_res_coloc <- vector('list', length(algs))
for(i in seq_along(algs))
  tmp_res_coloc[[i]] <- predict(ew, Method=algs[i], Subset=FinalDataset[,1:2])
names(tmp_res_coloc) <- algs

tmp_res <- c(tmp_res, tmp_res_coloc)
## TEMPORARY STUFF
tmp_resm <- matrix(NA_real_, nrow=nrow(FinalDataset), ncol=length(tmp_res))
for(j in seq_along(tmp_res)){
  test <- tmp_res[[j]]
  for(i in seq_len(nrow(FinalDataset))){
    tmp_resm[i,j] <- test[FinalDataset[i,1], FinalDataset[i,2]]
  }
}
colnames(tmp_resm) <- names(tmp_res)
#for(i in seq_along(algs)){
#  RawResults[,algs[i]] <- tmp_resm[,i]
#}
RawResults <- tmp_resm
res <- RawResults

trueanswer <- FinalDataset$isTP

res[is.na(res)] <- 0
roc_list <- lapply(seq_len(ncol(res)), \(i) vcheckans(res[,i], trueanswer))
names(roc_list) <- colnames(res)
save(RawResults, FinalDataset, roc_list, file=file.path(basedir, 'data','newFinalData','CORUM','CORUMNuclearPredictions.RData'))
