#indir <- '/Volumes/OrthologData/CORUM/blastres/'
fs <- list.files(indir, full.names=TRUE)
HSAProfiles <- list()
for(i in seq_along(fs)){
  cat(i, '/', length(fs), '\r')
  load(fs[i])
  org <- RETURNDATA$Org

  # human gene id, bitscore
  # don't need organism gene id
  if(nrow(RETURNDATA$BBH) == 0) next
  res <- RETURNDATA$BBH[,c(1,12)]
  res2 <- tapply(res[,2], res[,1], max)
  n <- names(res2)
  r <- data.frame(Genome=org, Bitscore=0.0)
  for(j in seq_along(res2)){
    g <- n[j]
    r$Bitscore[1] <- res2[j]
    if(is.null(HSAProfiles[[g]])){
      HSAProfiles[[g]] <- r
    } else {
      HSAProfiles[[g]] <- rbind(HSAProfiles[[g]], r)
    }
  }
}

save(HSAProfiles, file='Data/SupplementalData/CORUM/NormalizedBLASTProfiles.RData')
