datadir <- "Data/SupplementalData/Homology"

f <- c("ComplexPairwise_PID.RData", "ComplexPairwise_BlastPID.RData",
       "ModulePairwise_PID.RData", "ModulePairwise_BlastPID.RData",
       "ModulePairwiseFP_PID.RData", "ModulePairwiseFP_BlastPID.RData")

n <- c("Complex", "ComplexBlast",
       "Module", "ModuleBlast",
       "ModuleFP", "ModuleFPBlast")

PIDs <- list()
for(i in seq_len(6)){
  load(file.path(datadir, f[i]))
  if(i < 5){
    l <- list(Positive=PID_positive, Negative=PID_negative)
  } else {
    l <- list(Positive=PID_falsepositive)
  }
  PIDs[[n[i]]] <- l
}

save(PIDs, file=file.path(datadir, "AllPIDs.RData"))
