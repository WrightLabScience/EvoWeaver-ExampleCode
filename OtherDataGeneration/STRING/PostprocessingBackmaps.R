load('StringToKEGG.RData') # generated from readingSTRINGdata.R

StringSub <- StringToKEGG[StringToKEGG$KOGroup!='',]
cogs <- tapply(StringSub$KOGroup, StringSub$orthogroup, list)

maj <- numeric(length(cogs))
for(i in seq_along(maj)){
  maj[i] <- max(table(cogs[[i]]) / length(cogs[[i]]))
}

# Trying a consensus annotation
cutoff <- 1
subCogs <- cogs[maj >= cutoff]
SubCogIds <- character(length(subCogs))
for(i in seq_along(SubCogIds)){
  tab <- table(subCogs[[i]]) / length(subCogs[[i]])
  SubCogIds[i] <- names(tab)[which.max(tab)]
}
names(SubCogIds) <- names(subCogs)

# 0.5 cutoff
Cog50 <- SubCogIds
# 1.0 cutoff
Cog100 <- SubCogIds
save(cogs, consensus_prob, Cog50, Cog100, file='StringCOGs.RData')
