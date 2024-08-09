# this doesn't include dereplication -- I had to deduplicate some entries
load(file='StringCOGs.RData')
load('ModulesMulticlassData.RData', v=T)
load('ModulePredsAllPairs.RData')

Cogset <- Cog100

# Map cogs to modules
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
save(MappedMods, file='String_forEW_100.RData')
