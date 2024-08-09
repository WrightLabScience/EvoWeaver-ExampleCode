## This script generates Data/OtherData/KEGGOrganismBreakdown.xlsx
## This is an excel doc with the following sheets:
##  - "KEGG Organism Breakdown": taxonomy for all KEGG organisms used, as well as
##      accession IDs used in the CORUM benchmark. Note that Modules/Complexes/Multiclass
##      used location data drawn directly from KEGG.
##  - "KEGG Taxa Counts": hierarchically arranged taxonomy, with counts for each group
##  - "KEGG Taxa Counts for Genomes Missing Location Data": same as above, subset
##      to genomes that were missing location data in KEGG

library(xlsx)
BaseDir <- NULL
if(is.null(BaseDir)) stop("Please set the base directory!")
outfile <- file.path(BaseDir, "Data", "OtherData", "KEGGOrganismBreakdown.xlsx")
load(file.path(BaseDir, "Data", "OtherData", "AllKEGGCDSs.RData"))
load(file.path(BaseDir, "Data", "OtherData", "KEGGTaxonomy.RData"))
load(file.path(BaseDir, "Data", "Modules", "ModulesSpeciesTree.RData"))
load(file.path(BaseDir, "Data", "OtherData", "AllAccessionNums.RData"))
not_used_corum <- read.table(file.path(BaseDir, "Data", "OtherData", "UnusedForCORUM.txt"))[,1]

## AllCDSs is the set of all KEGG calls
allSpecs <- labels(SpecTree)

## verify that all are in AllCDSs
any(match(allSpecs, names(AllCDSs), nomatch=0) == 0)

MissingCalls <- logical(length(AllCDSs))
for(i in seq_along(AllCDSs)){
  tmp <- AllCDSs[[i]]
  if(all(is.na(tmp$start) | tmp$start=="Unknown" | is.na(tmp$end) | tmp$end=="Unknown")){
    MissingCalls[i] <- TRUE
  }
}

matching_pos <- match(KEGGTaxonomy[,2], names(AllCDSs), nomatch=0)
matching_pos <- which(matching_pos != 0)
KEGGTaxonomy <- KEGGTaxonomy[matching_pos,]

all_taxa <- strsplit(KEGGTaxonomy$taxonomy, ';')
for(i in seq_along(all_taxa)){
  if(length(all_taxa[[i]]) == 3){
    all_taxa[[i]] <- c(all_taxa[[i]], 'Unknown')
  }
}
any(lengths(all_taxa) != 4)

all_taxa_df <- do.call(rbind, all_taxa)
colnames(all_taxa_df) <- paste("Taxa Level", 1:4)

to_write <- KEGGTaxonomy
colnames(to_write) <- paste("KEGG", colnames(to_write))
to_write <- cbind(to_write[,-4], all_taxa_df)

to_write[["Has Gene Location Data"]] <- !MissingCalls

acc_ids <- character(nrow(to_write))
names(acc_ids) <- to_write[,2]
corum_accs <- vapply(SubAccessions, paste, character(1L), collapse=', ')
in_common <- intersect(names(acc_ids), names(corum_accs))
acc_ids[in_common] <- corum_accs
acc_ids[not_used_corum] <- ''
to_write[["Accession(s) searched for CORUM Benchmark"]] <- acc_ids

write.xlsx2(to_write, file=outfile,
           sheetName="KEGG Organism Breakdown",
           row.names=FALSE)

# Get counts by grouping
atd_sorted <- all_taxa_df[order(all_taxa_df[,1],
                                all_taxa_df[,2],
                                all_taxa_df[,3],
                                all_taxa_df[,4]),]

## get counts by grouping
total_rows <- length(unique(KEGGTaxonomy$taxonomy))

## first taxa level
m <- matrix('', nrow=nrow(KEGGTaxonomy), ncol=length(tax_levels))
nc <- length(tax_levels)
for(i in seq_len(4)){
  tmp <- apply(atd_sorted[,seq_len(i),drop=FALSE], 1, paste, collapse=';')
  all_entries <- unique(tmp)
  tab <- table(match(tmp, all_entries))
  all_entries <- vapply(strsplit(all_entries, ';'), .subset, character(1L), i)
  names(tab) <- all_entries
  nrs <- c(1,tab)
  nrs <- cumsum(nrs)
  for(j in seq_along(nrs[-1])){
    m[nrs[j],i] <- paste0(all_entries[j], " (", tab[j], ')')
  }
}
m <- m[apply(m, 1, \(x) !all(x=='')),]

colnames(m) <- paste("Taxa Level", 1:4)
m <- as.data.frame(m)
write.xlsx2(m, file=outfile,
            sheetName="KEGG Taxa Counts",
            row.names=FALSE,
            append=TRUE)

# Get counts by grouping, but only on organisms missing entries
all_taxa_df <- all_taxa_df[MissingCalls,]
atd_sorted <- all_taxa_df[order(all_taxa_df[,1],
                                all_taxa_df[,2],
                                all_taxa_df[,3],
                                all_taxa_df[,4]),]
## get counts by grouping
total_rows <- length(unique(KEGGTaxonomy$taxonomy))

## first taxa level
m <- matrix('', nrow=nrow(KEGGTaxonomy), ncol=length(tax_levels))
nc <- length(tax_levels)
for(i in seq_len(4)){
  tmp <- apply(atd_sorted[,seq_len(i),drop=FALSE], 1, paste, collapse=';')
  all_entries <- unique(tmp)
  tab <- table(match(tmp, all_entries))
  all_entries <- vapply(strsplit(all_entries, ';'), .subset, character(1L), i)
  names(tab) <- all_entries
  nrs <- c(1,tab)
  nrs <- cumsum(nrs)
  for(j in seq_along(nrs[-1])){
    m[nrs[j],i] <- paste0(all_entries[j], " (", tab[j], ')')
  }
}
m <- m[apply(m, 1, \(x) !all(x=='')),]

colnames(m) <- paste("Taxa Level", 1:4)
m <- as.data.frame(m)
write.xlsx2(m, file=outfile,
            sheetName="KEGG Taxa Counts for Genomes Missing Location Data",
            row.names=FALSE,
            append=TRUE)
