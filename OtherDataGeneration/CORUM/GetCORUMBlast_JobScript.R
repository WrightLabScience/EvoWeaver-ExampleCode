## Orgname should be the KEGG organism code
## Running this generates two files:
##  1. [orgname].fa: proteome
##  2. ResultBlast_[orgname].RData: saves BLAST hits for all human genes against this organism
##
## All the `ResultBlast` files are collated in `CORUM_BLAST_Results.RData`
## All the proteomes are available in CORUM_proteomes.zip
Args <- commandArgs(trailingOnly=TRUE)
Orgname <- Args[1L]
load("AllAccessionNums.RData")
SLEEPTIMER <- 10L

suppressPackageStartupMessages(library(SynExtend))
library(httr)

Accs <- SubAccessions[[Orgname]]
outfile_seqs <- paste0(Orgname, '.fa')

get_proteome_for_accession <- function(acc, outname=paste0(acc, '.fa')){
  # First make sure the record isn't suppressed or moved
  cmd1 <- paste0("elink -target protein -db nuccore -id ", acc)
  cmd1 <- paste0(cmd1, " | esummary")
  cmd1 <- paste0(cmd1, " | xtract -pattern DocumentSummary -element Gi")
  cat("Getting GIs...\n")
  failedlookup <- FALSE
  Sys.sleep(sample(SLEEPTIMER, 1))
  all_geneids <- system(cmd1, intern=TRUE)
  sanitycheck <- strsplit(all_geneids, '\n')
  genomeid <- acc
  if(length(sanitycheck) < 1 || length(sanitycheck[[1]]) < 1){
    # Record likely suppressed or updated
    # first find the organism name
    failedlookup <- TRUE
    orgname <- system(paste0("esearch -db nuccore -query ", acc,
                             " | esummary",
                             " | xtract -pattern DocumentSummary -element Organism"), intern=TRUE)
    if(nchar(orgname) == 0) stop("Couldn't find organism")
    cat("Missing assembly, searching for '", orgname, "'\n", sep='')
    Sys.sleep(sample(SLEEPTIMER, 1))
    genomeid <- system(paste0("esearch -db nuccore -query '(\"",
                            orgname, '"[Organism]) AND "complete genome"[All Fields]\'',
                            " | esummary | xtract -pattern DocumentSummary -element Id"), intern=TRUE)
    cat("Found ids:", paste(genomeid, collapse=', '), '\n')
    cat("Getting GIs again...\n")
    all_geneids <- character(0L)
    for(id in genomeid){
      cat("\t-", id, '\n')
      Sys.sleep(sample(SLEEPTIMER, 1))
      tmp <- system(paste0("elink -target protein -db nuccore -id ", genomeid,
                                   " | esummary | xtract -pattern DocumentSummary -element Gi"), intern=TRUE)
      all_geneids <- c(all_geneids, tmp)
    }
  }
  cat("Found", length(all_geneids), "GIs.\n")
  tf <- tempfile()
  writeLines(all_geneids, tf)
  cat("Fetching fasta...\n")
  cmd2 <- paste0("efetch -db protein -input ", tf, ' -format fasta >> ', outname)
  Sys.sleep(sample(SLEEPTIMER, 1))
  system(cmd2)
  cat("Done!\n")
  rm(tf)
  return(list(outname=outname, IDs=genomeid, failedlookup=failedlookup))
}

blast_against_proteome <- function(query, subj){
  makecmd <- paste0("makeblastdb -dbtype prot -in ", subj, ' -out blastdb')
  system(makecmd)
  querycmd <- paste0("blastp -query ", query, ' -db blastdb -outfmt 6 -evalue 0.01 -out blastresults.tsv')
  system(querycmd)
  res <- read.table("blastresults.tsv", sep='\t', header=FALSE)

  # query is CORUM, subject is genome/proteome
  colnames(res) <- c('CORUM_seqid', 'Subj_seqid', 'local pident', 'length', 'mismatch',
                     'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
  u1 <- unique(res[,1])
  tmp <- vector('list', length(u1))
  for(i in seq_along(u1)){
    cur <- res[res[,1] == u1[i],]
    tmp[[i]] <- cur[which.min(cur[,11L]),]
  }
  res <- do.call(rbind, tmp)
  system('rm blastdb*')
  system('rm blastresults.tsv')
  res
}

all_ids <- character(0L)
for(acc in Accs){
  r <- get_proteome_for_accession(acc, outfile_seqs)
  all_ids <- c(all_ids, r$IDs)
  if(r$failedlookup) break
}

if(length(Accs) == 0L){
  system(paste0("touch ", outfile_seqs))
  system(paste0("touch NOACC_", Orgname))
  BlastHits <- data.frame()
} else {
  BlastHits <- blast_against_proteome("AllHumanGenes.fa", outfile_seqs)
}

RETURNDATA <- list(BBH=BlastHits, GenomeIDs=all_ids, Org=Orgname)
save(RETURNDATA, file=paste0("ResultBlast_", Orgname, '.RData'))
