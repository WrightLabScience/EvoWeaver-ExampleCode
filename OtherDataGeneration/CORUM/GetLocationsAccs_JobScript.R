## This job gets the gene indices and the accessions for a KEGG organism
## Takes as input ResultBlast_[orgname].RData (see other JobScript)
## returns as output ResultLoc_[orgname].RData

Args <- commandArgs(trailingOnly=TRUE)

Process <- Args[1L]
Orgname <- Args[2L]
SLEEPTIMER <- 0L

library(httr)

load(paste0("ResultBlast_", Orgname, ".RData"))

# Get all genes from a chromosome:
# elink -target protein -db nuccore -id NC_002745.2 | esummary | xtract -pattern DocumentSummary -element AccessionVersion

# get assembly name for gene:
# esearch -db Assembly -query NC_003140.1 | esummary | xtract -pattern DocumentSummary -element AssemblyAccession

# query datasets for accession:
# GCF_000009645.1

safe_request_get <- function(q, ntry){
  ctr <- 0L
  r <- GET(q)
  while(r$status_code == 429 || r$status_code == 502){
    ctr <- ctr + 1L
    if(ctr == ntry) stop("Failed after ", ntry, " requests :(")
    print("Too many requests / bad gateway, pausing...")
    #Sys.sleep(sample(SLEEPTIMER,1)*10L)
    r <- GET(q)
  }
  if(r$status_code != 200){
    stop("Could not get data, error code ", r$status_code)
  }
  content(r)
}

get_correct_accession <- function(Acc){
  if(!grepl("^GC[FA]", Acc)){
    #Sys.sleep(sample(SLEEPTIMER, 1)*10)
    # First get the correct accession id for the gene
    q <- paste0("esearch -db Assembly -query ",
                Acc,
                " | esummary | xtract -pattern DocumentSummary -element AssemblyAccession")
    Acc <- system(q, intern=TRUE)
  }
  Acc
}

get_gene_locs <- function(Acc){
  require(httr)
  API_KEY <- "51f197ed70f7c06fe306ceca8250ed70ad09"
  cat("Accession is", Acc, '\n')

  ## Double check that the accession is valid
  NCBI_API <- "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"
  q <- paste(NCBI_API, Acc, "check", sep='/')
  q <- paste0(q, "?api_key=", API_KEY)
  r <- safe_request_get(q, 10L)
  if(is.null(r$valid_assemblies) || !(Acc %in% r$valid_assemblies)){
    stop("Assembly ", Acc, " was invalid")
  }

  ## Get all the gene calls
  q <- paste(NCBI_API, Acc, "annotation_report", sep='/')
  q <- paste0(q, "?api_key=", API_KEY)
  next_page <- ''
  total_count <- 0
  num_found <- 0
  req_data <- NULL
  while(!is.null(next_page)){
    if(next_page!='')
      qs <- paste0(q, "&page_token=", next_page)
    else
      qs <- q
    r <- safe_request_get(qs, 10)
    if(total_count == 0) total_count <- r$total_count
    num_found <- num_found + length(r$reports)
    cat("\t- Found", num_found, '/', total_count, "genes\n")
    next_page <- r$next_page_token
    more_data <- lapply(r$reports, \(x){
      m <- matrix('', nrow=0, ncol=5)
      ranges <- x$annotation$genomic_regions
      proteins <- x$annotation$proteins
      assemblies <- x$annotation$annotations
      protein_names <- vapply(proteins, \(y) y$accession_version, character(1L))
      if(length(protein_names) < length(ranges)){
        p <- character(length(ranges))
        p[seq_along(protein_names)] <- protein_names
      }
      assemblies <- vapply(assemblies, \(y) y$assembly_accession, character(1L))
      for(ii in seq_along(ranges)){
        p <- ranges[[ii]]$gene_range$range
        for(jj in seq_along(p))
          m <- rbind(m, c(p[[jj]]$begin, p[[jj]]$end, p[[jj]]$orientation, protein_names[ii], assemblies[ii]))
      }
      m
    })
    more_data <- do.call(rbind, more_data)
    if(is.null(req_data)){
      req_data <- as.data.frame(more_data)
    } else {
      req_data <- rbind(req_data, more_data)
    }
  }

  colnames(req_data) <- c("start", "end", "orientation", "accession", "chromosome")
  req_data <- req_data[order(as.integer(req_data$start), as.integer(req_data$end)),]
  req_data
}

ids <- RETURNDATA$GenomeIDs
if(length(ids) > 100) quit('no')
print(ids)
print(Process)
print(Orgname)
for(i in seq_along(ids)){
  id <- get_correct_accession(ids[i])
  if(length(id) == 1){
    ids[i] <- id
  } else {
    ids[i] <- ''
  }
}
ids <- unique(ids)
ids <- ids[ids != '']
# if(length(ids) > 0){
# for(i in seq_along(ids)){
#   tmp <- get_gene_locs(ids[i])
#   if(i == 1){
#     gene_locs <- tmp
#   } else {
#     gene_locs <- rbind(gene_locs, tmp)
#   }
# }
# gene_locs$chromosome_index <- match(gene_locs$chromosome, unique(gene_locs$chromosome))
# }
bbh <- RETURNDATA$BBH
orgname <- RETURNDATA$Org

unique_genes <- unique(bbh[,2])
p <- vapply(unique_genes, \(x){
  which.max((1-bbh$evalue) * (bbh$Subj_seqid==x))
}, integer(1L))
bbh <- bbh[p,]

# if(length(ids) > 0){
# new_names <- vector('list', nrow(bbh))
# for(i in seq_len(nrow(bbh))){
#   gene <- bbh[i,2]
#   pos <- which(gene_locs$accession == gene)
#   r <- character(length(pos))
#   for(j in seq_along(pos)){
#     p <- pos[j]
#     r[j] <- paste(orgname,
#                   gene_locs$chromosome_index[p], # chromosome index
#                   ifelse(gene_locs$orientation[p]=="plus", "0", "1"), # direction
#                   p, # index
#                   sep="_")
#   }
#   new_names[[i]] <- r
# }
# } else {
# new_names <- NULL
# }

save(ids, bbh, Process, Orgname, file=paste0("ResultLoc_", Orgname, ".RData"))
