library(httr)
#Args <- commandArgs(trailingOnly=TRUE)
#Sys.sleep(sample(10,1))
#org <- Args[2]
outdir <- '/Volumes/OrthologData/RefSeqDb/Proteomes'
#system(paste0("touch ", outf))

retry_until_fail <- function(link, maxTries=10L, returnErr=FALSE){
  require(httr)
  WAIT_TIME_SECS <- 30L
  r <- GET(link)
  ctr <- 0
  while (r$status_code %/% 100 != 2 && r$status_code != 404 ){
    if (ctr == maxTries){
      if (returnErr)
        return(r)
      else
        return(NULL)
    }
    ctr <- ctr + 1
    Sys.sleep(WAIT_TIME_SECS + sample(3,1))
    r <- GET(link)
  }

  if (!returnErr && r$status_code==404){
    return(NULL)
  }

  return(r)
}

for(j in seq_along(orgs)){
  org <- orgs[j]
  cat(org, ':\n')
  outf <- file.path(outdir, paste0(org, "_prot.fasta"))
  if(file.exists(outf)) next
  q <- paste0("rest.kegg.jp/list/", org)
  r <- retry_until_fail(q)
  if(r$status_code == 400) stop("initial query failed")
  genes <- read.table(text=content(r), sep='\t')
  genes <- genes[genes[,2]=="CDS",1]
  Sys.sleep(1)
  if(length(genes) > 0){
    queries <- character(length(genes)/10 + 1L)
    for(i in seq_along(queries)){
      cat('\t', i, '/', length(queries), '\r')
      g <- genes[seq((i-1)*10+1, i*10)]
      g <- paste(g, collapse='+')
      q <- paste0("rest.kegg.jp/get/", g, "/aaseq")
      r <- retry_until_fail(q)
      if(r$status_code==400) stop("protein query", i, "failed")
      cat(content(r), file=outf, append=TRUE)
      if(i%%3 == 0) Sys.sleep(1)
      #Sys.sleep(sample(5,1))
    }
  }
  cat('\n')
}



