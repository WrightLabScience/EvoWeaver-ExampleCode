## Assorted functions to parse data from KEGG

library(DECIPHER)
library(httr)

BASE_STRING <- "http://rest.kegg.jp"

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
    Sys.sleep(WAIT_TIME_SECS)
    r <- GET(link)
  }

  if (!returnErr && r$status_code==404){
    return(NULL)
  }

  return(r)
}

query_db <- function(db_name){
  r <- GET(paste(BASE_STRING, 'list', db_name, sep='/'))
  stop_for_status(r)

  tex <- strsplit(content(r, 'text'), '\n')[[1]]
  identifiers <- gsub('(md:[A-Z0-9]*)\t.*', '\\1', tex)

  return(identifiers)
}

find_from_id <- function(db, id){
  req <- paste(BASE_STRING, 'find', db, id, sep='/')
  r <- retry_until_fail(req, returnErr=TRUE)
  if (stopWithError) stop_for_status(r)

  return(r)
}

get_from_id <- function(id, stopWithError=TRUE){
  req <- paste(BASE_STRING, 'get', id, sep='/')
  r <- retry_until_fail(req, returnErr=TRUE)
  if (stopWithError) stop_for_status(r)

  return(r)
}

get_module <- function(id){
  r <- get_from_id(id)
  r <- content(r, 'text')

  name <- trimws(gsub('.*NAME(.*)DEFINITION.*', '\\1', r))

  defn <- trimws(gsub('.*DEFINITION(.*)ORTHOLOGY.*', '\\1', r))
  ortho <- trimws(gsub('.*ORTHOLOGY(.*)CLASS.*', '\\1', r))
  ortho <- strsplit(ortho, '\n')[[1]]
  ortho <- gsub('[ ]*([K0-9,]*).*', '\\1', ortho)
  ortho <- strsplit(ortho, ',')

  return(list(name=name, defn=defn, ortho=ortho))
}

get_module_table <- function(id, pathway=FALSE){
  require(httr)
  require(rvest)
  if (!pathway){
    id <- gsub('md:(.*)', '\\1', id)
    query_url <- paste0('https://www.kegg.jp/kegg-bin/view_ortholog_table?md=',id)
  } else {
    id <- gsub('map([0-9]*)', '\\1', id)
    query_url <- paste0('https://www.kegg.jp/kegg-bin/view_ortholog_table?map=',id)
  }
  cat('Pulling module table...\n')
  site <- read_html(query_url)
  cat('Reading headers...\n')
  fulltable <- as.data.frame(html_table(site)[[1]])[1:2,]
  fulltable[2,] <- gsub('(K[0-9]*).*', '\\1', fulltable[2,])

  cells <- html_elements(site, 'tr')

  fulldf <- data.frame(matrix(NA, nrow=length(cells), ncol=ncol(fulltable)))
  fulldf[1:2,] <- fulltable
  cat('Reading cells...\n')
  pb <- txtProgressBar(max=nrow(fulldf), style=3)
  for (i in 3:nrow(fulldf)){
    cur <- html_elements(cells[i], 'td')
    outrow <- rep(NA_character_, ncol(fulldf))
    outrow[1:3] <- html_text2(cur[1:3])
    for ( j in 4:length(cur) ){
      entry <- as.character(cur[[j]])
      v <- gsub('.*<a href="(.*)" .*', 'https://www.kegg.jp\\1', entry)
      outrow[j] <- ifelse(v=='<td></td>', NA_character_, v)
    }
    if (grepl('M .*', outrow[3]))
      outrow[3] <- gsub('M (.*)', '\\1', outrow[3])
    fulldf[i,] <- outrow
    setTxtProgressBar(pb, i)
  }
  cat('\n')
  return(fulldf)
}

find_assembly_for_name <- function(to_search){
  link <- paste0('https://rest.kegg.jp/find/genome/', to_search)

  r <- retry_until_fail(link)
  if(is.null(r)) return(r)

  ex <- gsub('(.*)\\t(.*);.*', '\\1 \\2', strsplit(content(r), '\n')[[1]])
  res <- ''
  for (entry in ex){
    e <- strsplit(entry, ' ')[[1]]
    if(e[2] == to_search){
      res <- e[1]
      break
    }
  }

  genomelink <- paste0('https://rest.kegg.jp/get/', res)
  r <- retry_until_fail(genomelink)
  if(is.null(r)) return(NULL)
  r <- content(r)
  asn <- NULL
  if (grepl('Assembly', r))
    asn <- gsub(".*Assembly:([0-9A-Za-z_.]*).*", '\\1', r)
  return(asn)
}

find_genes_for_ko <- function(entry, verbose=TRUE){
  r <- get_from_id(entry, stopWithError = FALSE)
  stop_for_status(r)
  if (verbose){
    message_for_status(r)
    cat('\n')
  }
  r <- content(r, 'text')
  r <- strsplit(r, '\n')[[1]]
  genestart <- which(substr(r,1, 5) == 'GENES')
  r <- r[(genestart):length(r)]
  genend <- which(substr(r, 1, 1)!= ' ')
  genend <- ifelse(length(genend)==1, length(r)+1, genend[2])
  r <- r[1:(genend-1)]

  genes <- gsub("(GENES)? *([A-Z]{2,4}): ([A-Z0-9a-z_ ]*).*", '\\2:\\3', r)
  extras <- c()
  for (i in seq_along(genes)){
    gene <- genes[i]
    if(grepl(' ', gene)){
      splitgene <- strsplit(gene, ' ')[[1]]
      tmp <- strsplit(splitgene[1], ':')[[1]]
      genes[i] <- splitgene[1]
      splitgene <- splitgene[-1]
      id <- tmp[1]
      splitgene <- vapply(splitgene, \(x) paste0(id, ':', x),
                          character(1), USE.NAMES=FALSE)
      extras <- c(extras, splitgene)
    }
  }
  genes <- c(genes, extras)
  names(genes) <- NULL
  genes <- vapply(genes, \(x){
    spg <- strsplit(x, ':')[[1]]
    paste0(tolower(spg[1]), ':', spg[2])
  }, character(1), USE.NAMES=FALSE)

  return(genes)
}

get_ko_description <- function(entry, verbose=TRUE){
  r <- get_from_id(entry, stopWithError = FALSE)
  stop_for_status(r)
  if (verbose){
    message_for_status(r)
    cat('\n')
  }
  r <- content(r, 'text')
  r <- strsplit(r, '\n')[[1]]

  symline <- r[grepl("SYMBOL", r)][1]
  symbol <- gsub('SYMBOL *([^ ].*)', '\\1', symline)

  defnline <- r[grepl('NAME', r)][1]
  defn <- gsub('NAME *([^ ].*)', '\\1', defnline)

  retval <- c(symbol, defn)
  names(retval) <- c("Symbol", "Definition")
  return(retval)
}

get_ko_modpaths <- function(entry, verbose=TRUE){
  r <- get_from_id(entry, stopWithError = FALSE)
  stop_for_status(r)
  if (verbose){
    message_for_status(r)
    cat('\n')
  }
  r <- content(r, 'text')
  #r <- strsplit(r, '\n')[[1]]

  allmodules <- character(0L)
  allpaths <- character(0L)

  if(grepl("MODULE", r)){
    allmodules <- unlist(regmatches(r,
                                    gregexpr(" (M[0-9]{5}) ", r)))
    allmodules <- substring(allmodules, first=2, last=nchar(allmodules)-1)
  }
  if(grepl("PATHWAY", r)){
    allpaths <- unlist(regmatches(r,
                                    gregexpr(" map[0-9]{5} ", r)))
    allpaths <- substring(allpaths, first=2, last=nchar(allpaths)-1)
  }

  retval <- list(Pathways=allpaths, Modules=allmodules)
  return(retval)
}

get_seq_for_dbentry <- function(entry, position=FALSE, uselink=FALSE){
  if (uselink){
    require(httr)
    entry <- gsub('https://www.kegg.jp/entry(.*)', 'http://rest.kegg.jp/get\\1', entry)
    r <- retry_until_fail(entry)
    r <- content(r, 'text')
  } else {
    r <- content(get_from_id(entry), 'text')
  }

  retval <- list()
  AAseq <- NULL
  NTseq <- NULL
  if (grepl('AASEQ', r)){
    AAseq <- gsub('.*AASEQ(.*)NTSEQ.*', '\\1', r)
    AAseq <- paste0(strsplit(AAseq, '\n[ ]*')[[1]][-1], collapse='')
    AAseq <- AAString(AAseq)
  }
  if (grepl('NTSEQ', r)){
    NTseq <- gsub('.*NTSEQ +(.*)', '\\1', r)
    NTseq <- paste0(strsplit(NTseq, '\n[ ///]*')[[1]][-1], collapse='')
    NTseq <- DNAString(NTseq)
  }
  retval$AA <- AAseq
  retval$NT <- NTseq
  if (position && grepl('POSITION', r)){
    if (grepl('POSITION +Un', r)){
      pos <- NULL
    } else if (grepl('POSITION +[0-9A-Z:]*complement', r)){
      pos <- gsub('.*POSITION +([0-9A-Z:.]*complement\\([0-9:.]*).*', '\\1', r)
      pos <- gsub("(.*)complement\\((.*)", 'complement \\1\\2', pos)
    } else {
      pos <- gsub('.*POSITION +([0-9A-Z:.]*).*', '\\1', r)
    }
    retval$pos <- pos
  }
  retval$errored <- FALSE
  return(retval)
}

seqset_for_cog <- function(acctable, colnum, AA=T){
  orgs <- acctable[, colnum]
  l <- length(orgs)

  if(AA){
    check <- 'AAStringSet'
    sets <- AAStringSet(character(l))
  }
  else {
    check <- 'DNAStringSet'
    sets <- DNAStringSet(character(l))
  }

  for (i in seq_len(l)){
    print(paste('Seq', i, 'of', l))
    v <- orgs[i]
    vals <- get_seq_for_dbentry(v, AA)
    if (!is.null(vals$seq))
      sets[[i]] <- vals$seq
    if (!is.null(vals$errored) && vals$errored)
      orgs[i] <- paste0("ERROR ", v)
    #stopifnot(is(vals, check))
    vals <- NA
  }

  names(sets) <- orgs
  return(sets)
}

get_multiple_seqs_from_ids <- function(ids, AA=FALSE, partition=10L, verbose=TRUE){
  tf <- tempfile(tmpdir = file.path(getwd(), "temp"))
  num_iter <- length(ids) %/% partition + ifelse(length(ids) %% partition == 0, 0, 1)
  if (verbose) pb <- txtProgressBar(style=3,
                                    max=num_iter)
  cantfind <- c()
  for (i in seq_len(num_iter)){
    s <- (i-1)*partition + 1
    end <- min(i*partition, length(ids))
    ss <- ids[s:end]
    req <- paste(BASE_STRING, 'get',
                 paste0(ss, collapse='+'),
                 ifelse(AA, 'aaseq', 'ntseq'), sep='/')
    r <- retry_until_fail(req)
    if (!is.null(r)){
      r <- content(r, 'text')
      seqs_only <- strsplit(r, '>')[[1]]
      seqs_only <- paste0(seqs_only[2:length(seqs_only)], collapse='\n>')
      seqs_only <- paste0('>', seqs_only)
      cat(seqs_only, file=tf, append=TRUE)
    } else {
      cantfind <- c(cantfind, ss)
    }
    if(verbose) setTxtProgressBar(pb, i)
  }

  if(length(cantfind > 0)){
    warning("Possible error evaluating the following sets: \n- ", paste(cantfind, collapse='\n- '))
  }
  if (file.exists(tf)){
    if (AA){
      seqs <- readAAStringSet(tf, format='fasta')
    } else {
      seqs <- readDNAStringSet(tf, format='fasta')
    }
  } else {
    if (AA){
      seqs <- AAStringSet()
    } else {
      seqs <- DNAStringSet()
    }
  }
  return(seqs)
}

seqset_for_dbcog <- function(entries, pos=TRUE){
  if (grepl("M[0-9]*", entries[2])){
    entries <- entries[-2]
    entries[2] <- gsub('.*(K[0-9]*).*', '\\1', entries[2])
  }
  entries <- entries[-(1:2)] # remove header
  entries <- entries[!is.na(entries)]
  if (length(entries) < 1){
    return(list())
  }
  nolink <- gsub('https://www.kegg.jp/entry/(.*)', '\\1', entries)
  orgs <- gsub('(.*):.*', '\\1', nolink)
  l <- length(entries)
  retval <- list()
  AAset <- AAStringSet(character(l))
  NTset <- DNAStringSet(character(l))
  if (pos){
    positions <- data.frame(Organism=orgs, ID=nolink, link=entries,
                            chromosome=rep(NA_integer_, l),
                            start=rep(NA_integer_, l),
                            end=rep(NA_integer_, l))
  }
  pb <- txtProgressBar(max=l, style=3)
  for (i in seq_len(l)){
    v <- entries[i]
    vals <- get_seq_for_dbentry(v, position=pos, uselink=TRUE)
    if (!is.null(vals$AA))
      AAset[[i]] <- vals$AA
    if (!is.null(vals$NT))
      NTset[[i]] <- vals$NT
    if (pos){
      posvals <- vals$pos
      outpos <- rep(NA_integer_, 3)
      if (!is.null(posvals)){
        COMP <- FALSE
        if (grepl('complement', posvals)){
          posvals <- gsub('complement (.*)', '\\1', posvals)
          COMP <- TRUE
        }
        if (grepl(':', posvals)){
          posvals <- gsub('(\\d*):(\\d*)\\.\\.(\\d*)', '\\1 \\2 \\3', posvals)
          posvals <- as.integer(strsplit(posvals, ' ')[[1]])
        } else {
          posvals <- gsub('(\\d*)\\.\\.(\\d*)', '\\1 \\2', posvals)
          posvals <- c(1, as.integer(strsplit(posvals, ' ')[[1]]))
        }
        if (COMP){
          posvals[1] <- posvals[1] * -1
        }
      } else {
        posvals <- rep(NA_integer_, 3)
      }
      outpos[1:length(posvals)] <- posvals
      positions[i,4:6] <- outpos
    }
    #stopifnot(is(vals, check))
    vals <- NULL
    setTxtProgressBar(pb, i)
  }

  names(AAset) <- orgs
  names(NTset) <- orgs
  retval$AA <- AAset
  retval$NT <- NTset
  if (pos) retval$positions <- positions
  return(retval)
}

pathset_for_cog <- function(cogidx, cogs, silent=F){
  cog <- COGs[[cogidx]]
  ids <- which(IDs %in% cog)
  annots <- unique(Annotations[ids])
  annots <- annots[annots!='PNACT']
  if (length(annots) == 0){
    if (!silent)
      warning('No identified functions (COG is only PNACT)')
    return(c('PNACT'))
  }

  kgroups <- gsub('(K[0-9A-Z]*) .*', '\\1', annots)

  pathways <- vector('list', length(kgroups))
  for (i in seq_along(kgroups)){
    k <- kgroups[i]
    r <- get_from_id(k)
    r <- content(r, 'text')

    if (grepl('PATHWAY', r)){
      entry <- regmatches(r, gregexpr('map[0-9]{5}', r))[[1]]
    } else {
      entry <- character(0)
    }

    pathways[[i]] <- entry
  }

  return(list(paths=unique(unlist(pathways)), annotations=annots))
}

build_paths_set <- function(){
  cogpaths <- vector('list', length(COGs))
  for (i in seq_along(COGs)){
    cat('COG', i, 'of', length(COGs), '\r')
    p <- pathset_for_cog(i, T)
    if (!is.null(p))
      cogpaths[[i]] <- p

    if (i %% 100 == 0)
      save(cogpaths, file='/Users/aidan/Nextcloud/RStudioSync/streptomyces_data/COG_KEGG_Pathways.RData')
  }
  cat('Done.              \n')
  save(cogpaths, file='/Users/aidan/Nextcloud/RStudioSync/streptomyces_data/COG_KEGG_Pathways.RData')
  return(cogpaths)
}

genome_length_from_link <- function(name){
  require(httr)
  entry <- paste0('http://rest.kegg.jp/get/', name)
  r <- GET(entry)
  ctr <- 0
  while (r$status_code == 403){
    if (ctr == 10) stop('Retried 10 times with failures.')
    ctr <- ctr + 1
    Sys.sleep(30)
    r <- GET(entry)
  }
  #stop_for_status(r)
  if (r$status_code==404){
    return(list(errored=TRUE))
  }
  r <- content(r, 'text')
}

ko_for_gene <- function(entry){
  r <- get_from_id(entry, stopWithError = FALSE)
  if (r$status_code==404){
    return(NULL)
  }
  r <- content(r, 'text')
  if(!grepl("ORTHOLOGY", r)){
    return(NULL)
  }

  line <- gsub('.*ORTHOLOGY *(K[0-9]{5}).*', '\\1', r)
  return(line)
}

### Beast of a function to parse out KEGG's modules into pairs
#
# - Complicated
# - Run get_path with a module name
# -   Example: get_path('M00001')

# - Default return value is pairs of mutual connections
# - asAdjMatrix=T returns adjacency matrix instead of pairs
get_module_path <- function(mod, asAdjMatrix=F, verbose=F){
  split_str <- function(modstr, useCommas=F){
    if (useCommas) splitChar <- ','
    else splitChar <- ' '
    chars <- strsplit(modstr, '')[[1]]
    vals <- rep(0, length(chars))
    vals[chars=='('] <- 1
    vals[chars==')'] <- -1
    cs <- cumsum(vals)

    split_points <- c(0, intersect((which(cs==0) + 1), which(chars==splitChar)), length(chars)+1)
    subs <- character(length(split_points)-1)
    i <- 1
    for (i in seq_along(subs)){
      subs[i] <- paste(chars[(split_points[i]+1):(split_points[i+1]-1)], collapse='')
    }

    return(subs)
  }

  deparse_module <- function(modstring, AdjMatrix=F){
    # remove newlines
    modstring <- gsub('\\n', ' ', modstring)

    # I'm not really sure how to handle +, this seems like a good approach
    #modstring <- gsub('((K[0-9]{5}\\+)+K[0-9]{5})', '\\(\\1\\)', modstring)

    modstring <- gsub(' *\\+ *', ' ', modstring)
    modstring <- gsub(' *\\- *', ' ', modstring)
    # initialize a list
    nodes <- list()
    nodes[[1]] <- modstring

    incoming <- outgoing <- list()
    incoming[[1]] <- outgoing[[1]] <- 0

    MAX_ITER <- 100
    k <- 0
    while ( k < MAX_ITER && any(sapply(nodes, function(x) nchar(x) > 6)) ){
      k <- k + 1
      newnodes <- nodes
      newinc <- incoming
      newout <- outgoing
      for ( i in seq_along(nodes) ){
        space_split <- split_str(nodes[[i]], useCommas=F)

        for ( j in seq_along(space_split) ){
          if ( j == 1 ){
            newnodes[[i]] <- space_split[[j]]
            if (length(space_split) > 1){
              newout[[i]] <- length(newnodes)+1
            }
          } else {
            if ( j > 2 )
              newout[[length(newout)]] <- length(newout) + 1
            newinc <- append(newinc, length(newnodes))
            newnodes <- append(newnodes, space_split[[j]])
            newout[[length(newnodes)]] <- outgoing[[i]]
          }
        }

        if (length(space_split) > 1){
          idxs_to_change <- which(sapply(newinc, function(x) i %in% x))
          idxs_to_change <- idxs_to_change[idxs_to_change != (length(newnodes)-length(space_split)+2)]
          for (idx in idxs_to_change){
            v <- newinc[[idx]]
            newinc[[idx]] <- c(v[v!=i], length(newnodes))
          }
        }
      }
      nodes <- newnodes
      outgoing <- newout
      incoming <- newinc

      nodes <- lapply(nodes, function(x){
        dp <- strsplit(x, '')[[1]]
        if (dp[1] == '(' && dp[length(dp)] == ')'){
          dp <- dp[-c(1, length(dp))]
        }
        return(paste(dp, collapse=''))
      })


      for ( i in seq_along(nodes) ){
        comm_split <- split_str(nodes[[i]], useCommas=T)
        prev_length <- length(newnodes)
        for ( j in seq_along(comm_split) ){
          if ( j == 1 ){
            newnodes[[i]] <- comm_split[[j]]
            append_to_all <- F
          } else {
            append_to_all <- T
            newnodes <- append(newnodes, comm_split[[j]])
            newinc <- append(newinc, incoming[[i]])
            newout <- append(newout, outgoing[[i]])
          }
        }
        idxs_to_change <- sapply(newout, function(x) (i %in% x))
        if ( append_to_all & any(idxs_to_change) ){
          for ( val in which(idxs_to_change) )
            newout[[val]] <- c(newout[[val]], (prev_length+1):(prev_length+j-1))
        }
        idxs_to_change <- sapply(newinc, function(x) (i %in% x))
        if ( append_to_all & any(idxs_to_change) ){
          for ( val in which(idxs_to_change) )
            newinc[[val]] <- c(newinc[[val]], (prev_length+1):(prev_length+j-1))
        }
      }

      nodes <- newnodes
      outgoing <- newout
      incoming <- newinc
    }

    if (k == MAX_ITER){
      stop(paste0('Function hanging. Terminating to prevent infinite loop.\n \
                   \rProblematic string: \n', modstring))
    }

    nodes <- unlist(nodes)


    outgoing <- lapply(outgoing, function(x) nodes[x[x!=0]])
    incoming <- lapply(incoming, function(x) nodes[x[x!=0]])

    names(outgoing) <- names(incoming) <- nodes

    if (AdjMatrix){
      A <- matrix(0, ncol=length(nodes), nrow=length(nodes))
      rownames(A) <- colnames(A) <- nodes
      for (node in nodes){
        A[outgoing[[node]], node] <- 1
        A[node,incoming[[node]]] <- 1
      }
      return(A)
    } else {
      pairs <- NULL
      for (node in nodes){
        m <- rbind(expand.grid(node, incoming[[node]]), expand.grid(node, outgoing[[node]]))
        pairs <- rbind(pairs, m)
      }
      colnames(pairs) <- c("KO_1", 'KO_2')
      pairs <- as.matrix(pairs)
      if (nrow(pairs) == 0)
        return(NULL)
      for (i in 1:nrow(pairs)){
        if (pairs[i,1] < pairs[i,2]){
          t <- pairs[i,1]
          pairs[i,1] <- pairs[i,2]
          pairs[i,2] <- t
        }
      }
      return(unique(pairs))
    }
  }

  mod <- trimws(mod)
  stopifnot('Malformed Pattern'=grepl('^M[0-9]{5}$', mod))
  val <- content(get_from_id(mod))
  if (grepl('DEFINITION', val)){
    charstring <- gsub('.*DEFINITION  ([K0-9()+ ,-\\n]*).*ORTHOLOGY.*', '\\1', val)
    charstring <- gsub('( *\\-\\- *)+', ' ', charstring)
    charstring <- trimws(charstring)
    elems <- split_str(charstring, T)
    elems <- sapply(elems, function(x){
      s <- strsplit(x, '')[[1]]
      v <- rep(0, length(s))
      v[s=='('] <- 1
      v[s==')'] <- -1
      cs <- cumsum(v)
      if (length(which(cs==0)) == 0){
        to_move_in <- min(cs)
        return(paste(s[(1+to_move_in):(length(s)-to_move_in)]), collapse=T)
      }
      else
        return(x)
    }, USE.NAMES=F)

    if (!asAdjMatrix){
      allrows <- NULL
      for (elem in elems)
        allrows <- rbind(allrows, deparse_module(elem))
      if (!is.null(allrows)){
        if (verbose) cat(nrow(allrows), 'Total Pairs \n')
        return(allrows)
      } else {
        if (verbose) cat('No pairs found\n')
        return(matrix(nrow=0, ncol=2))
      }
    }
    else {
      outlst <- vector('list', length=length(elems))
      for (i in seq_along(elems))
        outlst[[i]] <- deparse_module(elems[i], T)
      return(outlst)
    }
  }
}
