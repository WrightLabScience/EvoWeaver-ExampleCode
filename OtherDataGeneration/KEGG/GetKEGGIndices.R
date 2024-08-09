library(httr)

get_all_cds <- function(orgname){
	url_path <- paste0("rest.kegg.jp/list/", orgname)
	r <- GET(url_path)
	#tf <- tempfile()
	#writeLines(content(r), con=tf)
	tab <- read.delim(text=content(r), header=FALSE)
	colnames(tab) <- c("id", "type", "position", "description")
	tab$chromosome <- tab$position
	tab$direction <- grepl("complement", tab$position)
	tab$position[tab$direction] <- gsub("complement\\(([0-9]+\\.\\.[0-9]+)\\)", '\\1', tab$position[tab$direction])
	if(!any(grepl(":", tab$chromosome))){
		tab$chromosome <- '1'
	} else {
		tab$chromosome <- vapply(strsplit(tab$chromosome, ':'), .subset, character(1L), 1L)
		tab$position <- vapply(strsplit(tab$position, ':'), .subset, character(1L), 2L)
	}
	tab$start <- vapply(strsplit(tab$position, '..', fixed=TRUE), .subset, character(1L), 1L)
	tab$end <- vapply(strsplit(tab$position, '..',fixed=TRUE), .subset, character(1L), 2L)
	tab$is_complement <- as.integer(tab$direction)
	tab <- tab[,c('id', 'type', 'chromosome', "is_complement", 'start', 'end', 'description')]
	tab
}

load('ModulesSpeciesTree.RData')
load('ModuleIndex_KO.RData', v=TRUE)

all_orgs <- labels(SpecTree)
ll <- length(all_orgs)
times <- rep(Sys.time(), 3L)
for(i in seq_along(all_orgs)){
	cat(i, '/', ll, '\r')
	Orgname <- all_orgs[i]
	fname <- file.path("KEGGIndices", paste0("AllCDS_", Orgname, ".RData"))
	if(file.exists(fname)) next
	cds_hits <- get_all_cds(Orgname)
	save(Orgname, cds_hits, file=fname)
	times[((i-1) %% 3) + 1] <- Sys.time()
	difference_time <- difftime(times[((i-3) %%3)+1], times[((i-1)%%3) + 1], units='secs')[[1]]
	if(difference_time < 1){
		Sys.sleep(1L)
	}
}

AllCDSs <- list()
for(i in seq_along(all_orgs)){
	cat(i, '/', ll, '\r')
	Orgname <- all_orgs[i]
	fname <- file.path("KEGGIndices", paste0("AllCDS_", Orgname, ".RData"))
	res <- NULL
	if(file.exists(fname)){
		load(fname)
		hits <- cds_hits
		hits <- hits[!apply(hits, 1, \(x) all(is.na(x))),]
		hits$start <- gsub("^[^0-9]*([0-9]+)[^0-9]?.*", '\\1', hits$start)
		hits$end <- gsub("^[^0-9]*([0-9]+)[^0-9]?.*", '\\1', hits$end)
	}
	AllCDSs[[Orgname]] <- hits
}

for(i in seq_along(AllCDSs)){
	cat(i, '/', ll, '\r')
	hits <- AllCDSs[[i]]
	hits$start[grepl(',', hits$start)] <- vapply(hits$start[grepl(',', hits$start)])
	hits$start <- gsub('[^0-9]', '', hits$start)
	hits$end <- gsub('^[0-9]', '', hits$end)
	AllCDSs[[i]] <- hits
}

save(AllCDSs, file='AllKEGGCDSs.RData')

# load(file='AllKEGGCDSs.RData')
KOgroup_dir <- 'KOGroups'
#KOgroup_dir <- 'KOGroups_noseqs'
outpath <- 'KOGroupsWithIndices'
HAS_SEQS <- TRUE
USE_CDS_ONLY <- TRUE
l <- list.files(KOgroup_dir)
for(i in seq_along(l)){
	cat(i, '/', length(l), '\n')
	if(file.exists(file.path(outpath, l[i]))) next
	load(file.path(KOgroup_dir, l[i]))
	if(HAS_SEQS) n <- names(seqs)
	orign <- n
	specs <- vapply(strsplit(orign, ":"), .subset, character(1L), 1)
	for(j in seq_along(specs)){
		cat('\t', j, '/', length(specs), '\r')
		hits <- AllCDSs[[specs[j]]]
		hits <- hits[hits$chromosome != "Unknown" & hits$start != "Unknown" & hits$end != "Unknown",]
#		hits <- hits[!apply(hits, 1, \(x) all(is.na(x))),]
		if(USE_CDS_ONLY){
			hits <- hits[hits$type=="CDS",]
		}
		id <- which(hits$id == orign[j])
		if(length(id) == 0){
			n[j] <- ""
			next
		}
		chrom <- hits[id,3]
		iscomp <- hits[id,4]
		pos <- as.integer(hits[id,5])
		hits <- hits[hits$chromosome==chrom,]
		hits$start <- as.integer(hits$start)
		pos <- sum(hits$start < pos, na.rm=TRUE)
		n[j] <- paste(specs[j], chrom, iscomp, pos, sep="_")
	}
	cat('\n')
	save(n, file=file.path(outpath, l[i]))
}
