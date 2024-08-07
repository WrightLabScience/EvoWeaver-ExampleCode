## This script was used to generate AllPairs data for Modules on distributed compute

## We will need the following:
##	- Pairings.RData: to get the job id
##	- EvoWeaver object
##	- inds object
##  - original species tree
##	- KEGG species tree
NUM_PAIRS <- 100L

PKG_NAME <- "SynExtend_1.15.4.tar.gz"
LIB <- "./libraries/"
dir.create(LIB)
Sys.setenv(TMPDIR=getwd())
install.packages(PKG_NAME, lib=LIB)
.libPaths(c(LIB, .libPaths()))
library(SynExtend)

Args <- commandArgs(trailingOnly=TRUE)
subval <- as.integer(Args[1L])
start_row <- subval*NUM_PAIRS + 1L
end_row <- (subval+1L) * NUM_PAIRS

load('Pairings.RData') # 0MB, Pairings
load("ModulesEvoWeaverObject.RData") # 650MB, pw -
load("ModsWithPositions.RData") # 15MB, inds
load('KEGGSpeciesTree.RData') # tree -
load("ModulesSpeciesTree.RData") # SpecTree -

end_row <- ifelse(end_row > nrow(Pairings), nrow(Pairings), end_row)

## normalize species tree heights
h1 <- attr(SpecTree, 'height')
h2 <- attr(tree, 'height')

SpecTree <- dendrapply(SpecTree, \(x){attr(x, "height") <- attr(x,'height')/h1;x})
tree <- dendrapply(tree, \(x){attr(x, "height") <- attr(x,'height')/h2;x})

## first create subpw objects to make processing faster
Subpairings <- Pairings[start_row:end_row,]
modules_used <- unique(unlist(Subpairings))
coloc_modules <- coloc_modules[modules_used]
subpw <- unclass(pw[modules_used])


algs_nocoloc <- c("PAJaccard", "PAOverlap", "GLDistance", "GLMI",
									"TreeDistance", "RPMirrorTree", "RPContextTree")
									#"SequenceInfo", "GeneVector")

algs_coloc <- c("GeneDistance", "MoransI", "OrientationMI")

n_total_algs <- length(algs_coloc) + length(algs_nocoloc)

pw_nocoloc_orig <- EvoWeaver(subpw, MySpeciesTree=SpecTree, NoWarn=TRUE)
pw_coloc_orig <- EvoWeaver(coloc_modules, MySpeciesTree=SpecTree, NoWarn=TRUE)
pw_nocoloc_kegg <- EvoWeaver(subpw, MySpeciesTree=tree, NoWarn=TRUE)
pw_coloc_kegg <- EvoWeaver(coloc_modules, MySpeciesTree=tree, NoWarn=TRUE)

print("Original Tree, no coloc")
set.seed(635L)
p1 <- predict(pw_nocoloc_orig, Method=algs_nocoloc, Verbose=TRUE, Subset=Subpairings)
print("Colocalization")
p2 <- predict(pw_coloc_orig, Method=algs_coloc, Verbose=TRUE, Subset=Subpairings, ReturnDataFrame=FALSE)
print("KEGG tree, no coloc")
set.seed(635L)
p3 <- predict(pw_nocoloc_kegg, Method=algs_nocoloc, Verbose=TRUE, Subset=Subpairings, ReturnDataFrame=FALSE)
print("KEGG tree, coloc")
p4 <- predict(pw_coloc_kegg, Method=algs_coloc, Verbose=TRUE, Subset=Subpairings, ReturnDataFrame=FALSE)

print("parsing results...")
results_orig <- results_kegg <- matrix(NA_real_, nrow=length(start_row:end_row), ncol=n_total_algs)
colnames(results_orig) <- colnames(results_kegg) <- c(algs_nocoloc, algs_coloc)
for(i in seq_along(algs_nocoloc)){
	preds_orig <- p1[[i]]
	preds_kegg <- p3[[i]]

	for(j in seq_len(nrow(Subpairings))){
		vs <- unlist(Subpairings[j,])
		results_orig[j,i] <- preds_orig[vs[1], vs[2]]
		results_kegg[j,i] <- preds_kegg[vs[1], vs[2]]
	}
}

offset <- length(algs_nocoloc)
for(i in seq_along(algs_coloc)){
	preds_orig <- p2[[i]]
	preds_kegg <- p4[[i]]

	for(j in seq_len(nrow(Subpairings))){
		vs <- unlist(Subpairings[j,])
		results_orig[i+offset,j] <- preds_orig[vs[1], vs[2]]
		results_kegg[i+offset,j] <- preds_kegg[vs[1], vs[2]]
	}
}

# p1 <- vapply(p1, \(x) x[[2]], numeric(1L))
# p2 <- vapply(p2, \(x) x[[2]], numeric(1L))
# p3 <- vapply(p3, \(x) x[[2]], numeric(1L))
# p4 <- vapply(p4, \(x) x[[2]], numeric(1L))

# p_orig <- c(p1,p2)
# p_kegg <- c(p3,p4)
# p_orig[is.na(p_orig)] <- 0
# p_kegg[is.na(p_kegg)] <- 0

# results <- lapply(1:2, \(x){v<-numeric(n_total_algs); names(v)<-c(algs_nocoloc,algs_coloc);v})
# names(results) <- c("Original", "KEGG")
# results$Original[] <- c(p1,p2)
# results$KEGG[] <- c(p3,p4)

# RETURNDATA <- list(rnum=subval, predictions=results)
RETURNDATA <- list(start=start_row, end=end_row, pairs=Subpairings, original=results_orig, kegg=results_kegg)
save(RETURNDATA, file=paste0("ModResult_", subval-1L, ".RData"))
