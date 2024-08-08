load(file.path(datadir, "Modules", "ModulePredsAllPairs.RData"))
load(file.path(datadir, "Multiclass", 'MulticlassModuleData.RData'))

#cols <- c("#530FAB", "#83F1B0", "#DAB29D")
cols <- rep("#DAB29D", 3L)
USE_WEIGHT <- FALSE

outfile <- "6_FigCaseStudy.pdf"
N_STUDIES <- 6
N_GRAPHS <- 6
pdf(file=file.path(figdir, "MainFigures", outfile), onefile = TRUE, height=7.5, width=6.08)
#layout(matrix(seq_len(N_STUDIES*N_GRAPHS), byrow=TRUE, nrow=N_STUDIES))
layout_m <- matrix(0, nrow=N_STUDIES*2, ncol=(N_GRAPHS-4L)*2+2)
for(i in seq_len(N_STUDIES)){
  layout_m[2*i-1,1:2] <- (i-1)*6 + 1:2
  layout_m[2*i,1:2] <- (i-1)*6 + 3:4
  layout_m[c(2*i-1,2*i), 3:4] <- (i-1)*6 + 5
  layout_m[c(2*i-1,2*i), 5:6] <- (i-1)*6 + 6
}
#layout(layout_m, heights=rep(1,nrow(layout_m)-1, 0.01))
layout(layout_m, widths=c(1,1,0.75,0.75,1,1.1))
#layout(layout_m)

if(USE_WEIGHT){
  BuildCaseStudyForInputNetwork <- BuildCaseStudyForInputNetworkWeight
} else {
  BuildCaseStudyForInputNetwork <- BuildCaseStudyForInputNetworkNoWeight
}
par(mai=c(0,0,0,0), mgp=c(0,0,0), oma=c(3,3,3,0))

find_all_blocks <- function(to_find){
  p1 <- unlist(lapply(to_find, \(x) which(grepl(x, AllPairs$Mod1Orig, fixed=TRUE))))
  p2 <- unlist(lapply(to_find, \(x) which(grepl(x, AllPairs$Mod2Orig, fixed=TRUE))))
  subp <- AllPairs[unique(p1), c(1,19)]
  subp2 <- AllPairs[unique(p2), c(2,20)]
  colnames(subp) <- colnames(subp2) <- c("Encoded", "Original")
  subp <- rbind(subp, subp2)
  subp <- subp[!duplicated(subp[,1]),]
  mapping <- vapply(subp[,2], \(x) paste(sort(x), collapse=','), character(1L))
  names(mapping) <- subp[,1]
  mapping <- sort(mapping)
  kegg_blocks <- unique(mapping)
  gene_names <- character(length(kegg_blocks))
  names(gene_names) <- kegg_blocks
  names(kegg_blocks) <- NULL
  encoded_blocks <- names(mapping)
  list(EncodeToKEGG=mapping, KEGGToGene=gene_names,
       Encoded=encoded_blocks, Blocks=kegg_blocks)
}

reorder_all_blocks <- function(ModulesObject, neworder){
  ModulesObject$KEGGToGene <- ModulesObject$KEGGToGene[neworder]
  ModulesObject$Blocks <- ModulesObject$Block[neworder]
  ModulesObject
}

## Lysine Biosynthesis
.plotLysine <- function(plotHeader=FALSE){
  cat("\tPlotting Lysine case study...\n")
  modules_of_interest <- c("M00030", "M00031")
  all_blocks <- find_all_blocks(modules_of_interest)
  all_blocks
  all_blocks$KEGGToGene[] <- c(
    "LYS21", "ACO2", "LYS4", "ARO8", "LYS2", "LYS9", "LYS1",
    "LysX", "LysZ", "LysY", "LysJ", "LysK"
  )
  fills <- rep(cols[1:2], times=c(7,5))
  names(fills) <- all_blocks$KEGGToGene
  ActualGraph <- matrix(c(
    1:6, 8:11, 4, 7,
    2:7, 9:12, 8, 12
  ), ncol=2)

  title <- "Lysine\nBiosynthesis"
  BuildCaseStudyForInputNetwork(AllPairs, subpreds, allpredictions,
                                all_blocks, fills,
                                ActualGraph, title, plotHeader)
}

## Penicillin Biosynthesis
.plotPenicillin <- function(plotHeader=FALSE){
  cat("\tPlotting Penicillin case study...\n")
  modules_of_interest <- c("M00672", "M00673")
  all_blocks <- find_all_blocks(modules_of_interest)
  all_blocks
  all_blocks$KEGGToGene[] <- c(
    "PCBAB", "PCBC", "PENDE",
    "CefD", "CefE", "CefF", "CmcH", "CmcI", "CmcJ"
  )
  fills <- rep(cols[1:2], times=c(3,6))
  names(fills) <- all_blocks$KEGGToGene
  ActualGraph <- matrix(c(
    1:2, 2,4:8,
    2:3, 4,5:9
  ), ncol=2)

  title <- "Penicillin &\nCephalosporin\nBiosynthesis"
  BuildCaseStudyForInputNetwork(AllPairs, subpreds, allpredictions,
                                all_blocks, fills,
                                ActualGraph, title, plotHeader)
}

## Acarbose and Validamycin Biosynthesis
.plotAcarboseValidamycin <- function(plotHeader=FALSE){
  cat("\tPlotting Acarbose/Validamycin case study...\n")
  modules_of_interest <- c("M00814", "M00815")
  all_blocks <- find_all_blocks(modules_of_interest)
  all_blocks
  all_blocks$KEGGToGene[] <- c(
    "AcbC","AcbM","AcbO","AcbL","AcbN","AcbU","AcbR","AcbS",
    "CetB","ValK","ValC","ValN","ValM","ValL","VldH","ValG"
  )
  fills <- rep(cols[1:2], times=c(8,8))
  names(fills) <- all_blocks$KEGGToGene
  ActualGraph <- matrix(c(
    1:7, 1, 9:15, 11, 11, 7,
    2:8, 9, 10:16, 6, 7, 14
  ), ncol=2)

  title <- "Acarbose &\nValidamycin\nBiosynthesis"
  BuildCaseStudyForInputNetwork(AllPairs, subpreds, allpredictions,
                                all_blocks, fills,
                                ActualGraph, title, plotHeader)
}

## Acarbose Biosynthesis
.plotAcarbose <- function(plotHeader=FALSE){
  cat("\tPlotting Acarbose case study...\n")
  modules_of_interest <- c("M00814")
  all_blocks <- find_all_blocks(modules_of_interest)
  all_blocks
  all_blocks$KEGGToGene[] <- c(
    "AcbC","AcbM","AcbO","AcbL","AcbN","AcbU","AcbR","AcbS"
  )
  fills <- rep(cols[1:2], times=c(8,8))
  names(fills) <- all_blocks$KEGGToGene
  ActualGraph <- matrix(c(
    1:7,
    2:8
  ), ncol=2)

  title <- "Acarbose\nBiosynthesis"
  BuildCaseStudyForInputNetwork(AllPairs, subpreds, allpredictions,
                                all_blocks, fills,
                                ActualGraph, title, plotHeader)
}

## Urea Cycle and Ornithine Biosynthesis
.plotUreaOrnithine <- function(plotHeader=FALSE){
  cat("\tPlotting Urea/Ornithine case study...\n")
  modules_of_interest <- c("M00763", "M00029")
  all_blocks <- find_all_blocks(modules_of_interest)
  all_blocks <- reorder_all_blocks(all_blocks, c(1:5,10,6:9))
  all_blocks
  all_blocks$KEGGToGene[] <- c(
    "CPS1", "OTC", "ArgG", "ArgH", "RocF",
    "ArgX", "ArgB", "ArgC", "ArgD", "ArgE"
  )
  fills <- rep(cols[1:2], times=c(5,5))
  names(fills) <- all_blocks$KEGGToGene
  ActualGraph <- matrix(c(
    1:4,5, 6:10,   6,
    2:5,2, 7:10,2, 1
  ), ncol=2)
  title <- "Urea Cycle &\nOrnithine\nBiosynthesis"
  BuildCaseStudyForInputNetwork(AllPairs, subpreds, allpredictions,
                                all_blocks, fills,
                                ActualGraph, title, plotHeader)
}

## Biotin metabolism
.plotBiotin <- function(plotHeader=FALSE){
  cat("\tPlotting Biotin case study...\n")
  modules_of_interest <- c("M00950", "M00573", "M00572")
  all_blocks <- find_all_blocks(modules_of_interest)
  all_blocks
  all_blocks$KEGGToGene[] <- c(
    "BioC", "FabBF", "FabG", "FabZ", "FabI", "BioH", #572
    "BioI", "BioF", "BioA", "BioD", "BioB", # 573
    "BioU" # 950
  )
  fills <- rep(cols, times=c(6,5,1))
  names(fills) <- all_blocks$KEGGToGene
  ActualGraph <- matrix(c(
    1:5,5, 7:10, 6, 8, 12,
    2:6,2, 8:11, 8, 12, 10
  ), ncol=2)
  title <- "Biotin\nBiosynthesis"
  BuildCaseStudyForInputNetwork(AllPairs, subpreds, allpredictions,
                                all_blocks, fills,
                                ActualGraph, title, plotHeader)
}

## Ubiquinone biosynthesis
.plotUbiquinone <- function(plotHeader=FALSE){
  cat("\tPlotting Ubiquinone case study...\n")
  modules_of_interest <- c("M00117", "M00128")
  all_blocks <- find_all_blocks(modules_of_interest)
  reorder_all_blocks(all_blocks, c(2:6,1,7:12))
  all_blocks
  all_blocks$KEGGToGene[] <- c(
    "UbiC", "UbiA", "UbiI", "UbiG", "UbiH", "UbiE", "UbiF",
    "COQ2", "COQ6", "COQ3", "COQ5", "COQ7"  # 128
  )
  fills <- rep(cols, times=c(7,5,0))
  names(fills) <- all_blocks$KEGGToGene
  ActualGraph <- matrix(c(
    1:6,7, 8:11,12, #2, 10, 5:7,
    2:7,4, 9:12,10 #9, 5,  11:12,10
  ), ncol=2)
  title <- "Ubiquinone\nBiosynthesis"
  BuildCaseStudyForInputNetwork(AllPairs, subpreds, allpredictions,
                                all_blocks, fills,
                                ActualGraph, title, plotHeader)
}

.plotCymeneCumate <- function(plotHeader=FALSE){
  cat("\tPlotting Cymene/Cumate case study...\n")
  modules_of_interest <- c("M00419", "M00539")
  all_blocks <- find_all_blocks(modules_of_interest)
  all_blocks
  all_blocks$KEGGToGene[] <- c(
    "CymB", "CymC", "CmtB", "CmtC", "CmtD", "CmtE"
  )
  fills <- rep(cols, times=c(2,4,0))
  names(fills) <- all_blocks$KEGGToGene
  ActualGraph <- matrix(c(
    1:5,
    2:6
  ), ncol=2)
  title <- "Cymene & Cumate\nDegradation"
  BuildCaseStudyForInputNetwork(AllPairs, subpreds, allpredictions,
                                all_blocks, fills,
                                ActualGraph, title, plotHeader)
}

.plotPectin <- function(plotHeader=FALSE){
  cat("\tPlotting Pectin case study...\n")
  modules_of_interest <- c("M00081", "M00061")
  all_blocks <- find_all_blocks(modules_of_interest)
  all_blocks <- reorder_all_blocks(all_blocks, c(6:8,1:5))
  all_blocks
  all_blocks$KEGGToGene[] <- c(
    "E3.1.1.11", "E3.2.1.15", "E3.2.1.67",
    "uxaC", "uxaB", "uxaA", "kdgK", "dgaF"
  )
  fills <- rep(cols, times=c(3,5,0))
  names(fills) <- all_blocks$KEGGToGene
  ActualGraph <- matrix(c(
    1:7,
    2:8
  ), ncol=2)
  title <- "Pectin\nDegradation"
  BuildCaseStudyForInputNetwork(AllPairs, subpreds, allpredictions,
                                all_blocks, fills,
                                ActualGraph, title, plotHeader)
}

.plotAcarbose(TRUE)
.plotCymeneCumate()
.plotUreaOrnithine()
.plotPectin()
.plotBiotin()
.plotPenicillin()
#.plotUbiquinone()
#.plotLysine()

cat("\tPlotting legend...\n")
cols_legend <- c('#D81B60', '#E0A608', '#2B6DA8', '#45A649', 'grey40', 'black')
## legend for component algorithms
LEGEND_X <- -8.75
LEGEND_Y <- -1.0
LEGEND_CEX <- 1
legend(x=LEGEND_X, y=LEGEND_Y, cex=LEGEND_CEX, x.intersp=c(1,0.25,0.75,-1.5),
       legend=c("Phylogenetic Profiling", "Phylogenetic Structure",
                "Gene Organization", "Sequence Level"),
       bty='n', xpd=NA,
       text.col=cols_legend[c(1,3,2,4)], text.font=2,
       horiz=TRUE)
legend(x=LEGEND_X+0.225, y=LEGEND_Y-0.35, horiz=TRUE, cex=LEGEND_CEX, text.font=2,
       x.intersp=0.5,
       legend="Edge Differs from KEGG",
       col=cols_legend[6], lty=2, lwd=1.5, bty='n', xpd=NA, seg.len=2)
legend(x=LEGEND_X+3.75, y=LEGEND_Y-0.35, horiz=TRUE, cex=LEGEND_CEX, text.font=2,
       x.intersp=0.5,
       legend="Edge Matches KEGG",
       col=cols_legend[6], lty=1, lwd=1.5, bty='n', xpd=NA, seg.len=2)
legend(x=LEGEND_X+6.975, y=LEGEND_Y-0.35, horiz=TRUE, cex=LEGEND_CEX, text.font=2,
       x.intersp=0.5,
       legend="Direction in KEGG",
       col=cols_legend[6], lty=1, lwd=1.5, bty='n', xpd=NA, seg.len=2)
polygon(x=rep(LEGEND_X+7.485,3)+c(0,0,0.075),
        y=rep(LEGEND_Y-0.605,3)+c(0.04,-0.04,0),
        col='black', xpd=NA)
dev.off(dev.list()['pdf'])
cat("\tDone!\n")
