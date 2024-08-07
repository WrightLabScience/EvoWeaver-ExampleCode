load(file.path(datadir, "Modules", "ModulePredsAllPairs.RData"))
load(file.path(datadir, "Multiclass", 'MulticlassModuleData.RData'))

BuildCaseStudyForInputNetwork <- BuildCaseStudyForInputNetworkNoWeight
cols <- c("#530FAB", "#83F1B0", "#DAB29D")
pdf(file.path(figdir, "SupplFigures", "SXX_ExtendedFig5.pdf"),
    width=6.25,height=2.75, onefile=TRUE)
M_offset <- 0.25
N_STUDIES <- 2L
N_GRAPHS <- 6L
layout_m <- matrix(0, nrow=N_STUDIES*2, ncol=(N_GRAPHS-4L)*2+2)
for(i in seq_len(N_STUDIES)){
  layout_m[2*i-1,1:2] <- (i-1)*6 + 1:2
  layout_m[2*i,1:2] <- (i-1)*6 + 3:4
  layout_m[c(2*i-1,2*i), 3:4] <- (i-1)*6 + 5
  layout_m[c(2*i-1,2*i), 5:6] <- (i-1)*6 + 6
}
#layout(layout_m, heights=rep(1,nrow(layout_m)-1, 0.01))
par(mai=c(0,0,0,0), mgp=c(0,0,0), oma=c(3,2,2,1))
layout(layout_m)
modules_of_interest <- c("M00070", "M00071", "M00072", "M00073", "M00074", "M00075")
all_blocks <- find_all_blocks(modules_of_interest)
all_blocks <- reorder_all_blocks(all_blocks, c(1:10,12))
all_blocks
all_blocks$KEGGToGene[] <- c(
  "B3GNT5", "B3GALT", "B4GALT",
  "MOGS", "GANAB", "MAN1B", "MAN1A_C",
  "OCH1", "MNN2",
  "FUT8", "ST6GAL1"
)
#all_blocks$KEGGToGene[] <- gsub("M000(7[0-9].*)", "\\1", sapply(strsplit(all_blocks$Blocks, ','), .subset, 1))
fills <- rep(cols, times=c(3,8,0))
names(fills) <- all_blocks$KEGGToGene
ActualGraph <- matrix(c(
  1,1,4:9,10,3,
  2,3,5:10,3,11
), ncol=2)
title <- "ST6GAL/B3GNT5\nOther Connections"
BuildCaseStudyForInputNetwork(AllPairs, subpreds, allpredictions,
                              all_blocks, fills,
                              ActualGraph, title, TRUE, wrong_lty=3,
                              VLC=0.8, highlight_labels=c("B3GNT5", "ST6GAL1"),
                              TITLE_OFFSET=1, MARGIN_OFFSET=M_offset, VERTNAMESCENTER=TRUE)

modules_of_interest <- c("M00892", "M00055")
all_blocks <- find_all_blocks(modules_of_interest)
all_blocks <- reorder_all_blocks(all_blocks, c(1,6:11,2:5,13:17))
all_blocks
all_blocks$KEGGToGene[] <- c(
  "ALG7", "ALG1", "ALG2", "ALG11", "ALG3",
  "ALG9", "ALG12", "ALG5", "ALG6", "ALG8", "ALG10",
  "GPI", "GFPT", "GNPNAT1", "PGM3", "UAP1"
)
#all_blocks$KEGGToGene[] <- gsub("M00([0-9].*)", "\\1", sapply(strsplit(all_blocks$Blocks, ','), .subset, 1))
fills <- rep(cols, times=c(11,5,0))
names(fills) <- all_blocks$KEGGToGene
ActualGraph <- matrix(c(
  1:6,7,6,8:10,8,12:15,
  2:7,6,9,9:11,1,13:16
), ncol=2)
title <- "UDP-GlcNAc,\nN-glycans"
BuildCaseStudyForInputNetwork(AllPairs, subpreds, allpredictions,
                              all_blocks, fills,
                              ActualGraph, title, FALSE, VLC=0.8,
                              wrong_lty=3, MARGIN_OFFSET=M_offset,
                              GRAYLINE_VALS=c(-1.75,22), VERTNAMESCENTER=TRUE)

cat("\tPlotting legend...\n")
cols_legend <- c('#D81B60', '#E0A608', '#2B6DA8', '#45A649', 'grey40', 'black')
## legend for component algorithms
LEGEND_X <- -9.875
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
       col=cols_legend[6], lty=3, lwd=1.5, bty='n', xpd=NA, seg.len=2)
legend(x=LEGEND_X+4.000, y=LEGEND_Y-0.35, horiz=TRUE, cex=LEGEND_CEX, text.font=2,
       x.intersp=0.5,
       legend="Edge Matches KEGG",
       col=cols_legend[6], lty=1, lwd=1.5, bty='n', xpd=NA, seg.len=2)
legend(x=LEGEND_X+7.475, y=LEGEND_Y-0.35, horiz=TRUE, cex=LEGEND_CEX, text.font=2,
       x.intersp=0.5,
       legend="Direction in KEGG",
       col=cols_legend[6], lty=1, lwd=1.5, bty='n', xpd=NA, seg.len=2)
polygon(x=rep(LEGEND_X+8.025,3)+c(0,0,0.075),
        y=rep(LEGEND_Y-0.625,3)+c(0.04,-0.04,0),
        col='black', xpd=NA)

layout(1)
dev.off()
