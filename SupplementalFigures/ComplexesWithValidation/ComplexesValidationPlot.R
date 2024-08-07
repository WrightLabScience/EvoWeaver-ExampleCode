outfile <- 'SXX_ComplexesVerification.pdf'
infiles <- c("ExtendedComplexEnsembleMethods.RData",
              "ComplexStatistics_complexgeneholdouts.RData",
              "ComplexStatistics_fullcomplexholdouts.RData")
titles <- c("Complexes Benchmark with Extended Ensemble Models",
            "Complexes Benchmark with Gene Holdouts",
            "Complexes Benchmark with Complex Holdouts")

alg_names <- c("RandomForest25",
               "RandomForest2", "RandomForest10",
               "RandomForest100", "RandomForestInf",
               "NeuralNetwork1", "NeuralNetwork2",
               "NeuralNetwork3", "NeuralNetwork4",
               "NeuralNetworkBig1", "NeuralNetworkBig2",
               "NeuralNetworkBig3", "NeuralNetworkBig4",
               "Logit")
               #"SVMlinear", "SVMradial", "SVMpolynomial", "SVMsigmoid")
alg_newnames <- c("Random Forest (25 nodes)",
                  "Random Forest (2 nodes)",
                  "Random Forest (10 nodes)",
                  "Random Forest (100 nodes)",
                  "Random Forest (Inf nodes)",
                  "Logistic Regression",
                  "Neural Network (1x12 nodes)",
                  "Neural Network (2x12 nodes)",
                  "Neural Network (3x12 nodes)",
                  "Neural Network (4x12 nodes)",
                  "Neural Network (1x24 nodes)",
                  "Neural Network (2x24 nodes)",
                  "Neural Network (3x24 nodes)",
                  "Neural Network (4x24 nodes)"
                  # "SVM (kernel=Linear)",
                  # "SVM (kernel=Radial)",
                  # "SVM (kernel=Polynomial)",
                  # "SVM (kernel=Sigmoid)"
                  )

colList <- list(red=c('#E82B70','#C80B50','#A80030','#780000'),
                blue=c('#1E88E5','#0E68C5','#0048A5','#0028C5'),
                yellow=c('#FFC107','#DFA100','#BF8100', "#734d00"),
                other=c('#5CDC63','#3CBC83', '#7CFC83', "#036e08"),
                black=c('black', 'gray40', 'gray70'))

col_select <- c(5,1:4,5,5,2:4,1:4)
col_sub <-    c(1,1,1,1,1,2,3,3,3,3,4,4,4,4)
text_bold <-  c(2,1,1,1,1,2,2,1,1,1,1,1,1,1)
ltys_used <- c(1,5,4,3)
## three plots, put the legend in the fourth area

all_aurocs <- vector('list', 3)

pdf(file.path(figdir,outfile), onefile=TRUE, width=4.3*2, height=4.3*2)
par(mar=c(3,3,1,0.5)+0.1, mgp=c(1.5,0.5,0))

layout(matrix(1:4, nrow=2, byrow=TRUE))
all_cols <- vapply(seq_along(col_select),
                   \(x) colList[[col_select[x]]][col_sub[x]],
                   character(1L))
for(i in seq_along(infiles)){
  load(file.path(datadir, infiles[i]))
  to_plot <- EnsembleComplexStatistics[alg_names]
  plot(0:1,0:1, xlim=c(0,1), ylim=c(0,1), type='l', lty=2,
       main=titles[i], cex.main=1,
       xlab="False positive rate", ylab="True positive rate",
       xaxs='i', yaxs='i')
  for(j in seq_along(to_plot)){
    lines(x=to_plot[[j]]$FPR, y=to_plot[[j]]$TPR,
          col=all_cols[j],
          lty=ltys_used[col_sub[j]])
  }
  all_aurocs[[i]] <- vapply(to_plot, \(x) x$AUROC, numeric(1L))
}

## Make a legend in last plot
X_OFF <- -0.15
Y_OFF <- 0
LEGEND_CEX <- 0.9
plot(NULL, xlim=0:1, ylim=0:1,
     axes=FALSE, frame.plot=FALSE, main='',
     ylab='',xlab='')

alg_newnames <- c(alg_newnames[grepl("Random", alg_newnames)], '',
                  alg_newnames[grepl("Logistic", alg_newnames)], '',
                  alg_newnames[grepl("Neural", alg_newnames)], '',
                  "Random Guessing"
                  )
all_aurocs <- lapply(all_aurocs, \(x) c(x, 0.5))
ltys <- rep(NA, length(alg_newnames))
fonts <- rep(1, length(alg_newnames))
fonts[alg_newnames != ''] <- c(text_bold, 1)
ltys[alg_newnames != ''] <- c(ltys_used[col_sub], 2)
cols <- rep('black', length(alg_newnames))
cols[alg_newnames != ''] <- c(all_cols, "black")
legend('topleft',
       legend=alg_newnames,
       col=cols,
       lty=ltys,
       lwd=1,
       bty='n',
       title="\nAlgorithm", title.font=2, title.adj=0.275,
       inset=c(X_OFF,Y_OFF), xpd=NA,
       cex=LEGEND_CEX, text.font=fonts
       )

all_aurocs <- lapply(all_aurocs, \(x) sprintf("%.03f", x))
MULT <- 0.2
titles <- c("Complexes\nBenchmark", "Gene\nHoldout", "Complex\nHoldout")
for(i in seq_along(all_aurocs)){
  alg_aurocs <- rep("", length(alg_newnames))
  alg_aurocs[alg_newnames!=''] <- all_aurocs[[i]]
  fonts <- rep(1, length(alg_aurocs))
  tmp <- as.numeric(alg_aurocs)
  fonts[which(tmp == max(tmp, na.rm=TRUE))] <- 2
  legend('topleft',
         legend=alg_aurocs,
         bty='n',
         text.font=fonts, adj=0.25,
         title=titles[i], title.font=2, title.adj=0.5,
         inset=c(X_OFF +0.37+i*MULT+ifelse(i==2,0.0225,0), Y_OFF), xpd=NA,
         cex=LEGEND_CEX
         )
}

## Add inset of low FPR

# left, right, bottom, top
lrbt <- c(0.29, 0.47, 0.59, 0.77)
lrbt_coords <- list(lrbt,
                    lrbt + c(0.5,0.5,0,0),
                    lrbt + c(0,0,-0.5,-0.5))
for(i in seq_along(infiles)){
  load(file.path(datadir, infiles[i]))
  to_plot <- EnsembleComplexStatistics[alg_names]
  par(fig=lrbt_coords[[i]], new=TRUE, mar=c(0,0,0,0))
  plot(c(0,1), c(0,1), type='l', xaxs='i', yaxs='i',
       col='black', lty=2, lwd=1, xlab='', ylab='',
       ylim=c(0, 1), xlim=c(0, 0.02),
       main='', axes=FALSE, frame=TRUE,
       oma=c(0,0,0,0))
  # x axis
  axis(1, seq(0,0.02,0.005),
       labels = c("0.00", '', '0.01', '', '0.02'),
       cex.axis=1, tck=-0.05, mgp=c(3,0.5,0), cex.lab=1)
  axis(2, seq(0,1.0,0.2),
       labels= c("0.0", '', '0.4', '', '0.8', ''),
       cex.axis=1, tck=-0.05, mgp=c(3,0.5,0), cex.lab=1)
  for(j in seq_along(to_plot)){
    lines(x=to_plot[[j]]$FPR, y=to_plot[[j]]$TPR,
          col=all_cols[j],
          lty=ltys_used[col_sub[j]])
  }
}

dev.off(dev.list()['pdf'])
