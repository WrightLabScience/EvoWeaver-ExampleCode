out_fname <- "SXX_Misclasses.pdf"

pdf(file.path(figdir, out_fname), onefile=TRUE, width=4.3*2, height=4.3)
par(mgp=c(2,0.5,0))
tab <- read.xlsx(file.path(datadir, "Multiclass", "Top100Misclasses.xlsx"), 1,
                 colIndex = c(7,33,34))

tab[,1] <- as.numeric(tab[,1])
tab[,2] <- as.numeric(tab[,2])
tab[is.infinite(tab[,2]),2] <- 20
pathways <- tab[,3]
all_paths <- table(pathways)
singletons <- names(all_paths)[all_paths==1]
pathways[pathways%in%singletons] <- "Other"
all_paths <- table(pathways)
cols <- palette.colors(n=length(all_paths), palette="Polychrome 36")#[c(seq(2,length(all_paths)), 1L)]
all_paths <- sort(all_paths, decreasing=TRUE)
op <- which(names(all_paths) == "Other")
all_paths <- all_paths[c(seq_len(op-1L), seq(op+1L, length(all_paths)), op)]
cols <- rev(cols)

plot(x=tab[,1], y=tab[,2], xlim=c(1,0.7), ylim=c(0,20),
     xlab='Ensemble confidence in "Direct Connection"',
     ylab="Degrees of separation in KEGG",
     #main=paste0("Top ", nrow(tab), " Misclasses"),
     pch=19, col=cols[match(pathways, names(all_paths))], yaxt='n')

axis(side = 2, at=seq(0,20,by=5), labels=c(as.character(seq(0,15,5)), "Inf"))

LINE_W <- 0.0025
X_BPOS <- 1.012 + c(-1,1)*LINE_W
Y_BPOS <- c(18.3,18.7)
lines(x=rep(mean(X_BPOS),2), y=c(mean(Y_BPOS)-0.075, mean(Y_BPOS)-0.425), col='white', lwd=3, xpd=NA)
lines(x=X_BPOS, y=Y_BPOS, xpd=NA, lwd=1)
lines(x=X_BPOS, y=Y_BPOS-0.5, xpd=NA, lwd=1)
## top connection is BioC-BioD
xp <- tab[1,1]
yp <- tab[1,2]
arrows(x1=xp*0.9981, y1=yp*1.015,
       x0=xp*0.975, y0=yp*1.225,
       length=0.05, lwd=2)
#text(x=xp*0.959, y=yp*0.865, font=2,
text(x=xp*0.959, y=yp*1.265, font=2,
     labels = "BioC-BioD", cex=0.75, adj=0.5)

legend("topleft", cex=0.75, inset=c(0.0125,0.045),
       title="Biosynthetic Pathway", title.font=2,
       pch=19, col=cols, legend=paste0(names(all_paths), ' (', all_paths, ')'))

dev.off(dev.list()['pdf'])
