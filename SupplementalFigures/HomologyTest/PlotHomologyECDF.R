load(file.path(datadir, "AllPIDs.RData"))
CPp <- PIDs$Complex$Positive
CPn <- PIDs$Complex$Negative
MPp <- PIDs$Module$Positive
MPn <- PIDs$Module$Negative
MFP3 <- PIDs$ModuleFP$Positive[,1]
MFP4 <- PIDs$ModuleFP$Positive[,2]
MFP5 <- PIDs$ModuleFP$Positive[,3]

legend_cex <- 0.8207
P_LWD <- 1.25
BITSCORE_XLIM <- 120
L_INSET <- c(0,0.01)
cols <- c('#45A649','#D81B60', '#1E88E5', '#E0A608', "#824484")

pdf(file=file.path(figdir, "SXX_ECDFHomology_v2.pdf"), width = 4.3*2, height=4.3, onefile=TRUE)
layout(matrix(1:6, nrow=2, byrow=TRUE))
#pdf(file=file.path(figdir, "SXX_ECDFHomology.pdf"), width = 4.3*2, height=4.3, onefile=TRUE)
#layout(matrix(1:2, nrow=1, byrow=TRUE))
par(mar=c(2.5,2.5,1.2,0.5)+0.1, mgp=c(1.5,0.5,0))

plot.ecdf(CPp, col=cols[1], pch='.', verticals=TRUE,
          xlab='Percent Identity (PID)', ylab="Cumulative Probability", main="KEGG Complexes",
          xlim=c(0,1)
          )
plot.ecdf(CPn, col=cols[2], add=TRUE, pch='.', verticals = TRUE)
for(i in c(0.2,0.4)){
  abline(v=i, col=ifelse(i==0.2, cols[3], cols[4]), lwd=P_LWD)
  p1 <- sum(CPp < i) / length(CPp)
  p2 <- sum(CPn < i) / length(CPn)
  points(x=c(i,i), y=c(p1,p2), col=cols[1:2], lwd=P_LWD)
  text(x=rep(i+0.075,2), y=c(p1,p2)-c(0.015, 0.025), font=2,
       labels=sprintf("%.3f", c(p1,p2)), col=cols[1:2], cex=legend_cex)
}
legend("bottomright", legend=c("Within KO Group", "Between KO Groups"),
       col=cols[1:2], inset=L_INSET, lty=1, cex=legend_cex, bty='n')

plot.ecdf(MPp, col=cols[1], pch='.', verticals=TRUE,
          xlab='Percent Identity (PID)', ylab="", main="KEGG Modules",
          xlim=c(0,1)
          )
plot.ecdf(MPn, col=cols[2], add=TRUE, pch='.', verticals = TRUE)
for(i in c(0.2,0.4)){
  abline(v=i, col=ifelse(i==0.2, cols[3], cols[4]), lwd=P_LWD)
  p1 <- sum(MPp < i) / length(MPp)
  p2 <- sum(MPn < i) / length(MPn)
  points(x=c(i,i), y=c(p1,p2), col=cols[1:2], lwd=P_LWD)
  text(x=rep(i+0.075,2), y=c(p1,p2)-c(0.015, 0.025), font=2,
       labels=sprintf("%.3f", c(p1,p2)), col=cols[1:2], cex=legend_cex)
}
legend("bottomright", legend=c("Within Module Block", "Between Module Blocks"),
       col=cols[1:2], inset=L_INSET, lty=1, cex=legend_cex, bty='n')

## "direct connection" misclasses from multiclass benchmark
plot.ecdf(MFP3, col=cols[1], pch='.', verticals=TRUE, xlim=c(0,1),
          xlab='Percent Identity (PID)', ylab="", main="'Direct Connection' Misclasses",
)
plot.ecdf(MFP4, col=cols[2], add=TRUE, pch='.', verticals = TRUE)
plot.ecdf(MFP5, col=cols[5], add=TRUE, pch='.', verticals = TRUE)
for(i in c(0.2,0.4)){
  abline(v=i, col=ifelse(i==0.2, cols[3], cols[4]), lwd=P_LWD)
  p1 <- sum(MFP3 < i) / length(MFP3)
  p2 <- sum(MFP4 < i) / length(MFP4)
  p3 <- sum(MFP5 < i) / length(MFP5)
  points(x=c(i,i,i), y=c(p1,p2,p3), col=cols[c(1:2,5)], lwd=P_LWD)
  # text(x=rep(i+0.065,3), y=c(p1,p2,p3)-c(0.015, 0.025, 0.015), font=2,
  #      labels=sprintf("%.3f", c(p1,p2,p3)), col=cols[c(1:2,5)], cex=legend_cex)
}
legend("bottomright", legend=c("Same Pathway", "Same Global Pathway", "No Connection"),
       col=cols[c(1,2,5)], inset=L_INSET, lty=1, cex=legend_cex, bty='n')

### Plotting the results from BLAST
xlab <- 'BLASTP Best Bitscore'
CPp <- PIDs$ComplexBlast$Positive
CPn <- PIDs$ComplexBlast$Negative
MPp <- PIDs$ModuleBlast$Positive
MPn <- PIDs$ModuleBlast$Negative
MFP3 <- PIDs$ModuleFPBlast$Positive[,1]
MFP4 <- PIDs$ModuleFPBlast$Positive[,2]
MFP5 <- PIDs$ModuleFPBlast$Positive[,3]

#xmax <- quantile(CPp, 0.995)
xmax <- BITSCORE_XLIM
plot.ecdf(CPp, col=cols[1], pch='.', verticals=TRUE,
          xlab=xlab, ylab="Cumulative Probability", main="",
          xlim=c(0, xmax), ylim=c(0,1),
)
plot.ecdf(CPn, col=cols[2], add=TRUE, pch='.', verticals = TRUE)
legend("right", legend=c("Within KO Group", "Between KO Groups"),
       col=cols[1:2], inset=L_INSET, lty=1, cex=legend_cex, bty='n')

#xmax <- quantile(MPp, 0.995)
xmax <- BITSCORE_XLIM
plot.ecdf(MPp, col=cols[1], pch='.', verticals=TRUE, main='',
          xlab=xlab, ylab="",
          xlim=c(0,xmax), ylim=c(0,1)
)
plot.ecdf(MPn, col=cols[2], add=TRUE, pch='.', verticals = TRUE)
legend("right", legend=c("Within Module Block", "Between Module Blocks"),
       col=cols[1:2], inset=L_INSET, lty=1, cex=legend_cex, bty='n')

## "direct connection" misclasses from multiclass benchmark
#xmax <- quantile(MFP4, 0.995)
xmax <- BITSCORE_XLIM
plot.ecdf(MFP3, col=cols[1], pch='.', verticals=TRUE, main='',
          xlab=xlab, ylab="",
          xlim=c(0,xmax), ylim=c(0,1)
)
plot.ecdf(MFP4, col=cols[2], add=TRUE, pch='.', verticals = TRUE)
plot.ecdf(MFP5, col=cols[5], add=TRUE, pch='.', verticals = TRUE)
## plot vertical line for MFP3
bsc <- 60
p3 <- sum(MFP3 <= bsc) / length(MFP3)
abline(v=bsc, col=cols[4], lwd=P_LWD)
points(x=bsc, y=p3, col=cols[1], lwd=P_LWD)
text(x=bsc+7.5, y=p3-0.05, font=2, labels=sprintf("%.03f", p3), col=cols[1], cex=legend_cex)
legend("bottomright", legend=c("Same Pathway", "Same Global Pathway", "No Connection"),
       col=cols[c(1,2,5)], inset=L_INSET, lty=1, cex=legend_cex, bty='n')




# plot.ecdf(COp, col=cols[1], pch='.', verticals=TRUE,
#           xlab='Percent Identity (PID)', ylab="Cumulative Probability",
#           main="CORUM",
#           mgp=c(2,0.5,0))
# plot.ecdf(COn, col=cols[2], add=TRUE, pch='.', verticals = TRUE)
# for(i in c(0.2,0.4)){
#   abline(v=i, col=ifelse(i==0.2, cols[3], cols[4]), lwd=P_LWD)
#   p1 <- sum(COp < i) / length(COp)
#   p2 <- sum(COn < i) / length(COn)
#   points(x=c(i,i), y=c(p1,p2), col=cols[1:2], lwd=P_LWD)
#   text(x=rep(i+0.065,2), y=c(p1,p2)-c(0.015, 0.025), font=2,
#        labels=sprintf("%.3f", c(p1,p2)), col=cols[1:2], cex=legend_cex)
# }
# legend("bottomright", legend=c("Within CORUM Orthogroup", "Between CORUM Orthogroup"),
#        col=cols[1:2], inset=c(0.025,0.075), lty=1, cex=legend_cex)
dev.off()

#sum(COn>=0.4)
sum(MPn>=0.4)
sum(CPn>=0.4)
