load(file.path(datadir, 'ComplexPairwise_PID.RData'))
CPp <- PID_positive
CPn <- PID_negative

load(file.path(datadir, "ModulePairwise_PID.RData"))
MPp <- PID_positive
MPn <- PID_negative

load(file.path(datadir, 'CorumPairwise_PID.RData'))
COp <- PID_positive
COn <- PID_negative

legend_cex <- 0.75
P_LWD <- 1.25
cols <- c('#45A649','#D81B60', '#1E88E5', '#E0A608', "#824484")

#pdf(file=file.path(figdir, "SXX_ECDFHomology.pdf"), width = 4.3*2, height=4.3*2, onefile=TRUE)
#layout(matrix(1:4, nrow=2, byrow=TRUE))
pdf(file=file.path(figdir, "SXX_ECDFHomology.pdf"), width = 4.3*2, height=4.3, onefile=TRUE)
par(mar=c(3,2.5,2,0.5)+0.1, mgp=c(1.5,0.5,0))
layout(matrix(1:2, nrow=1, byrow=TRUE))

plot.ecdf(CPp, col=cols[1], pch='.', verticals=TRUE,
          xlab='Percent Identity (PID)', ylab="Cumulative Probability", main="KEGG Complexes",
          )
plot.ecdf(CPn, col=cols[2], add=TRUE, pch='.', verticals = TRUE)
for(i in c(0.2,0.4)){
  abline(v=i, col=ifelse(i==0.2, cols[3], cols[4]), lwd=P_LWD)
  p1 <- sum(CPp < i) / length(CPp)
  p2 <- sum(CPn < i) / length(CPn)
  points(x=c(i,i), y=c(p1,p2), col=cols[1:2], lwd=P_LWD)
  text(x=rep(i+0.065,2), y=c(p1,p2)-c(0.015, 0.025), font=2,
       labels=sprintf("%.3f", c(p1,p2)), col=cols[1:2], cex=legend_cex)
}
legend("bottomright", legend=c("Within KO Group", "Between KO Groups"),
       col=cols[1:2], inset=c(0.025,0.075), lty=1, cex=legend_cex)

plot.ecdf(MPp, col=cols[1], pch='.', verticals=TRUE,
          xlab='Percent Identity (PID)', ylab="", main="KEGG Modules",
          )
plot.ecdf(MPn, col=cols[2], add=TRUE, pch='.', verticals = TRUE)
for(i in c(0.2,0.4)){
  abline(v=i, col=ifelse(i==0.2, cols[3], cols[4]), lwd=P_LWD)
  p1 <- sum(MPp < i) / length(MPp)
  p2 <- sum(MPn < i) / length(MPn)
  points(x=c(i,i), y=c(p1,p2), col=cols[1:2], lwd=P_LWD)
  text(x=rep(i+0.065,2), y=c(p1,p2)-c(0.015, 0.025), font=2,
       labels=sprintf("%.3f", c(p1,p2)), col=cols[1:2], cex=legend_cex)
}
legend("bottomright", legend=c("Within Module Block", "Between Module Blocks"),
       col=cols[1:2], inset=c(0.025,0.075), lty=1, cex=legend_cex)

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

sum(COn>=0.4)
sum(MPn>=0.4)
sum(CPn>=0.4)
