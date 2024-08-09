txtpath <- 'COG.mappings.v12.0.txt'

allcogmap <- read.table(txtpath, sep='\t')
pos_kegg <- which(grepl('KEGG', allcogmap[,5L]))

keggpos <- allcogmap[pos_kegg,]
kegg_annots <- keggpos[,5]
kegg_annots <- vapply(keggpos[,5], \(x){
  x <- strsplit(x, ' ')[[1]]
  p <- which(grepl('^[a-z]{3,4}:[^;]+$', x))
  if(length(p) > 0) return(x[p[1]])
  return('')
}, character(1L))
names(kegg_annots) <- NULL

keggpos2 <- keggpos[which(kegg_annots!=''),]
keggpos2 <- cbind(keggpos2[,1:4], kegg_annots[kegg_annots!=''])
colnames(keggpos2) <- c('protein', 'start_pos', 'end_pos', 'orthogroup', 'kegg_name')

StringToKEGG <- keggpos2
save(StringToKEGG, file='StringToKEGG.RData')
