basedir <- getwd()
p <- file.path(basedir, 'Data', 'OtherData')
f <- file.path(p, "AllHumanGenes.fa")
x <- system(paste0("blastp -query ", f, " -subject ", f, " -outfmt 6"), intern=TRUE)

tab <- read.table(text=x, sep='\t')

colnames(tab) <- c('Query_seqid', 'Subj_seqid', 'local pident', 'length', 'mismatch',
                      'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
BlastResultHuman <- tab
save(BlastResultHuman, file=file.path(p, "HumanSelfBlastResult.RData"))
subtab <- tab[tab[,1] == tab[,2],]

SelfHits <- subtab
save(BlastResultHuman, SelfHits, file=file.path(p, "HumanSelfBlastResult.RData"))


AllResultFileName <- file.path(basedir, "Data", "SupplementalData", "CORUM", "ExternalAlgorithmResultsCORUM.RData")
load(AllResultFileName, v=T)
n <- names(AllResults$bitscores[,'hsa'])
bb <- SelfHits$bitscore
names(bb) <- SelfHits$Query_seqid
AllResults$bitscores[,'hsa'] <- bb[n]
save(AllResults, file=AllResultFileName)
