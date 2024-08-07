load(file.path(datadir, "SupplementalData", "CORUM", "CORUMNuclearPredictions.RData"))

## Ensemble trained on CORUM
set.seed(454L)
EnsemblePreds <- matrix(nrow=nrow(FinalDataset), ncol=3)
colnames(EnsemblePreds) <- c("Logit", "RandomForest", "NeuralNetwork")
ppos <- sample(rep(1:5, length.out=sum(FinalDataset$isTP)))
pneg <- sample(rep(1:5, length.out=sum(!FinalDataset$isTP)))
allpos <- c(ppos,pneg)
res <- RawResults
res[is.na(res)] <- 0
trueanswer <- FinalDataset$isTP
rfmodels <- list()
for(j in seq_len(5)){
  sampv <- c(which(allpos!=j))
  testp <- c(which(allpos==j))
  traindf <- as.data.frame(res[sampv,])
  traindf$isTP <- as.integer(trueanswer[sampv])

  ## logistic regression
  model <- suppressWarnings(glm(isTP~., family='binomial', data=traindf))
  testdf <- as.data.frame(res[testp,])
  EnsemblePreds[testp,1] <- predict(model, testdf)

  ## random forest
  rfdf <- as.data.frame(res[sampv,])
  rfdf$isTP <- as.factor(trueanswer[sampv])
  x <- randomForest(isTP ~ ., data=rfdf, maxnodes=25, importance=TRUE)
  rfmodels[[j]] <- x
  EnsemblePreds[testp,2] <- predict(x, testdf, type='prob')[,"TRUE"]

  ## neuralnetwork
  x <- neuralnet(isTP ~ ., hidden=12, threshold=0.5, stepmax=1e6,
                 data=traindf, linear.output=FALSE)
  EnsemblePreds[testp,3] <- predict(x, testdf)[,1]
}
ens_rocs <- lapply(seq_len(ncol(EnsemblePreds)), \(x) vcheckans(EnsemblePreds[,x], trueanswer))
names(ens_rocs) <- colnames(EnsemblePreds)

## Transfer learning performance
load(file.path(datadir, "Modules", "EnsembleModels.RData"))
#colnames(res)[8] <- "GeneVector"
#colnames(res)[9] <- "SequenceInfo"
res <- as.data.frame(res)
TransferRes <- lapply(1:3, \(x) matrix(NA_real_, nrow=nrow(FinalDataset), ncol=5))
for(i in seq_len(5)){
  TransferRes[[1]][,i] <- predict(LRModels[[i]], res)
  TransferRes[[2]][,i] <- predict(RFModels[[i]], res, type='prob')[,'TRUE']
  TransferRes[[3]][,i] <- 1-predict(NNModels[[i]], res)[,1]
}

TransferRes <- lapply(TransferRes, \(x) vcheckans(rowMeans(x), trueanswer))
names(TransferRes) <- c("Logit", "RandomForest", "NeuralNetwork")

print(vapply(ens_rocs, \(x) x$AUROC, numeric(1)))
print(vapply(TransferRes, \(x) x$AUROC, numeric(1)))

## Save results
save(RawResults, FinalDataset, roc_list, ens_rocs, EnsemblePreds, rfmodels, TransferRes,
     file=file.path(datadir, "SupplementalData", "CORUM", 'CORUMNuclearPredictions.RData'))
