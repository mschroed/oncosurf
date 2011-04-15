###################################################
## load libraries
###################################################
library(survcomp)
library(Biobase)
library(genefu)
library(breastCancerNKI)
library(breastCancerTRANSBIG)
library(breastCancerUPP)
library(breastCancerUNT)
library(breastCancerVDX)
library(breastCancerMAINZ)

###################################################
### import data
###################################################
data(breastCancerData)
data(mainz,transbig,upp,unt,vdx,nki)
data(scmgene)

###################################################
### create variables
###################################################
gsList <- fData(mainz7g)[,"Gene.symbol"]
gidList <- fData(mainz7g)[,"Gene.ID"]
probesNKI <- as.character(fData(nki7g)[,"probe"])
probesAffy <- fData(mainz7g)[,"probe"]

cinfo <- colnames(pData(mainz7g))
data.all7g <- c("mainz7g"=mainz7g, "transbig7g"=transbig7g, "upp7g"=upp7g, "unt7g"=unt7g, "vdx7g"=vdx7g, "nki7g"=nki7g)
data.all <- c("mainz"=mainz, "transbig"=transbig, "upp"=upp, "unt"=unt, "vdx"=vdx, "nki"=nki)

###################################################
### get subtypes
###################################################
count <- 1
subtypeData <- as.list(NULL)
for(i in data.all){
  subtypeData[[count]] <- subtype.cluster.predict(sbt.model=scmgene,
    data=t(exprs(i)), annot=fData(i), do.mapping=TRUE, do.scale=TRUE,
    do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
  count <- count + 1
}
names(subtypeData) <- names(data.all)         


###################################################
### draw boxplots for each gene in each dataset
###################################################

geneSymbol <- fData(mainz7g)[,"Gene.symbol"]
probeAffy <- fData(mainz7g)[,"probe"]
probeAgi <- as.character(fData(nki7g)[,"probe"])
subtypesOrder <- c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")

pdf("subtype-boxplots-7genes.pdf")

for(k in 1:length(geneSymbol)){
  par(mfrow=c(2,3))
  patientsSubType <- as.list(NULL)
  count <- 1
  for(i in names(data.all)){
    dd <- gg <- as.list(NULL)
    ddMin <- ddMax <- NULL
    patientsSubType[[count]] <- as.character(subtypeData[[i]]$subtype2)
    for(j in 1:length(subtypesOrder)){
      gg[[j]] <- patientsSubType[[count]] == subtypesOrder[j]
      if(i != "nki"){
        dd[[j]] <- exprs(data.all[[i]])[probeAffy[k],gg[[j]]]
      } else {
        dd[[j]] <- exprs(data.all[[i]])[probeAgi[k],gg[[j]]]
      }
      if(is.null(ddMin)){
        ddMin <- min(dd[[j]], na.rm=TRUE)
        ddMax <- max(dd[[j]], na.rm=TRUE)
      }
      if(ddMin > min(dd[[j]], na.rm=TRUE)){
        ddMin <- min(dd[[j]], na.rm=TRUE)
      }
      if(ddMax < max(dd[[j]], na.rm=TRUE)){
        ddMax <- max(dd[[j]], na.rm=TRUE)
      }
    }
    names(dd) <- NULL
    
    boxplotplus2(x=dd, .las=3, .jit=0.75, .ylim=c(ddMin, ddMax), pt.cex=0.75,
      pt.col=c(rep("darkred", length(dd[[1]])), rep("darkgreen", length(dd[[2]])),
      rep("darkblue", length(dd[[3]])), rep("orange", length(dd[[4]]))), pt.pch=c(0, 9, 17),
      main=paste(geneSymbol[k], " in ", i), xlab="", yaxt="n", xaxt="n")
    labels <- subtypesOrder
    text(1:4, par("usr")[3] - 0.25, srt = 45, adj = 1,
      labels = labels, xpd = TRUE)                        
    count <- count + 1
  }
}
dev.off()




###################################################
### pairwise wilcox test
###################################################

geneSymbol <- fData(mainz7g)[,"Gene.symbol"]
probeAffy <- fData(mainz7g)[,"probe"]
probeAgi <- as.character(fData(nki7g)[,"probe"])
pairwiseTestLess <- pairwiseTestGreater <- pairwiseTestDatasetsLess <- pairwiseTestDatasetsGreater <- as.list(NULL)
subtypesOrder <- c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")

for(k in 1:length(geneSymbol)){
  patientsSubType <- pairwiseTestDatasetsLess <- pairwiseTestDatasetsGreater <- as.list(NULL)
  count <- 1
  for(i in names(data.all)){
    dd <- gg <- as.list(NULL)
    patientsSubType[[count]] <- as.character(subtypeData[[i]]$subtype2)
    for(j in 1:length(subtypesOrder)){
      gg[[j]] <- patientsSubType[[count]] == subtypesOrder[j]
      if(i != "nki"){
        dd[[j]] <- as.numeric(exprs(data.all[[i]])[probeAffy[k],gg[[j]]])
      } else {
        dd[[j]] <- as.numeric(exprs(data.all[[i]])[probeAgi[k],gg[[j]]])
      }
    }
    names(dd) <- subtypesOrder
    
    ## create factor with patients and their subtypes for dataset i
##    subtypesAllFactor <- as.list(NULL)
##    for(m in 1:length(subtypesOrder)){
##      subtypesAllFactor[[m]] <- rep(names(dd)[m], length(dd[[m]]))
##    }    
    
    #######################                 
    #### other approach (two vectors)
    #### #working
    exprsDataAll <- patientDataAll <- NULL
    if(i != "nki"){
      exprsDataAll <- as.numeric(exprs(data.all[[i]])[probeAffy[k],])
    } else {
      exprsDataAll <- as.numeric(exprs(data.all[[i]])[probeAgi[k],])
    }
    
    for(m in 1:length(subtypesOrder)){
      patientDataAll[gg[[m]]] <- subtypesOrder[m]
    }       
    levels(patientDataAll) <- subtypesOrder
    pairwiseTestDatasetsLess[[i]] <- pairwise.wilcox.test(x=exprsDataAll, g=as.factor(patientDataAll), alternative="less", p.adjust="none")
    pairwiseTestDatasetsGreater[[i]] <- pairwise.wilcox.test(x=exprsDataAll, g=as.factor(patientDataAll), alternative="greater", p.adjust="none")
    count <- count + 1
  }
  pairwiseTestLess[[k]] <-pairwiseTestDatasetsLess
  pairwiseTestGreater[[k]] <-pairwiseTestDatasetsGreater  
}
names(pairwiseTestLess) <- geneSymbol
names(pairwiseTestGreater) <- geneSymbol

## initiate results matrix
resultsMatrix <- rep(NA, length(subtypesOrder)^2)
dim(resultsMatrix) <- c(length(subtypesOrder),length(subtypesOrder))
colnames(resultsMatrix) <- subtypesOrder
rownames(resultsMatrix) <- subtypesOrder

combPvalueL <- rep(1, (length(subtypesOrder)-1)^2)
dim(combPvalueL) <- c(length(subtypesOrder)-1,length(subtypesOrder)-1)
colnames(combPvalueL) <- colnames(pairwiseTestLess[[1]][[1]]$p.value)
rownames(combPvalueL) <- rownames(pairwiseTestLess[[1]][[1]]$p.value)
combPvalueG <- combPvalueL

combinedPairwiseTestLess <- combinedPairwiseTestGreater <- as.list(NULL)
finalMatrix <- as.list(NULL)
for(i in geneSymbol){                              
  tmpPvaluesL <- tmpPvaluesG <- NULL
  for(j in names(data.all)){    
    tmpPvaluesL <- rbind(tmpPvaluesL, pairwiseTestLess[[i]][[j]]$p.value[lower.tri(pairwiseTestLess[[i]][[j]]$p.value,diag=TRUE)])    
    tmpPvaluesG <- rbind(tmpPvaluesG, pairwiseTestGreater[[i]][[j]]$p.value[lower.tri(pairwiseTestGreater[[i]][[j]]$p.value,diag=TRUE)])        
  }
  subtypePvaluesL <- subtypePvaluesG <- NULL
  for(m in 1:length(names(data.all))){
    subtypePvaluesL[m] <- combine.test(tmpPvaluesL[,m])
    subtypePvaluesG[m] <- combine.test(tmpPvaluesG[,m])
  }
  combPvalueL[lower.tri(combPvalue, diag=TRUE)] <- subtypePvaluesL
  combinedPairwiseTestLess[[i]] <- combPvalueL
  
  combPvalueG[lower.tri(combPvalue, diag=TRUE)] <- subtypePvaluesG
  combinedPairwiseTestGreater[[i]] <- combPvalueG
  
  orderPosL <- c(3,1,2,5,6,4)
  orderPosG <- c(3,1,5,2,6,4)
  resultsMatrix[lower.tri(resultsMatrix)] <- subtypePvaluesL[orderPosL]
  resultsMatrix[upper.tri(resultsMatrix)] <- subtypePvaluesG[orderPosG]
  finalMatrix[[i]] <- resultsMatrix
}



## xtable(finalMatrix[["ESR1"]], digits=3, display=c("s",rep("E",4)))

for(i in geneSymbols){
  print(xtable(finalMatrix[[i]], digits=3, display=c("s",rep("E",4))))
  flush.console()
}






