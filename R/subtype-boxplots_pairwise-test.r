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
pairwiseTest <- pairwiseTestDatasets <- as.list(NULL)
subtypesOrder <- c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")

for(k in 1:length(geneSymbol)){
  patientsSubType <- pairwiseTestDatasets <- as.list(NULL)
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
    pairwiseTestDatasets[[i]] <- pairwise.wilcox.test(x=exprsDataAll, g=as.factor(patientDataAll))                     
    count <- count + 1
  }
  pairwiseTest[[k]] <-pairwiseTestDatasets
}
names(pairwiseTest) <- geneSymbol




