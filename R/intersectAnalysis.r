#######################################
### VENN DIAGRAM
#######################################
library(survcomp)
library(Biobase)

library(breastCancerNKI)
library(breastCancerTRANSBIG)
library(breastCancerUPP)
library(breastCancerUNT)
library(breastCancerVDX)
library(breastCancerMAINZ)

###################################################
## load data
###################################################
data(mainz,transbig,upp,unt,vdx, nki)

##load cindices first
##load("Rdata/cindices.Rdata")

## define number of probes from each dataset (top x and top x of highest/lowest cindices)
nprobes <- 10000

getAnnotgid <- function(x) {
  rr <- NULL
  for(i in x){
    tt <- fData(mainz)[ ,colnames(fData(mainz))=='EntrezGene.ID']==i
    tt[is.na(tt)]<-F
##  fData(mainz)[tt,c("probe", "EntrezGene.ID", "Gene.symbol", "Gene.title", "GO.Function", "GO.Process", "GO.Component")]
    rr <- rbind(rr,fData(mainz)[tt, ])
  }
  return(rr)
}

getAnnot <- function(x, pName) {
  if(pName == "mainz") {
    rr <- NULL
    paste("In mainz")
    for(i in x){
      paste("in loop")
      rr <- rbind(rr,fData(mainz)[fData(mainz)[i,"probe"],])
    }
    return(rr)
  } else if(pName == "transbig"){
    rr <- NULL
    for(i in x){
      rr <- rbind(rr,fData(transbig)[fData(transbig)[i,"probe"],])
    }
    return(rr)
  } else if(pName == "upp"){
    rr <- NULL
    for(i in x){
      rr <- rbind(rr,fData(upp)[fData(upp)[i,"probe"],])
    }
    return(rr)
  } else if(pName == "unt"){
    rr <- NULL
    for(i in x){
      rr <- rbind(rr,fData(unt)[fData(unt)[i,"probe"],])
    }
    return(rr)
  } else if(pName == "vdx"){
    rr <- NULL
    for(i in x){
      rr <- rbind(rr,fData(vdx)[fData(vdx)[i,"probe"],])
    }
    return(rr)
  } else if(pName == "nki"){
    rr <- NULL
    for(i in x){
      rr <- rbind(rr,fData(nki)[fData(nki)[i,"probe"],])
    }
    return(rr)
  }      
}

#######################
## top low/high risk genes MAINZ
#######################
##qm <- quantile(cindexall.mainz[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE)
##cindexall.mainz[cindexall.mainz[,"cindex"] > qm[[1]],]
##low.risk.probes.mainz <- as.data.frame(cindexall.mainz[cindexall.mainz[,"cindex"] > qm[[2]],])
##high.risk.probes.mainz <- as.data.frame(cindexall.mainz[cindexall.mainz[,"cindex"] < qm[[1]],])
low.risk.probes.mainz.full <- as.data.frame(cindexall.mainz[order(cindexall.mainz[,"cindex"], decreasing=TRUE),])
low.risk.probes.mainz <- low.risk.probes.mainz.full[1:nprobes,]
rNames <- rownames(low.risk.probes.mainz)
rowAnnot <- getAnnot(rNames,"mainz")
low.risk.probes.mainz <- cbind(low.risk.probes.mainz,rowAnnot)

high.risk.probes.mainz.full <- as.data.frame(cindexall.mainz[order(cindexall.mainz[,"cindex"]),])
high.risk.probes.mainz <- high.risk.probes.mainz.full[1:nprobes,]
rNames <- rownames(high.risk.probes.mainz)
rowAnnot <- getAnnot(rNames,"mainz")
high.risk.probes.mainz <- cbind(high.risk.probes.mainz,rowAnnot)



#######################
## top low/high risk genes TRANSBIG
#######################
low.risk.probes.transbig.full <- as.data.frame(cindexall.transbig[order(cindexall.transbig[,"cindex"], decreasing=TRUE),])
low.risk.probes.transbig <- low.risk.probes.transbig.full[1:nprobes,]
rNames <- rownames(low.risk.probes.transbig)
rowAnnot <- getAnnot(rNames,"transbig")
low.risk.probes.transbig <- cbind(low.risk.probes.transbig,rowAnnot)

high.risk.probes.transbig.full <- as.data.frame(cindexall.transbig[order(cindexall.transbig[,"cindex"]),])
high.risk.probes.transbig <- high.risk.probes.transbig.full[1:nprobes,]
rNames <- rownames(high.risk.probes.transbig)
rowAnnot <- getAnnot(rNames,"transbig")
high.risk.probes.transbig <- cbind(high.risk.probes.transbig,rowAnnot)



#######################
## top low/high risk genes UPP
#######################
low.risk.probes.upp.full <- high.risk.probes.upp.full <- low.risk.probes.upp <- high.risk.probes.upp <- rr <- NULL

low.risk.probes.upp.full <- as.data.frame(cindexall.upp[order(cindexall.upp[,"cindex"], decreasing=TRUE),])
low.risk.probes.upp <- low.risk.probes.upp.full[1:nprobes,]
rNames <- rownames(low.risk.probes.upp)
rowAnnot <- getAnnot(rNames,"upp")
low.risk.probes.upp <- cbind(low.risk.probes.upp,rowAnnot)

high.risk.probes.upp.full <- as.data.frame(cindexall.upp[order(cindexall.upp[,"cindex"]),])
high.risk.probes.upp <- high.risk.probes.upp.full[1:nprobes,]
rNames <- rownames(high.risk.probes.upp)
rowAnnot <- getAnnot(rNames,"upp")
high.risk.probes.upp <- cbind(high.risk.probes.upp,rowAnnot)

#######################
## top low/high risk genes UNT
#######################
low.risk.probes.unt.full <- high.risk.probes.unt.full <- low.risk.probes.unt <- high.risk.probes.unt <- rr <- NULL

low.risk.probes.unt.full <- as.data.frame(cindexall.unt[order(cindexall.unt[,"cindex"], decreasing=TRUE),])
low.risk.probes.unt <- low.risk.probes.unt.full[1:nprobes,]
rNames <- rownames(low.risk.probes.unt)
rowAnnot <- getAnnot(rNames,"unt")
low.risk.probes.unt <- cbind(low.risk.probes.unt,rowAnnot)

high.risk.probes.unt.full <- as.data.frame(cindexall.unt[order(cindexall.unt[,"cindex"]),])
high.risk.probes.unt <- high.risk.probes.unt.full[1:nprobes,]
rNames <- rownames(high.risk.probes.unt)
rowAnnot <- getAnnot(rNames,"unt")
high.risk.probes.unt <- cbind(high.risk.probes.unt,rowAnnot)


#######################
## top low/high risk genes VDX
#######################
low.risk.probes.vdx.full <- high.risk.probes.vdx.full <- low.risk.probes.vdx <- high.risk.probes.vdx <- rr <- NULL

low.risk.probes.vdx.full <- as.data.frame(cindexall.vdx[order(cindexall.vdx[,"cindex"], decreasing=TRUE),])
low.risk.probes.vdx <- low.risk.probes.vdx.full[1:nprobes,]
rNames <- rownames(low.risk.probes.vdx)
rowAnnot <- getAnnot(rNames,"vdx")
low.risk.probes.vdx <- cbind(low.risk.probes.vdx,rowAnnot)

high.risk.probes.vdx.full <- as.data.frame(cindexall.vdx[order(cindexall.vdx[,"cindex"]),])
high.risk.probes.vdx <- high.risk.probes.vdx.full[1:nprobes,]
rNames <- rownames(high.risk.probes.vdx)
rowAnnot <- getAnnot(rNames,"vdx")
high.risk.probes.vdx <- cbind(high.risk.probes.vdx,rowAnnot)


#######################
## top low/high risk genes NKI
#######################
low.risk.probes.nki.full <- high.risk.probes.nki.full <- low.risk.probes.nki <- high.risk.probes.nki <- rr <- NULL

low.risk.probes.nki.full <- as.data.frame(cindexall.nki[order(cindexall.nki[,"cindex"], decreasing=TRUE),])
low.risk.probes.nki <- low.risk.probes.nki.full[1:nprobes,]
rNames <- rownames(low.risk.probes.nki)
rowAnnot <- getAnnot(rNames,"nki")
low.risk.probes.nki <- cbind(low.risk.probes.nki,rowAnnot)

high.risk.probes.nki.full <- as.data.frame(cindexall.nki[order(cindexall.nki[,"cindex"]),])
high.risk.probes.nki <- high.risk.probes.nki.full[1:nprobes,]
rNames <- rownames(high.risk.probes.nki)
rowAnnot <- getAnnot(rNames,"nki")
high.risk.probes.nki <- cbind(high.risk.probes.nki,rowAnnot)


## plot cindices distribution
plot(high.risk.probes.mainz.full[,"cindex"], ylim=c(0.2,0.8), xlim=c(0,25000), type="l")
lines(high.risk.probes.transbig.full[,"cindex"])
lines(high.risk.probes.upp.full[(1:length(high.risk.probes.upp.full)%%2)==1,"cindex"])
lines(high.risk.probes.unt.full[(1:length(high.risk.probes.unt.full)%%2)==1,"cindex"])
lines(high.risk.probes.vdx.full[,"cindex"])
lines(high.risk.probes.nki.full[,"cindex"])
abline(h=0.5)


#######################
## VENN DIAGRAM data
#######################
## compare with entrezgene id's
gLowMainz <- as.numeric(na.omit(low.risk.probes.mainz[,"EntrezGene.ID"]))
gLowTransbig <- as.numeric(na.omit(low.risk.probes.transbig[,"EntrezGene.ID"]))
gLowUpp <- as.numeric(na.omit(low.risk.probes.upp[,"EntrezGene.ID"]))
gLowUnt <- as.numeric(na.omit(low.risk.probes.unt[,"EntrezGene.ID"]))
gLowVdx <- as.numeric(na.omit(low.risk.probes.vdx[,"EntrezGene.ID"]))
gLowNki <- as.numeric(na.omit(low.risk.probes.nki[,"EntrezGene.ID"]))

gHighMainz <- as.numeric(na.omit(high.risk.probes.mainz[,"EntrezGene.ID"]))
gHighTransbig <- as.numeric(na.omit(high.risk.probes.transbig[,"EntrezGene.ID"]))
gHighUpp <- as.numeric(na.omit(high.risk.probes.upp[,"EntrezGene.ID"]))
gHighUnt <- as.numeric(na.omit(high.risk.probes.unt[,"EntrezGene.ID"]))
gHighVdx <- as.numeric(na.omit(high.risk.probes.vdx[,"EntrezGene.ID"]))
gHighNki <- as.numeric(na.omit(high.risk.probes.nki[,"EntrezGene.ID"]))


## or with probe id's
gLowMainz <- low.risk.probes.mainz[,"probe"]
gLowTransbig <- low.risk.probes.transbig[,"probe"]
gLowUpp <- low.risk.probes.upp[,"probe"]
gLowUnt <- low.risk.probes.unt[,"probe"]
gLowVdx <- low.risk.probes.vdx[,"probe"]
gLowNki <- as.character(low.risk.probes.nki[,"probe"])

gHighMainz <- high.risk.probes.mainz[,"probe"]
gHighTransbig <- high.risk.probes.transbig[,"probe"]
gHighUpp <- high.risk.probes.upp[,"probe"]
gHighUnt <- high.risk.probes.unt[,"probe"]
gHighVdx <- high.risk.probes.vdx[,"probe"]
gHighNki <- as.character(high.risk.probes.nki[,"probe"])



gLow <- cbind(gLowMainz,gLowTransbig,gLowUpp,gLowUnt,gLowVdx,gLowNki)
gHigh <- cbind(gHighMainz,gHighTransbig,gHighUpp,gHighUnt,gHighVdx,gHighNki)

colnames(gLow) <- c("MainzLow","TransbigLow", "UppLow", "UntLow", "VdxLow", "NkiLow")
colnames(gHigh) <- c("MainzHigh", "TransbigHigh", "UppHigh", "UntHigh", "VdxHigh", "NkiHigh")

write.table(cbind(gLow,gHigh), file="vennData-probes.csv", sep=" ")

## get the intersect of all datasets
geneIdAll <- intersect(gHighMainz,intersect(gHighTransbig,intersect(gHighUpp,intersect(gHighUnt,intersect(gHighVdx,gHighNki)))))
geneIdMinusVdx <- intersect(gHighMainz,intersect(gHighTransbig,intersect(gHighUpp,intersect(gHighUnt,gHighNki))))
geneIdMinusNki <- intersect(gHighMainz,intersect(gHighTransbig,intersect(gHighUpp,intersect(gHighUnt,gHighVdx))))
geneIdMinusMainz <- intersect(gHighNki,intersect(gHighTransbig,intersect(gHighUpp,intersect(gHighUnt,gHighVdx))))
geneIdMinusTransbig <- intersect(gHighMainz,intersect(gHighNki,intersect(gHighUpp,intersect(gHighUnt,gHighVdx))))
geneIdMinusUpp <- intersect(gHighMainz,intersect(gHighTransbig,intersect(gHighNki,intersect(gHighUnt,gHighVdx))))
geneIdMinusUnt <- intersect(gHighMainz,intersect(gHighTransbig,intersect(gHighUpp,intersect(gHighNki,gHighVdx))))


xx <- high.risk.probes.mainz[,"EntrezGene.ID"]==11394

iData



iData <- c(5961,1539,690)
xx <- rr <- tt <- NULL
for(i in iData) {
  xx <- (as.numeric(na.omit(high.risk.probes.mainz[,"EntrezGene.ID"])))==i
  rr <- low.risk.probes.mainz[xx,]
  tt <- rbind(tt,rr)
}  


intersect(low.risk.probes.mainz[1:500,"EntrezGene.ID"], intersect(low.risk.probes.transbig[1:500,"EntrezGene.ID"],intersect(low.risk.probes.upp[1:500,"EntrezGene.ID"],intersect(low.risk.probes.unt[1:500,"EntrezGene.ID"], intersect(low.risk.probes.vdx[1:500,"EntrezGene.ID"], low.risk.probes.nki[1:500,"EntrezGene.ID"])))))


intersect(low.risk.probes.mainz[1:500,"Gene.symbol"], intersect(low.risk.probes.transbig[1:500,"Gene.symbol"],intersect(low.risk.probes.upp[1:500,"Gene.symbol"],intersect(low.risk.probes.unt[1:500,"Gene.symbol"], intersect(low.risk.probes.vdx[1:500,"Gene.symbol"], low.risk.probes.nki[1:500,"NCBI.gene.symbol"])))))


intersect(high.risk.probes.mainz[1:500,"Gene.symbol"], intersect(high.risk.probes.transbig[1:500,"Gene.symbol"],intersect(high.risk.probes.upp[1:500,"Gene.symbol"],intersect(high.risk.probes.unt[1:500,"Gene.symbol"], intersect(high.risk.probes.vdx[1:500,"Gene.symbol"], high.risk.probes.nki[1:500,"NCBI.gene.symbol"])))))

#
########################
### ranking the genes in each dataset
########################
#rankHighMainz <- cbind(high.risk.probes.mainz[,"Gene.symbol"], rank(high.risk.probes.mainz[,"cindex"]))
#rankHighMainz <- as.data.frame(unique(rankHighMainz[complete.cases(rankHighMainz),]))
#uniqueRank <- unique(rankHighMainz[,1])
#rownames(rankHighMainz) <- rankHighMainz[,1]
#
#rankHighTransbig <- cbind(high.risk.probes.transbig[,"Gene.symbol"], rank(high.risk.probes.transbig[,"cindex"]))
#rankHighTransbig <- rankHighTransbig[complete.cases(rankHighTransbig),]
#
#rankHighUpp <- cbind(high.risk.probes.upp[,"Gene.symbol"], rank(high.risk.probes.upp[,"cindex"]))
#rankHighUpp <- rankHighUpp[complete.cases(rankHighUpp),]
#
#rankHighUnt <- cbind(high.risk.probes.unt[,"Gene.symbol"], rank(high.risk.probes.unt[,"cindex"]))
#rankHighUnt <- rankHighUnt[complete.cases(rankHighUnt),]
#
#rankHighVdx <- cbind(high.risk.probes.vdx[,"Gene.symbol"], rank(high.risk.probes.vdx[,"cindex"]))
#rankHighVdx <- rankHighVdx[complete.cases(rankHighVdx),]
#
#rankHighNki <- cbind(high.risk.probes.nki[,"HUGO.gene.symbol"], rank(high.risk.probes.nki[,"cindex"]))
#rankHighNki <- rankHighNki[complete.cases(rankHighNki),]
#
#geneNames <- intersect(rankHighMainz[,1],intersect(rankHighTransbig[,1],intersect(rankHighUpp[,1],rankHighUnt[,1])))
#
#head(rankHighMainz[order(as.numeric(rankHighMainz[,2])),])
#
#
#

#######################
## ranking the genes in each dataset
#######################
library(survcomp)
data(breastCancerData)

mainzfData <- fData(mainz7g)
cindex.mainz <- t(apply(X=exprs(mainz7g), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainz7g)[ ,"t.dmfs"], z=pData(mainz7g)[ ,"e.dmfs"]))
mainzfData <- cbind(mainzfData,cindexall.mainz,"rankHighRisk"=rank(cindex.mainz[,1]),"rankLowRisk"=rank(-cindex.mainz[,1]))

transbigfData <- fData(transbig7g)
cindex.transbig <- t(apply(X=exprs(transbig7g), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbig7g)[ ,"t.dmfs"], z=pData(transbig7g)[ ,"e.dmfs"]))
transbigfData <- cbind(transbigfData,cindex.transbig,"rankHighRisk"=rank(cindex.transbig[,1]),"rankLowRisk"=rank(-cindex.transbig[,1]))

uppfData <- fData(upp7g)
cindex.upp <- t(apply(X=exprs(upp7g), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(upp7g)[ ,"t.rfs"], z=pData(upp7g)[ ,"e.rfs"]))
uppfData <- cbind(uppfData, cindex.upp,"rankHighRisk"=rank(cindex.upp[,1]),"rankLowRisk"=rank(-cindex.upp[,1]))

untfData <- fData(unt7g)
cindex.unt <- t(apply(X=exprs(unt7g), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(unt7g)[ ,"t.dmfs"], z=pData(unt7g)[ ,"e.dmfs"]))
untfData <- cbind(untfData, cindex.unt,"rankHighRisk"=rank(cindex.unt[,1]),"rankLowRisk"=rank(-cindex.unt[,1]))

vdxfData <- fData(vdx7g)
cindex.vdx <- t(apply(X=exprs(vdx7g), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdx7g)[ ,"t.dmfs"], z=pData(vdx7g)[ ,"e.dmfs"]))
vdxfData <- cbind(vdxfData, cindex.vdx,"rankHighRisk"=rank(cindex.vdx[,1]),"rankLowRisk"=rank(-cindex.vdx[,1]))

nkifData <- fData(nki7g)
cindex.nki <- t(apply(X=exprs(nki7g), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nki7g)[ ,"t.dmfs"], z=pData(nki7g)[ ,"e.dmfs"]))
nkifData <- cbind(nkifData, cindex.nki,"rankHighRisk"=rank(cindex.nki[,1]),"rankLowRisk"=rank(-cindex.nki[,1]))

names <- as.character(intersect(mainzfData$Gene.symbol, intersect(transbigfData$Gene.symbol, intersect(uppfData$Gene.symbol, intersect(untfData$Gene.symbol, intersect(vdxfData$Gene.symbol, nkifData$HUGO.gene.symbol))))))



riskDF <- data.frame("genes"=names)
for(i in 1:length(names)){
    riskDF[i,2] <- mainzfData[mainzfData[,"Gene.symbol"]==names[i], "rankHighRisk"]
    riskDF[i,3] <- transbigfData[transbigfData[,"Gene.symbol"]==names[i], "rankHighRisk"]
    riskDF[i,4] <- uppfData[uppfData[,"Gene.symbol"]==names[i], "rankHighRisk"]
    riskDF[i,5] <- untfData[untfData[,"Gene.symbol"]==names[i], "rankHighRisk"]
    riskDF[i,6] <- vdxfData[vdxfData[,"Gene.symbol"]==names[i], "rankHighRisk"]
    riskDF[i,7] <- nkifData[nkifData[,"HUGO.gene.symbol"]==names[i], "rankHighRisk"]
    riskDF[i,8] <- mean(as.numeric(riskDF[i,2:7]))
}
colnames(riskDF) <- c("genes","MAINZ","TRANSBIG","UPP","UNT","VDX","NKI", "rankMean")

high.risk <- riskDF[order(riskDF$rankMean),]
low.risk <- riskDF[order(-riskDF$rankMean),]

                  
##################################################
## full datasets
##################################################

mainzfData <- fData(mainz)
transbigfData <- fData(transbig)
uppfData <- fData(upp)
untfData <- fData(unt)
vdxfData <- fData(vdx)
nkifData <- fData(nki)

mainzfData <- cbind(mainzfData, cindexall.mainz, "rankHighRisk"=rank(cindexall.mainz[,1]), "rankLowRisk"=rank(-cindexall.mainz[,1]))

transbigfData <- cbind(transbigfData, cindexall.transbig, "rankHighRisk"=rank(cindexall.transbig[,1]), "rankLowRisk"=rank(-cindexall.transbig[,1]))

uppfData <- cbind(uppfData, cindexall.upp, "rankHighRisk"=rank(cindexall.upp[,1]), "rankLowRisk"=rank(-cindexall.upp[,1]))

untfData <- cbind(untfData, cindexall.unt, "rankHighRisk"=rank(cindexall.unt[,1]), "rankLowRisk"=rank(-cindexall.unt[,1]))

vdxfData <- cbind(vdxfData, cindexall.vdx, "rankHighRisk"=rank(cindexall.vdx[,1]), "rankLowRisk"=rank(-cindexall.vdx[,1]))

nkifData <- cbind(nkifData, cindexall.nki, "rankHighRisk"=rank(cindexall.nki[,1]), "rankLowRisk"=rank(-cindexall.nki[,1]))


names <- na.omit(as.character(intersect(mainzfData$Gene.symbol, intersect(transbigfData$Gene.symbol, intersect(uppfData$Gene.symbol, intersect(untfData$Gene.symbol, intersect(vdxfData$Gene.symbol, nkifData$HUGO.gene.symbol)))))))


riskDF <- data.frame("genes"=names)
for(i in 1:length(names)){
    riskDF[i,2] <- mainzfData[na.omit(mainzfData[,"Gene.symbol"]==names[i]), "rankHighRisk"][1]
    riskDF[i,3] <- transbigfData[na.omit(transbigfData[,"Gene.symbol"]==names[i]), "rankHighRisk"][1]
    riskDF[i,4] <- uppfData[na.omit(uppfData[,"Gene.symbol"]==names[i]), "rankHighRisk"][1]
    riskDF[i,5] <- untfData[na.omit(untfData[,"Gene.symbol"]==names[i]), "rankHighRisk"][1]
    riskDF[i,6] <- vdxfData[na.omit(vdxfData[,"Gene.symbol"]==names[i]), "rankHighRisk"][1]
    riskDF[i,7] <- nkifData[na.omit(nkifData[,"HUGO.gene.symbol"]==names[i]), "rankHighRisk"][1]
    riskDF[i,8] <- mean(as.numeric(riskDF[i,2:7]))
    riskDF[i,9] <- sum(as.numeric(riskDF[i,2:7]))
}
colnames(riskDF) <- c("genes","MAINZ","TRANSBIG","UPP","UNT","VDX","NKI", "rankMean", "rankSum")

high.riskM <- riskDF[order(riskDF$rankMean),]
low.riskM <- riskDF[order(-riskDF$rankMean),]

high.riskS <- riskDF[order(riskDF$rankSum),]
low.riskS <- riskDF[order(-riskDF$rankSum),]




## union of genes
names <- na.omit(as.character(union(mainzfData$Gene.symbol, union(transbigfData$Gene.symbol, union(uppfData$Gene.symbol, union(untfData$Gene.symbol, union(vdxfData$Gene.symbol, nkifData$HUGO.gene.symbol)))))))

highRiskDF <- data.frame("genes"=names)
for(i in 1:length(names)){
    highRiskDF[i,2] <- min(mainzfData[na.omit(mainzfData[,"Gene.symbol"]==names[i]), "rankHighRisk"])
    highRiskDF[i,3] <- min(transbigfData[na.omit(transbigfData[,"Gene.symbol"]==names[i]), "rankHighRisk"])
    highRiskDF[i,4] <- min(uppfData[na.omit(uppfData[,"Gene.symbol"]==names[i]), "rankHighRisk"])
    highRiskDF[i,5] <- min(untfData[na.omit(untfData[,"Gene.symbol"]==names[i]), "rankHighRisk"])
    highRiskDF[i,6] <- min(vdxfData[na.omit(vdxfData[,"Gene.symbol"]==names[i]), "rankHighRisk"])
    highRiskDF[i,7] <- min(nkifData[na.omit(nkifData[,"HUGO.gene.symbol"]==names[i]), "rankHighRisk"])
    highRiskDF[i,8] <- mean(as.numeric(highRiskDF[i,2:7]), na.rm=TRUE)
    highRiskDF[i,9] <- sum(as.numeric(highRiskDF[i,2:7]), na.rm=TRUE)
    highRiskDF[i,10] <- sum(sort(as.numeric(highRiskDF[i,2:7]))[-6], na.rm=TRUE)
}
colnames(highRiskDF) <- c("genes","MAINZ","TRANSBIG","UPP","UNT","VDX","NKI", "rankMean", "rankSum", "rankSumMinusMax")

lowRiskDF <- data.frame("genes"=names)
for(i in 1:length(names)){
    lowRiskDF[i,2] <- min(mainzfData[na.omit(mainzfData[,"Gene.symbol"]==names[i]), "rankLowRisk"])
    lowRiskDF[i,3] <- min(transbigfData[na.omit(transbigfData[,"Gene.symbol"]==names[i]), "rankLowRisk"])
    lowRiskDF[i,4] <- min(uppfData[na.omit(uppfData[,"Gene.symbol"]==names[i]), "rankLowRisk"])
    lowRiskDF[i,5] <- min(untfData[na.omit(untfData[,"Gene.symbol"]==names[i]), "rankLowRisk"])
    lowRiskDF[i,6] <- min(vdxfData[na.omit(vdxfData[,"Gene.symbol"]==names[i]), "rankLowRisk"])
    lowRiskDF[i,7] <- min(nkifData[na.omit(nkifData[,"HUGO.gene.symbol"]==names[i]), "rankLowRisk"])
    lowRiskDF[i,8] <- mean(as.numeric(lowRiskDF[i,2:7]), na.rm=TRUE)
    lowRiskDF[i,9] <- sum(as.numeric(lowRiskDF[i,2:7]), na.rm=TRUE)
    lowRiskDF[i,10] <- sum(sort(as.numeric(lowRiskDF[i,2:7]))[-6], na.rm=TRUE)
}
colnames(lowRiskDF) <- c("genes","MAINZ","TRANSBIG","UPP","UNT","VDX","NKI", "rankMean", "rankSum", "rankSumMinusMax")


for(i in 1:length(names)){
#    lowRiskDF[i,10] <- sum(sort(as.numeric(lowRiskDF[i,2:7]))[-6], na.rm=TRUE)
#    highRiskDF[i,10] <- sum(sort(as.numeric(highRiskDF[i,2:7]))[-6], na.rm=TRUE)
}    



high.riskM <- highRiskDF[order(highRiskDF$rankMean),]
low.riskM <- lowRiskDF[order(lowRiskDF$rankMean),]

high.riskS <- highRiskDF[order(highRiskDF$rankSum),]
low.riskS <- lowRiskDF[order(lowRiskDF$rankSum),]

high.riskX <- highRiskDF[order(highRiskDF$rankSumMinusMax),]
low.riskX <- lowRiskDF[order(lowRiskDF$rankSumMinusMax),]

































#######################
## VENN DIAGRAM real data
#######################
gmainz <- fData(mainzSample)[,"EntrezGene.ID"]
gtransbig <- fData(transbigSample)[,"EntrezGene.ID"]
gupp <- fData(uppSample)[,"EntrezGene.ID"]
g <- cbind(gmainz,gtransbig,gupp)
c <- vennCounts(g)
vennDiagram(c)


c <- na.omit(c)

ab <- intersect(a,b)
ac <- intersect(a,c)
abc <- intersect(ab,ac)
abc



