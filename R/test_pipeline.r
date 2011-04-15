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
data(mainz,transbig,upp,unt,vdx,nki,breastCancerData)
data(scmgene)


###################################################
### save everything
###################################################
pdf("test-pipeline-run-results.pdf")


###################################################
### create variables
###################################################
subtypesOrder <- c("ER-/HER2-", "HER2+", "ER+/HER2- High Prolif", "ER+/HER2- Low Prolif")
data.all <- c("mainz"=mainz, "transbig"=transbig, "upp"=upp, "unt"=unt, "vdx"=vdx, "nki"=nki)
tc <- 10 * 365

## replace with your probes of interest
## must be the same order in all vectors
geneSymbols <- fData(mainz7g)[,"Gene.symbol"]
probesAffy <- fData(mainz7g)[,"probe"]
probesAgi <- as.character(fData(nki7g)[,"probe"])


###################################################
### get subtypes
###################################################
count <- 1
subtypeData <- as.list(NULL)
for(i in data.all){
  subtypeData[[count]] <- subtype.cluster.predict(sbt.model=scmgene,
    data=t(exprs(i)), annot=fData(i), do.mapping=TRUE, do.scale=TRUE,
    do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
  title(names(data.all)[count])
  count <- count + 1
}
names(subtypeData) <- names(data.all)     


###################################################
### draw boxplots for each gene in each dataset
###################################################
#pdf("subtype-boxplots.pdf")

for(k in 1:length(geneSymbols)){
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
        dd[[j]] <- exprs(data.all[[i]])[probesAffy[k],gg[[j]]]
      } else {
        dd[[j]] <- exprs(data.all[[i]])[probesAgi[k],gg[[j]]]
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
      main=paste(geneSymbols[k], " in ", i), xlab="", yaxt="n", xaxt="n")
    labels <- subtypesOrder
    text(1:4, par("usr")[3] - 0.25, srt = 45, adj = 1,
      labels = labels, xpd = TRUE)                        
    count <- count + 1
  }
}
#dev.off()



###################################################
### survival curves using subtypes, all datasets
###################################################

## create vectors for survival time and event
surv.time.all <- surv.event.all <- NULL
for(i in names(data.all)){
  if(i != "upp"){
    surv.time.all <- c(surv.time.all, pData(data.all[[i]])[ ,"t.dmfs"])
    surv.event.all <- c(surv.event.all, pData(data.all[[i]])[ ,"e.dmfs"])
  } else {
    surv.time.all <- c(surv.time.all, pData(data.all[[i]])[ ,"t.rfs"])
    surv.event.all <- c(surv.event.all, pData(data.all[[i]])[ ,"e.rfs"])
  }
}

## creating vector containing subtypes of all patients
subtypes.all <- NULL     
subtypes.all.length <- as.list(NULL)
for(i in subtypeData){
  subtypes.all <- c(subtypes.all, as.character(i$subtype2))
  for(j in subtypesOrder){
    tmpSubtypeList <- i$subtype2[as.character(i$subtype2) == j]
    tmpSubtypeList <- tmpSubtypeList[!is.na(tmpSubtypeList)]
    subtypes.all.length[[j]] <- c(subtypes.all.length[[j]], length(tmpSubtypeList))
  }
}
subtypes.all[is.na(subtypes.all)] <- F

## create KM survival curves for each gene
for(k in 1:length(geneSymbols)){    
  ## get gene expression for gene k
  geneExprs <- NULL
  for(m in names(data.all)){
    if(m != "nki"){
      geneExprs <- c(geneExprs, exprs(data.all[[m]])[probesAffy[k],])
    } else {
      geneExprs <- c(geneExprs, exprs(data.all[[m]])[probesAgi[k],])
    }
  }
  
  par(mfrow=c(2,2))
  ## go through all subtypes
  for(i in subtypesOrder){    
    stype <- timeTmp <- eventTmp <- exprsTmp <- NULL    
    stype <- subtypes.all==i                            
    timeTmp <- surv.time.all[stype]
    eventTmp <- surv.event.all[stype]    
    exprsTmp <- geneExprs[stype]                    
  
    pos <- 0
    mygroup <- NULL
    for(j in subtypes.all.length[[i]]){
      qq <- exprsTmp[(pos+1):(pos+j)]
      myq <- quantile(qq, probs=c(0.33, 0.66), na.rm=TRUE)
      qq[exprsTmp[(pos+1):(pos+j)] < myq[1]] <- 1
      qq[exprsTmp[(pos+1):(pos+j)] >= myq[1] & exprsTmp[(pos+1):(pos+j)] < myq[2]] <- 2
      qq[exprsTmp[(pos+1):(pos+j)] > myq[2]] <- 3
      qq <- factor(x=qq, levels=1:3)
      mygroup <- c(mygroup,qq)
      pos <- pos + j
    }                     
  
    surv.data <- censor.time(surv.time=timeTmp / 365, surv.event=eventTmp, time.cens=tc / 365)  
    dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
        
    km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title=c("Probability of survival for ",geneSymbols[k], i), sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))
  }
}

#dev.off()


###################################################
### compute cindices for different subtype groups
###################################################
for(subT in subtypesOrder){
  
###################################################
### compute concordance indices
###################################################
tmpST <- as.character(subtypeData[["mainz"]]$subtype2) == subT
cindexMainz <- t(apply(X=exprs(mainz)[probesAffy,tmpST], MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(mainz)[tmpST,"t.dmfs"], z=pData(mainz)[tmpST,"e.dmfs"]))

tmpST <- as.character(subtypeData[["transbig"]]$subtype2) == subT
cindexTransbig <- t(apply(X=exprs(transbig)[probesAffy,tmpST], MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(transbig)[tmpST,"t.dmfs"], z=pData(transbig)[tmpST,"e.dmfs"]))

tmpST <- as.character(subtypeData[["upp"]]$subtype2) == subT
cindexUpp <- t(apply(X=exprs(upp)[probesAffy,tmpST], MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(upp)[tmpST,"t.rfs"], z=pData(upp)[tmpST,"e.rfs"]))

tmpST <- as.character(subtypeData[["unt"]]$subtype2) == subT
cindexUnt <- t(apply(X=exprs(unt)[probesAffy,tmpST], MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(unt)[tmpST,"t.dmfs"], z=pData(unt)[tmpST,"e.dmfs"]))

tmpST <- as.character(subtypeData[["vdx"]]$subtype2) == subT
cindexVdx <- t(apply(X=exprs(vdx)[probesAffy,tmpST], MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(vdx)[tmpST,"t.dmfs"], z=pData(vdx)[tmpST,"e.dmfs"]))

tmpST <- as.character(subtypeData[["nki"]]$subtype2) == subT
cindexNki <- t(apply(X=exprs(nki)[probesAgi,tmpST], MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(nki)[tmpST,"t.dmfs"], z=pData(nki)[tmpST,"e.dmfs"]))                         
    

###################################################
## combining the cindices for each gene from all datasets
###################################################
tt <- as.data.frame(NULL)
for(i in 1:length(geneSymbols)){
  tt <- rbind(
    tt,combine.est(x=cbind(   cindexMainz[i,"cindex"],
                              cindexTransbig[i,"cindex"],
                              cindexUpp[i,"cindex"],
                              cindexUnt[i,"cindex"],
                              cindexVdx[i,"cindex"],
                              cindexNki[i,"cindex"]),
                   x.se=cbind(cindexMainz[i,"cindex"],
                              cindexTransbig[i,"cindex.se"],
                              cindexUpp[i,"cindex.se"],
                              cindexUnt[i,"cindex.se"],
                              cindexVdx[i,"cindex.se"],
                              cindexNki[i,"cindex.se"]),)
              )
}
tt$lower <- tt$estimate + qnorm(0.025, lower.tail=TRUE) * tt$se
tt$upper <- tt$estimate + qnorm(0.025, lower.tail=FALSE) * tt$se
rownames(tt) <- geneSymbols
colnames(tt) <- c("cindex","cindex.se","lower","upper")
ccindex <- tt

dev.off()

 
 
 
           
#######################
## CINDEX: every SINGLE gene and the different datasets as forestplot
#######################
##pdf("cindex_forestplot-for-each-gene-showing-all-datasets.pdf")
for(i in 1:length(dataset.list)) {
  myspace <- " "
  tt <- rbind(cindexall.mainz.small[i,],
              cindexall.transbig.small[i,],
              cindexall.upp.small[i,],
              cindexall.unt.small[i,],
              cindexall.vdx.small[i,],
              cindexall.nki.small[i,],
              as.numeric(ccindex[i,]))

  rownames(tt) <- dataset.list
  tt <- as.data.frame(tt)
  labeltext <- cbind(c("Dataset",dataset.list),c(rep(myspace,length(dataset.list))))
  bs <- rep(0.5, nrow(labeltext))

  r.mean <- c(NA,tt$cindex)
  r.lower <- c(NA,tt$cindex + qnorm(0.025, lower.tail=TRUE) * tt$cindex.se)
  r.upper <- c(NA,tt$cindex + qnorm(0.025, lower.tail=FALSE) * tt$cindex.se)
  
  x.ticks.lower <- (floor((min(r.mean,na.rm=TRUE) - 0.1) * 10)/10)
  x.ticks.upper <- (floor((max(r.mean,na.rm=TRUE) + 0.2) * 10)/10)
   
  forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(x.ticks.lower,x.ticks.upper,0.05), xlab=paste(gs.list[i], myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(0.3,0.8),is.summary=(c(rep(FALSE,7),TRUE)))
  title(paste("cindex forestplot, ", gs.list[i]))
}
##dev.off()
