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
### loadind data
###################################################
data(breastCancerData)
data(mainz,transbig,upp,unt,vdx,nki)
data(scmgene)

###################################################
### create variables
###################################################
gsList <- tolower(fData(mainz7g)[,"Gene.symbol"])
gidList <- fData(mainz7g)[,"Gene.ID"]
probesNKI <- as.character(fData(nki7g)[,"probe"])
probesAffy <- fData(mainz7g)[,"probe"]
datasetList <- c("MAINZ","TRANSBIG","UPP","UNT","VDX","NKI","","Overall")
myspace <- " "
mybigspace <- "    "
tc <- 10 * 365  


###################################################
### get subtypes
###################################################
scmgene.mainz <- subtype.cluster.predict(sbt.model=scmgene, data=t(exprs(mainz)), annot=fData(mainz), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
scmgene.transbig <- subtype.cluster.predict(sbt.model=scmgene, data=t(exprs(transbig)), annot=fData(transbig), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
scmgene.upp <- subtype.cluster.predict(sbt.model=scmgene, data=t(exprs(upp)), annot=fData(upp), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
scmgene.unt <- subtype.cluster.predict(sbt.model=scmgene, data=t(exprs(unt)), annot=fData(unt), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
scmgene.vdx <- subtype.cluster.predict(sbt.model=scmgene, data=t(exprs(vdx)), annot=fData(vdx), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
scmgene.nki <- subtype.cluster.predict(sbt.model=scmgene, data=t(exprs(nki)), annot=fData(nki), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)


###################################################
### save plots
###################################################
pdf("KMPlot-STAT1-all.pdf")


###################################################
### survival curves using subtypes, all datasets
###################################################
surv.time.all <- c(pData(mainz)[ ,"t.dmfs"], pData(transbig)[ ,"t.dmfs"], pData(unt)[ ,"t.dmfs"], pData(upp)[ ,"t.rfs"], pData(vdx)[ ,"t.dmfs"], pData(nki)[ ,"t.dmfs"])  
surv.event.all <- c(pData(mainz)[ ,"e.dmfs"], pData(transbig)[ ,"e.dmfs"], pData(unt)[ ,"e.dmfs"], pData(upp)[ ,"e.rfs"], pData(vdx)[ ,"e.dmfs"], pData(nki)[ ,"e.dmfs"]) 

mygroup <- c(as.character(scmgene.mainz$subtype2), as.character(scmgene.transbig$subtype2), as.character(scmgene.upp$subtype2), as.character(scmgene.unt$subtype2), as.character(scmgene.vdx$subtype2), as.character(scmgene.nki$subtype2))
mygroup <- as.factor(mygroup)

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=as.numeric(mygroup))

km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Prob. of survival for subgroups of patients according to subtypes", sub.title="", leg.text=c(levels(mygroup)), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred", "orange"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))


###################################################
### survival curves STAT1, quantiles, all datasets
###################################################
stat1Gs <- "stat1"
stat1Gid <- 6772
stat1Paf <- "209969_s_at"
stat1Pagi <- "NM_007315"

surv.time.all <- c(pData(mainz)[ ,"t.dmfs"], pData(transbig)[ ,"t.dmfs"], pData(unt)[ ,"t.dmfs"], pData(upp)[ ,"t.rfs"], pData(vdx)[ ,"t.dmfs"], pData(nki)[ ,"t.dmfs"])
surv.event.all <- c(pData(mainz)[ ,"e.dmfs"], pData(transbig)[ ,"e.dmfs"], pData(unt)[ ,"e.dmfs"], pData(upp)[ ,"e.rfs"], pData(vdx)[ ,"e.dmfs"], pData(nki)[ ,"e.dmfs"])
stat1.exprs <- c(exprs(mainz)[stat1Paf,], exprs(transbig)[stat1Paf,], exprs(unt)[stat1Paf,], exprs(upp)[stat1Paf,], exprs(vdx)[stat1Paf,], exprs(nki)[stat1Pagi,])
stat1.exprs.length <- c(length(exprs(mainz)[stat1Paf,]), length(exprs(transbig)[stat1Paf,]), length(exprs(unt)[stat1Paf,]), length(exprs(upp)[stat1Paf,]), length(exprs(vdx)[stat1Paf,]), length(exprs(nki)[stat1Pagi,]))

pos <- 0
mygroup <- NULL
for(i in stat1.exprs.length){
  qq <- stat1.exprs[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.33, 0.66), na.rm=TRUE)
  qq[stat1.exprs[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[stat1.exprs[(pos+1):(pos+i)] >= myq[1] & stat1.exprs[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[stat1.exprs[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup <- c(mygroup,qq)
  pos <- pos + i
}                     

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)

km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Probability of survival for STAT1", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))



###################################################
### survival curves STAT1 using subtypes, all datasets
###################################################
stat1Gs <- "stat1"
stat1Gid <- 6772
stat1Paf <- "209969_s_at"
stat1Pagi <- "NM_007315"

subtypes.all <- c(as.character(scmgene.mainz$subtype2), as.character(scmgene.transbig$subtype2), as.character(scmgene.upp$subtype2), as.character(scmgene.unt$subtype2), as.character(scmgene.vdx$subtype2), as.character(scmgene.nki$subtype2))

subtypes.names <- levels(as.factor(subtypes.all))
subtypes.all[is.na(subtypes.all)] <- F

subtypes.all.length <- NULL
for(i in subtypes.names){
  subtypes.all.length <- c(length(na.omit(scmgene.mainz$subtype2[scmgene.mainz$subtype2==i])),
                         length(na.omit(scmgene.transbig$subtype2[scmgene.transbig$subtype2==i])),
                         length(na.omit(scmgene.upp$subtype2[scmgene.upp$subtype2==i])),
                         length(na.omit(scmgene.unt$subtype2[scmgene.unt$subtype2==i])),
                         length(na.omit(scmgene.vdx$subtype2[scmgene.vdx$subtype2==i])),
                         length(na.omit(scmgene.nki$subtype2[scmgene.nki$subtype2==i])))

  stype <- subtypes.all==i

  surv.time.all <- c(pData(mainz)[ ,"t.dmfs"], pData(transbig)[ ,"t.dmfs"], pData(unt)[ ,"t.dmfs"], pData(upp)[ ,"t.rfs"], pData(vdx)[ ,"t.dmfs"], pData(nki)[ ,"t.dmfs"])
  surv.time.all <- surv.time.all[stype]

  surv.event.all <- c(pData(mainz)[ ,"e.dmfs"], pData(transbig)[ ,"e.dmfs"], pData(unt)[ ,"e.dmfs"], pData(upp)[ ,"e.rfs"], pData(vdx)[ ,"e.dmfs"], pData(nki)[ ,"e.dmfs"]) 
  surv.event.all <- surv.event.all[stype]

  stat1.exprs <- c(exprs(mainz)[stat1Paf,], exprs(transbig)[stat1Paf,], exprs(unt)[stat1Paf,], exprs(upp)[stat1Paf,], exprs(vdx)[stat1Paf,], exprs(nki)[stat1Pagi,])
  stat1.exprs <- stat1.exprs[stype]                    

  pos <- 0
  mygroup <- NULL
  for(j in subtypes.all.length){
    qq <- stat1.exprs[(pos+1):(pos+j)]
    myq <- quantile(qq, probs=c(0.33, 0.66), na.rm=TRUE)
    qq[stat1.exprs[(pos+1):(pos+j)] < myq[1]] <- 1
    qq[stat1.exprs[(pos+1):(pos+j)] >= myq[1] & stat1.exprs[(pos+1):(pos+j)] < myq[2]] <- 2
    qq[stat1.exprs[(pos+1):(pos+j)] > myq[2]] <- 3
    qq <- factor(x=qq, levels=1:3)
    mygroup <- c(mygroup,qq)
    pos <- pos + j
  }                     

  surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)

  dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
  
  km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title=c("Probability of survival for STAT1, ", i), sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))
}

dev.off()
