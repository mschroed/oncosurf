###################################################
## load libraries
###################################################
library(survcomp)
library(Biobase)
library(genefu)
library(Hmisc)
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


###################################################
## specify genes
###################################################
##Oncotype DX
##gs.list <- c("ACTB","GAPDH","GUS","RPLPO","TFRC","GRB7","ER","HER2","PGR","BCL2","SCUBE2","Survivin","KI67","MYBL2","CCNB1","STK15","CTSL2","MMP11","GSTM1","BAG1","CD68")
gid.list <- c(7037,60,5241,4320,4605,573,2597,968,891,2944,1515,2886,596,57758)
geneid.list <- c("geneid.7037","geneid.60","geneid.5241","geneid.4320","geneid.4605","geneid.573","geneid.2597","geneid.891","geneid.2944","geneid.1515","geneid.2886","geneid.596","geneid.57758")
gs.list <- c("TFRC","PGR","MMP11","MYBL2","BAG1","GAPDH","CCNB1","GSTM1","CTSL2","GRB7","BCL2","SCUBE2")
## gs.list <- c("TFRC","ACTB","PGR","MMP11","MYBL2","BAG1","GAPDH","<NA>","CCNB1","GSTM1","CTSL2","GRB7","BCL2","SCUBE2")
## geneid.968 has no gene symbol in the annotation data (for mainz)

# Gene Symbols Oncotype DX
# KI67 STK15 BIRC5 CCNB1 MYBL2 MMP11 CTSL2 GRB7 HER2 ER PGR BCL2 SCUBE2 GSTM1 BAG1 CD68 ACTB GAPDH RPLP0 GUS TFRC
odx <- c("KI67", "STK15", "BIRC5", "CCNB1", "MYBL2", "MMP11", "CTSL2", "GRB7", "HER2", "ER", "PGR", "BCL2", "SCUBE2", "GSTM1", "BAG1", "CD68", "ACTB", "GAPDH", "RPLP0", "GUS", "TFRC")





###################################################
## concordance.index computation
###################################################
cindexall.mainz <- t(apply(X=exprs(mainz), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainz)[ ,"t.dmfs"], z=pData(mainz)[ ,"e.dmfs"]))

cindexall.transbig <- t(apply(X=exprs(transbig), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbig)[ ,"t.dmfs"], z=pData(transbig)[ ,"e.dmfs"]))

cindexall.upp <- t(apply(X=exprs(upp), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(upp)[ ,"t.rfs"], z=pData(upp)[ ,"e.rfs"]))

cindexall.unt <- t(apply(X=exprs(unt), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(unt)[ ,"t.dmfs"], z=pData(unt)[ ,"e.dmfs"]))

cindexall.vdx <- t(apply(X=exprs(vdx), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdx)[ ,"t.dmfs"], z=pData(vdx)[ ,"e.dmfs"])) 

cindexall.nki <- t(apply(X=exprs(nki), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nki)[ ,"t.dmfs"], z=pData(nki)[ ,"e.dmfs"]))



###################################################
## D.index computation
###################################################
dindexall.mainz <- t(apply(X=exprs(mainz), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainz)[ ,"t.dmfs"], z=pData(mainz)[ ,"e.dmfs"]))

dindexall.transbig <- t(apply(X=exprs(transbig), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbig)[ ,"t.dmfs"], z=pData(transbig)[ ,"e.dmfs"]))

dindexall.upp <- t(apply(X=exprs(upp), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(upp)[ ,"t.rfs"], z=pData(upp)[ ,"e.rfs"]))

dindexall.unt <- t(apply(X=exprs(unt), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(unt)[ ,"t.dmfs"], z=pData(unt)[ ,"e.dmfs"]))

dindexall.vdx <- t(apply(X=exprs(vdx), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdx)[ ,"t.dmfs"], z=pData(vdx)[ ,"e.dmfs"]))

dindexall.nki <- t(apply(X=exprs(nki), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nki)[ ,"t.dmfs"], z=pData(nki)[ ,"e.dmfs"]))



###################################################
### rescale function
###################################################
rescale <- function(x, na.rm=FALSE, q=0.05) {
  ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
  mi <- quantile(x, probs=q/2, na.rm=na.rm)
  x <- (x - mi) / (ma - mi)
  return((x - 0.5) * 2)
}


###################################################
### hazard ratio computation
###################################################
hratio.mainz <- t(apply(X=rescale(exprs(mainz) , q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(mainz)[ ,"t.dmfs"], z=pData(mainz)[ ,"e.dmfs"]))

hratio.transbig <- t(apply(X=rescale(exprs(transbig), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(transbig)[ ,"t.dmfs"], z=pData(transbig)[ ,"e.dmfs"]))

hratio.upp <- t(apply(X=rescale(exprs(upp), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(upp)[ ,"t.rfs"], z=pData(upp)[ ,"e.rfs"]))

hratio.unt <- t(apply(X=rescale(exprs(unt), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(unt)[ ,"t.dmfs"], z=pData(unt)[ ,"e.dmfs"]))

hratio.vdx <- t(apply(X=rescale(exprs(vdx), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(vdx)[ ,"t.dmfs"], z=pData(vdx)[ ,"e.dmfs"]))

hratio.nki <- t(apply(X=rescale(exprs(nki), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(nki)[ ,"t.dmfs"], z=pData(nki)[ ,"e.dmfs"]))






###################################################
### loadind data
###################################################
data(breastCancerData)
                
                
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
### compute C indices for subsets
###################################################
cindexall.mainz.small <- t(apply(X=exprs(mainz7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(mainz7g)[ ,"t.dmfs"], z=pData(mainz7g)[ ,"e.dmfs"]))


cindexall.transbig.small <- t(apply(X=exprs(transbig7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(transbig7g)[ ,"t.dmfs"], z=pData(transbig7g)[ ,"e.dmfs"]))

cindexall.vdx.small <- t(apply(X=exprs(vdx7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(vdx7g)[ ,"t.dmfs"], z=pData(vdx7g)[ ,"e.dmfs"])) 

cindexall.upp.small <- t(apply(X=exprs(upp7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(upp7g)[ ,"t.rfs"], z=pData(upp7g)[ ,"e.rfs"]))

cindexall.unt.small <- t(apply(X=exprs(unt7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(unt7g)[ ,"t.dmfs"], z=pData(unt7g)[ ,"e.dmfs"]))

cindexall.nki.small <- t(apply(X=exprs(nki7g), MARGIN=1, function(x, y, z) {
    tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
    return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(nki7g)[ ,"t.dmfs"], z=pData(nki7g)[ ,"e.dmfs"]))


###################################################
### compute D indices for subsets
###################################################
dindexall.mainz.small <- t(apply(X=exprs(mainz7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(mainz7g)[ ,"t.dmfs"], z=pData(mainz7g)[ ,"e.dmfs"]))

dindexall.transbig.small <- t(apply(X=exprs(transbig7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(transbig7g)[ ,"t.dmfs"], z=pData(transbig7g)[ ,"e.dmfs"]))

dindexall.upp.small <- t(apply(X=exprs(upp7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(upp7g)[ ,"t.rfs"], z=pData(upp7g)[ ,"e.rfs"]))

dindexall.unt.small <- t(apply(X=exprs(unt7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(unt7g)[ ,"t.dmfs"], z=pData(unt7g)[ ,"e.dmfs"]))

dindexall.vdx.small <- t(apply(X=exprs(vdx7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(vdx7g)[ ,"t.dmfs"], z=pData(vdx7g)[ ,"e.dmfs"]))

dindexall.nki.small <- t(apply(X=exprs(nki7g), MARGIN=1, function(x, y, z) {
    tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(nki7g)[ ,"t.dmfs"], z=pData(nki7g)[ ,"e.dmfs"]))


###################################################
### create rescaling function
###################################################
rescale <- function(x, na.rm=FALSE, q=0.05) {
  ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
  mi <- quantile(x, probs=q/2, na.rm=na.rm)
  x <- (x - mi) / (ma - mi)
  return(x) ##(x - 0.5) * 2)
}


###################################################
### compute hazard ratios for subsets
###################################################
hratio.mainz.small <- t(apply(X=rescale(exprs(mainz7g) , q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(mainz7g)[ ,"t.dmfs"], z=pData(mainz7g)[ ,"e.dmfs"]))

hratio.transbig.small <- t(apply(X=rescale(exprs(transbig7g), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(transbig7g)[ ,"t.dmfs"], z=pData(transbig7g)[ ,"e.dmfs"]))

hratio.upp.small <- t(apply(X=rescale(exprs(upp7g), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(upp7g)[ ,"t.rfs"], z=pData(upp7g)[ ,"e.rfs"]))

hratio.unt.small <- t(apply(X=rescale(exprs(unt7g), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(unt7g)[ ,"t.dmfs"], z=pData(unt7g)[ ,"e.dmfs"]))

hratio.vdx.small <- t(apply(X=rescale(exprs(vdx7g), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(vdx7g)[ ,"t.dmfs"], z=pData(vdx7g)[ ,"e.dmfs"]))

hratio.nki.small <- t(apply(X=rescale(exprs(nki7g), q=0.05, na.rm=TRUE), MARGIN=1, function(x, y, z) {
    tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
    return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); },
    y=pData(nki7g)[ ,"t.dmfs"], z=pData(nki7g)[ ,"e.dmfs"]))


###################################################
### combine Cindices
###################################################
tt <- as.data.frame(NULL)
for(i in 1:7){
  tt <- rbind(
    tt,combine.est(
    x=cbind(   cindexall.mainz.small[i,"cindex"],
        cindexall.transbig.small[i,"cindex"],
        cindexall.upp.small[i,"cindex"],
        cindexall.unt.small[i,"cindex"],
        cindexall.vdx.small[i,"cindex"],
        cindexall.nki.small[i,"cindex"]),
    x.se=cbind(cindexall.mainz.small[i,"cindex.se"],
        cindexall.transbig.small[i,"cindex.se"],
        cindexall.upp.small[i,"cindex.se"],
        cindexall.unt.small[i,"cindex.se"],
        cindexall.vdx.small[i,"cindex.se"],
        cindexall.nki.small[i,"cindex.se"]),)
    )
}
tt$lower <- tt$estimate + qnorm(0.025, lower.tail=TRUE) * tt$se
tt$upper <- tt$estimate + qnorm(0.025, lower.tail=FALSE) * tt$se
rownames(tt) <- gsList
colnames(tt) <- c("cindex","cindex.se","lower","upper")
ccindex.small <- tt


###################################################
### combine D indices
###################################################
tt <- as.data.frame(NULL)
for(i in 1:7){
  tt <- rbind(
    tt,combine.est(
    x=cbind( dindexall.mainz.small[i,"dindex"],
        dindexall.transbig.small[i,"dindex"],
        dindexall.upp.small[i,"dindex"],
        dindexall.unt.small[i,"dindex"],
        dindexall.vdx.small[i,"dindex"],
        dindexall.nki.small[i,"dindex"]),
    x.se=cbind(dindexall.mainz.small[i,"dindex.se"],
        dindexall.transbig.small[i,"dindex.se"],
        dindexall.upp.small[i,"dindex.se"],
        dindexall.unt.small[i,"dindex.se"],
        dindexall.vdx.small[i,"dindex.se"],
        dindexall.nki.small[i,"dindex.se"]),)
    )
}
tt$lower <- tt$estimate + qnorm(0.025, lower.tail=TRUE) * tt$se
tt$upper <- tt$estimate + qnorm(0.025, lower.tail=FALSE) * tt$se
rownames(tt) <- gsList
colnames(tt) <- c("dindex","dindex.se","lower","upper")
cdindex.small <- tt


###################################################
### chunk number 21: combineHratio
###################################################
#line 399 "survcomp.Rnw"
tt <- as.data.frame(NULL)
for(i in 1:7){
  tt <- rbind(
    tt,combine.est(
    x=cbind( hratio.mainz.small[i,"hratio"],
        hratio.transbig.small[i,"hratio"],
        hratio.upp.small[i,"hratio"],
        hratio.unt.small[i,"hratio"],
        hratio.vdx.small[i,"hratio"],
        hratio.nki.small[i,"hratio"]),
    x.se=cbind(hratio.mainz.small[i,"hratio.se"],
        hratio.transbig.small[i,"hratio.se"],
        hratio.upp.small[i,"hratio.se"],
        hratio.unt.small[i,"hratio.se"],
        hratio.vdx.small[i,"hratio.se"],
        hratio.nki.small[i,"hratio.se"]),)
    )
}
tt$lower <- tt$estimate + qnorm(0.025, lower.tail=TRUE) * tt$se
tt$upper <- tt$estimate + qnorm(0.025, lower.tail=FALSE) * tt$se
rownames(tt) <- gsList
colnames(tt) <- c("hratio","hratio.se","lower","upper")
chratio.small <- tt


#######################
## HRATIO: AURKA and different datasets as forestplot
#######################
##pdf("survcomp-app-note-figure_hratio.pdf")
tt <- rbind(hratio.mainz.small[3,],
            hratio.transbig.small[3,],
            hratio.upp.small[3,],
            hratio.unt.small[3,],
            hratio.vdx.small[3,],
            hratio.nki.small[3,],
            NA,
            as.numeric(chratio.small[3,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
##tt <- log2(tt)
labeltext <- cbind(c("Dataset",datasetList),c(rep(mybigspace,9)))
bs <- rep(0.5, nrow(labeltext))
r.mean <- c(NA, log2(tt$hratio))
r.lower <- c(NA, log2(tt$hratio + qnorm(0.025, lower.tail=TRUE) * tt$hratio.se))
r.upper <- c(NA, log2(tt$hratio + qnorm(0.025, lower.tail=FALSE) * tt$hratio.se))
  
forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,5.5,0.5), xlab=paste("Hazard Ratio", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.75,6),is.summary=(c(rep(FALSE,8),TRUE)))
##title(paste("log2 hratio forestplot, AURKA"))
dev.copy2eps(file="figure1c.eps")




###################################################
### forestplot all seven genes, each datasets
###################################################
### first (MAINZ and TRANSBIG)
tt <- rbind(cindexall.mainz.small[1,],
            cindexall.mainz.small[2,],
            cindexall.mainz.small[3,],
            cindexall.mainz.small[4,],
            cindexall.mainz.small[5,],                                    
            cindexall.mainz.small[6,],
            cindexall.mainz.small[7,],
            NA,
            cindexall.transbig.small[1,],            
            cindexall.transbig.small[2,],
            cindexall.transbig.small[3,],
            cindexall.transbig.small[4,],
            cindexall.transbig.small[5,],                                    
            cindexall.transbig.small[6,],
            cindexall.transbig.small[7,])

rownames(tt) <- c("MAINZ1", "MAINZ2", "MAINZ3", "MAINZ4","MAINZ5", "MAINZ6", "MAINZ7", "a",
                  "TRANSBIG1", "TRANSBIG2", "TRANSBIG3", "TRANSBIG4", "TRANSBIG5", "TRANSBIG6", "TRANSBIG7")
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",
                     "MAINZ", NA, NA, NA, NA, NA, NA, NA,
                     "TRANSBIG", NA, NA, NA, NA, NA, NA),
                   c("Gene",
                     rep(c(gsList,NA), length(datasetList)-7), c(gsList)),
                   c(rep(mybigspace,length(16))))
bs <- rep(0.5, nrow(labeltext))   
r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.3,0.8,0.1),
    col=meta.colors(line=c(rep(c(NA, rainbow(7)),7)), zero="firebrick",
    box=c(rep(c(NA,rainbow(7)),7))), box.size=bs, clip=c(0.3,0.8), is.summary=FALSE)##(c(rep(FALSE,length(rownames(tt))-6), c(rep(TRUE,6)))))
dev.copy2eps(file="forestplot-MT.eps")    
    
    
###########################
### second (UPP and UNT)
tt <- rbind(cindexall.upp.small[1,],            
            cindexall.upp.small[2,],
            cindexall.upp.small[3,],
            cindexall.upp.small[4,],
            cindexall.upp.small[5,],                                    
            cindexall.upp.small[6,],
            cindexall.upp.small[7,],
            NA,
            cindexall.unt.small[1,],            
            cindexall.unt.small[2,],
            cindexall.unt.small[3,],
            cindexall.unt.small[4,],
            cindexall.unt.small[5,],                                    
            cindexall.unt.small[6,],
            cindexall.unt.small[7,])

rownames(tt) <- c("UPP1", "UPP2", "UPP3", "UPP4", "UPP5", "UPP6", "UPP7", "c",
                  "UNT1", "UNT2", "UNT3", "UNT4", "UNT5", "UNT6", "UNT7")
                  
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",
                     "UPP", NA, NA, NA, NA, NA, NA, NA,
                     "UNT", NA, NA, NA, NA, NA, NA),
                   c("Gene",
                     rep(c(gsList,NA), length(datasetList)-7), c(gsList)),
                   c(rep(mybigspace,length(16))))
bs <- rep(0.5, nrow(labeltext))   
r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.3,0.8,0.1),
    col=meta.colors(line=c(rep(c(NA, rainbow(7)),7)), zero="firebrick",
    box=c(rep(c(NA,rainbow(7)),7))), box.size=bs, clip=c(0.3,0.8), is.summary=FALSE)  
dev.copy2eps(file="forestplot-UU.eps")        
    
    
#############################
## third (VDX and NKI)
tt <- rbind(cindexall.vdx.small[1,],            
            cindexall.vdx.small[2,],
            cindexall.vdx.small[3,],
            cindexall.vdx.small[4,],
            cindexall.vdx.small[5,],                                    
            cindexall.vdx.small[6,],
            cindexall.vdx.small[7,],
            NA,
            cindexall.nki.small[1,],            
            cindexall.nki.small[2,],
            cindexall.nki.small[3,],
            cindexall.nki.small[4,],
            cindexall.nki.small[5,],                                    
            cindexall.nki.small[6,],
            cindexall.nki.small[7,])

rownames(tt) <- c("VDX1", "VDX2", "VDX3", "VDX4", "VDX5", "VDX6", "VDX7", "e",
                  "NKI1", "NKI2", "NKI3", "NKI4", "NKI5", "NKI6", "NKI7")
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",
                     "VDX", NA, NA, NA, NA, NA, NA, NA,
                     "NKI", NA, NA, NA, NA, NA, NA),
                   c("Gene",
                     rep(c(gsList,NA), length(datasetList)-7), c(gsList)),
                   c(rep(mybigspace,length(16))))
bs <- rep(0.5, nrow(labeltext))   
r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.3,0.8,0.1),
    col=meta.colors(line=c(rep(c(NA, rainbow(7)),7)), zero="firebrick",
    box=c(rep(c(NA,rainbow(7)),7))), box.size=bs, clip=c(0.3,0.8), is.summary=FALSE)##(c(rep(FALSE,length(rownames(tt))-6), c(rep(TRUE,6)))))
dev.copy2eps(file="forestplot-VN.eps")     


################################
## overall (combined estimations)

tt <- rbind(as.numeric(ccindex.small[1,]), 
            as.numeric(ccindex.small[2,]),             
            as.numeric(ccindex.small[3,]), 
            as.numeric(ccindex.small[4,]), 
            as.numeric(ccindex.small[5,]),
            as.numeric(ccindex.small[6,]),                                   
            as.numeric(ccindex.small[7,]))

rownames(tt) <- c("ALL1", "ALL2", "ALL3", "ALL4", "ALL5", "ALL6", "ALL7")
tt <- as.data.frame(tt)
colnames(tt) <- c("cindex", "cindex.se", "lower", "upper")
labeltext <- cbind(c("Dataset",
                     "Overall", NA, NA, NA, NA, NA, NA),
                   c("Gene",
                     rep(c(gsList))),
                   c(rep(mybigspace,length(gsList)+1)))
bs <- rep(0.5, nrow(labeltext))   
r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.7,0.1),
    col=meta.colors(line=c(rep(c(NA, rainbow(7)),7)), zero="firebrick",
    box=c(rep(c(NA,rainbow(7)),7))), box.size=bs, clip=c(0.4,0.7), is.summary=TRUE)##(c(rep(FALSE,length(rownames(tt))-6), c(rep(TRUE,6)))))
dev.copy2eps(file="forestplot-overall.eps")
        


################################
## overall (combined estimations) WITH SPACE

tt <- rbind(as.numeric(ccindex.small[1,]),NA, 
as.numeric(ccindex.small[2,]),NA,             
as.numeric(ccindex.small[3,]),NA, 
as.numeric(ccindex.small[4,]),NA, 
as.numeric(ccindex.small[5,]),NA,
as.numeric(ccindex.small[6,]),NA,                                   
as.numeric(ccindex.small[7,]))

rownames(tt) <- c("ALL1","", "ALL2","", "ALL3","", "ALL4","", "ALL5","", "ALL6","", "ALL7")
tt <- as.data.frame(tt)
colnames(tt) <- c("cindex", "cindex.se", "lower", "upper")
labeltext <- cbind(c("Dataset",
"Overall","","","","","","", NA, NA, NA, NA, NA, NA),
c("Gene","esr1",NA,"erbb2",NA,"aurka",NA,"plau",NA,"vegfa",NA,"stat1",NA,"casp3"),
c(rep(mybigspace,14)))
bs <- rep(0.5, nrow(labeltext))   
r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.7,0.1),
col=meta.colors(line=c(rep(c(NA, rainbow(7)),7)), zero="firebrick",
box=c(rep(c(NA,rainbow(7)),7))), box.size=bs, clip=c(0.4,0.7), is.summary=TRUE)

dev.copy2eps(file="forestplot-overall-spaces.eps")





#forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
#    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.8,0.05), xlab=paste("AURKA and VEGF Concordance Index", myspace, sep=""),
#    col=meta.colors(line=c(rep(c(NA, "darkblue", "seagreen", "deeppink", "darkslategray", "khaki", "orangered", "maroon"),7)), zero="firebrick",
#    box=c(rep(c(NA," royalblue", "forestgreen", "hotpink", "slategray", "palegoldenrod", "tomato", "saddlebrown"),7))), box.size=bs,
#    clip=c(0.3,1), is.summary=FALSE)##(c(rep(FALSE,length(rownames(tt))-6), c(rep(TRUE,6)))))





######################################################################################################
## gsList <- c("ESR1", "ERBB2", "AURKA", "PLAU", "VEGF", "STAT1", "CASP3")
######################################################################################################
### forestplot all seven genes, each datasets
######################################################################################################
datasetListSmall <- c("MAINZ", "TRANSBIG", "UPP", "UNT", "VDX", "NKI")
### ESR1
tt <- rbind(cindexall.mainz.small[1,],
            cindexall.transbig.small[1,],
            cindexall.upp.small[1,],
            cindexall.unt.small[1,],
            cindexall.vdx.small[1,],
            cindexall.nki.small[1,])

rownames(tt) <- datasetListSmall
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetListSmall),c(mybigspace,rep(mybigspace,length(datasetListSmall))))
bs <- rep(0.5, nrow(labeltext))   

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            
                                   
forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.3,0.7,0.1),
    col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"),
    box.size=bs, clip=c(0.3,0.8), is.summary=FALSE, xlab="ESR1 Concordance Index")


dev.copy2eps(file="forestplot-ESR1.eps")
 

### ERBB2
tt <- rbind(cindexall.mainz.small[2,],
            cindexall.transbig.small[2,],
            cindexall.upp.small[2,],
            cindexall.unt.small[2,],
            cindexall.vdx.small[2,],
            cindexall.nki.small[2,])

rownames(tt) <- datasetListSmall
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetListSmall),c(mybigspace,rep(mybigspace,length(datasetListSmall))))
bs <- rep(0.5, nrow(labeltext))   

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            
                                   
forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.7,0.1),
    col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"),
    box.size=bs, clip=c(0.3,0.8), is.summary=FALSE, xlab="ERBB2 Concordance Index")


dev.copy2eps(file="forestplot-ERBB2.eps")
   

### AURKA
tt <- rbind(cindexall.mainz.small[3,],
            cindexall.transbig.small[3,],
            cindexall.upp.small[3,],
            cindexall.unt.small[3,],
            cindexall.vdx.small[3,],
            cindexall.nki.small[3,])

rownames(tt) <- datasetListSmall
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetListSmall),c(mybigspace,rep(mybigspace,length(datasetListSmall))))
bs <- rep(0.5, nrow(labeltext))   

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            
                                   
forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.8,0.1),
    col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"),
    box.size=bs, clip=c(0.3,0.8), is.summary=FALSE, xlab="AURKA Concordance Index")


dev.copy2eps(file="forestplot-AURKA.eps")
 


### PLAU
tt <- rbind(cindexall.mainz.small[4,],
            cindexall.transbig.small[4,],
            cindexall.upp.small[4,],
            cindexall.unt.small[4,],
            cindexall.vdx.small[4,],
            cindexall.nki.small[4,])

rownames(tt) <- datasetListSmall
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetListSmall),c(mybigspace,rep(mybigspace,length(datasetListSmall))))
bs <- rep(0.5, nrow(labeltext))   

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            
                                   
forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.3,0.7,0.1),
    col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"),
    box.size=bs, clip=c(0.3,0.8), is.summary=FALSE, xlab="PLAU Concordance Index")


dev.copy2eps(file="forestplot-PLAU.eps")
 

### VEGFA
tt <- rbind(cindexall.mainz.small[5,],
            cindexall.transbig.small[5,],
            cindexall.upp.small[5,],
            cindexall.unt.small[5,],
            cindexall.vdx.small[5,],
            cindexall.nki.small[5,])

rownames(tt) <- datasetListSmall
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetListSmall),c(mybigspace,rep(mybigspace,length(datasetListSmall))))
bs <- rep(0.5, nrow(labeltext))   

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            
                                   
forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.7,0.1),
    col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"),
    box.size=bs, clip=c(0.3,0.8), is.summary=FALSE, xlab="VEGFA Concordance Index")


dev.copy2eps(file="forestplot-VEGFA.eps")


### STAT1
tt <- rbind(cindexall.mainz.small[6,],
            cindexall.transbig.small[6,],
            cindexall.upp.small[6,],
            cindexall.unt.small[6,],
            cindexall.vdx.small[6,],
            cindexall.nki.small[6,])

rownames(tt) <- datasetListSmall
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetListSmall),c(mybigspace,rep(mybigspace,length(datasetListSmall))))
bs <- rep(0.5, nrow(labeltext))   

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            
                                   
forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.7,0.1),
    col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"),
    box.size=bs, clip=c(0.3,0.8), is.summary=FALSE, xlab="STAT1 Concordance Index")


dev.copy2eps(file="forestplot-STAT1.eps")
 


### CASP3
tt <- rbind(cindexall.mainz.small[7,],
            cindexall.transbig.small[7,],
            cindexall.upp.small[7,],
            cindexall.unt.small[7,],
            cindexall.vdx.small[7,],
            cindexall.nki.small[7,])

rownames(tt) <- datasetListSmall
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetListSmall),c(mybigspace,rep(mybigspace,length(datasetListSmall))))
bs <- rep(0.5, nrow(labeltext))   

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            
                                   
forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.7,0.1),
    col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"),
    box.size=bs, clip=c(0.3,0.8), is.summary=FALSE, xlab="CASP3 Concordance Index")


dev.copy2eps(file="forestplot-CASP3.eps")
 



### Overall
tt <- rbind(ccindex.small[1,],
            ccindex.small[2,],
            ccindex.small[3,],
            ccindex.small[4,],
            ccindex.small[5,],
            ccindex.small[6,],
            ccindex.small[7,])

rownames(tt) <- gsList
tt <- as.data.frame(tt)
labeltext <- cbind(c("Gene name",gsList),c(rep(mybigspace,length(gsList)+1)))
bs <- rep(0.5, nrow(labeltext))   

r.mean <- c(NA,tt$cindex)
r.lower <- c(NA,tt$lower)
r.upper <- c(NA,tt$upper)                            
                                   
forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
    align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.7,0.1),
    col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"),
    box.size=bs, clip=c(0.3,0.8), is.summary=TRUE, xlab="Overall Concordance Indices")


dev.copy2eps(file="forestplot-Overall.eps")
 











###################################################
### Kaplan Meier survival curves of all six datasets
###################################################
##pdf("KMPlot-all-datasets.pdf")

surv.data <- censor.time(surv.time=c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(vdx7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(nki7g)[ ,"t.dmfs"]) / 365, surv.event=c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(vdx7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(nki7g)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("vdx", nrow(pData(vdx7g))), rep("upp", nrow(pData(upp7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "vdx", "upp", "nki"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="", sub.title="", leg.pos="bottomright", leg.inset=0.05, o.text=FALSE, v.line=FALSE, h.line=FALSE, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkorange", "red", "darkblue", "darkgreen", "black", "brown"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE, ylim=c(0.1,1))

##dev.off()

dev.copy2eps(file="KMPlot-all-datasets.eps")

###################################################
### survival curves single genes
###################################################
## AURKA
##
##pdf("KMPlot-AURKA.pdf")

aurkaGs <- "aurka"
aurkaGid <- 6790
aurkaPaf <- "208079_s_at"
aurkaPagi <- "NM_003600"

surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(vdx7g)[ ,"t.dmfs"], pData(nki7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(vdx7g)[ ,"e.dmfs"], pData(nki7g)[ ,"e.dmfs"])            
aurka.exprs <- c(exprs(mainz7g)[aurkaPaf,], exprs(transbig7g)[aurkaPaf,], exprs(unt7g)[aurkaPaf,], exprs(upp7g)[aurkaPaf,], exprs(vdx7g)[aurkaPaf,], exprs(nki7g)[aurkaPagi,])
aurka.exprs.length <- c(length(exprs(mainz7g)[aurkaPaf,]), length(exprs(transbig7g)[aurkaPaf,]), length(exprs(unt7g)[aurkaPaf,]), length(exprs(upp7g)[aurkaPaf,]), length(exprs(vdx7g)[aurkaPaf,]), length(exprs(nki7g)[aurkaPagi,]))

pos <- 0
mygroup <- NULL
for(i in aurka.exprs.length){
  qq <- aurka.exprs[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.33, 0.66), na.rm=TRUE)
  qq[aurka.exprs[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[aurka.exprs[(pos+1):(pos+i)] >= myq[1] & aurka.exprs[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[aurka.exprs[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup <- c(mygroup,qq)
  pos <- pos + i
}

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("upp", nrow(pData(upp7g))), rep("vdx", nrow(pData(vdx7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Probability of survival for AURKA", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))

##dev.off()

dev.copy2eps(file="KMPlot-AURKA.eps")


###################################################
## VEGFA
##
##pdf("KMPlot-VEGFA.pdf")

aurkaGs <- "vegfa"
aurkaGid <- 7422
aurkaPaf <- "211527_x_at"
aurkaPagi <- "NM_003376"

surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(vdx7g)[ ,"t.dmfs"], pData(nki7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(vdx7g)[ ,"e.dmfs"], pData(nki7g)[ ,"e.dmfs"])            
aurka.exprs <- c(exprs(mainz7g)[aurkaPaf,], exprs(transbig7g)[aurkaPaf,], exprs(unt7g)[aurkaPaf,], exprs(upp7g)[aurkaPaf,], exprs(vdx7g)[aurkaPaf,], exprs(nki7g)[aurkaPagi,])
aurka.exprs.length <- c(length(exprs(mainz7g)[aurkaPaf,]), length(exprs(transbig7g)[aurkaPaf,]), length(exprs(unt7g)[aurkaPaf,]), length(exprs(upp7g)[aurkaPaf,]), length(exprs(vdx7g)[aurkaPaf,]), length(exprs(nki7g)[aurkaPagi,]))

pos <- 0
mygroup <- NULL
for(i in aurka.exprs.length){
  qq <- aurka.exprs[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.33, 0.66), na.rm=TRUE)
  qq[aurka.exprs[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[aurka.exprs[(pos+1):(pos+i)] >= myq[1] & aurka.exprs[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[aurka.exprs[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup <- c(mygroup,qq)
  pos <- pos + i
}

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("upp", nrow(pData(upp7g))), rep("vdx", nrow(pData(vdx7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Probability of survival for VEGFA", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))

##dev.off()

dev.copy2eps(file="KMPlot-VEGFA.eps")





###################################################
### ERBB2 KM curve
##pdf("KMPlot-ERBB2.pdf")

aurkaGs <- "erbb2"
aurkaGid <- 2064
aurkaPaf <- "216836_s_at"
aurkaPagi <- "NM_004448"

surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(vdx7g)[ ,"t.dmfs"], pData(nki7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(vdx7g)[ ,"e.dmfs"], pData(nki7g)[ ,"e.dmfs"])            
aurka.exprs <- c(exprs(mainz7g)[aurkaPaf,], exprs(transbig7g)[aurkaPaf,], exprs(unt7g)[aurkaPaf,], exprs(upp7g)[aurkaPaf,], exprs(vdx7g)[aurkaPaf,], exprs(nki7g)[aurkaPagi,])
aurka.exprs.length <- c(length(exprs(mainz7g)[aurkaPaf,]), length(exprs(transbig7g)[aurkaPaf,]), length(exprs(unt7g)[aurkaPaf,]), length(exprs(upp7g)[aurkaPaf,]), length(exprs(vdx7g)[aurkaPaf,]), length(exprs(nki7g)[aurkaPagi,]))

pos <- 0
mygroup <- NULL
for(i in aurka.exprs.length){
  qq <- aurka.exprs[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.33, 0.66), na.rm=TRUE)
  qq[aurka.exprs[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[aurka.exprs[(pos+1):(pos+i)] >= myq[1] & aurka.exprs[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[aurka.exprs[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup <- c(mygroup,qq)
  pos <- pos + i
}

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("upp", nrow(pData(upp7g))), rep("vdx", nrow(pData(vdx7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Probability of survival for ERBB2", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))

dev.copy2eps(file="KMPlot-ERBB2.eps")





###################################################
### ESR1 KM curve
##pdf("KMPlot-ESR1.pdf")

aurkaGs <- "esr1"
aurkaGid <- 2099
aurkaPaf <- "205225_at"
aurkaPagi <- "NM_000125"

surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(vdx7g)[ ,"t.dmfs"], pData(nki7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(vdx7g)[ ,"e.dmfs"], pData(nki7g)[ ,"e.dmfs"])            
aurka.exprs <- c(exprs(mainz7g)[aurkaPaf,], exprs(transbig7g)[aurkaPaf,], exprs(unt7g)[aurkaPaf,], exprs(upp7g)[aurkaPaf,], exprs(vdx7g)[aurkaPaf,], exprs(nki7g)[aurkaPagi,])
aurka.exprs.length <- c(length(exprs(mainz7g)[aurkaPaf,]), length(exprs(transbig7g)[aurkaPaf,]), length(exprs(unt7g)[aurkaPaf,]), length(exprs(upp7g)[aurkaPaf,]), length(exprs(vdx7g)[aurkaPaf,]), length(exprs(nki7g)[aurkaPagi,]))

pos <- 0
mygroup <- NULL
for(i in aurka.exprs.length){
  qq <- aurka.exprs[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.33, 0.66), na.rm=TRUE)
  qq[aurka.exprs[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[aurka.exprs[(pos+1):(pos+i)] >= myq[1] & aurka.exprs[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[aurka.exprs[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup <- c(mygroup,qq)
  pos <- pos + i
}

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("upp", nrow(pData(upp7g))), rep("vdx", nrow(pData(vdx7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Probability of survival for ESR1", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))

##dev.off()

dev.copy2eps(file="KMPlot-ESR1.eps")



###################################################
### PLAU KM curve
##pdf("KMPlot-PLAU.pdf")

aurkaGs <- "plau"
aurkaGid <- 5328
aurkaPaf <- "211668_s_at"
aurkaPagi <- "NM_002658"

surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(vdx7g)[ ,"t.dmfs"], pData(nki7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(vdx7g)[ ,"e.dmfs"], pData(nki7g)[ ,"e.dmfs"])            
aurka.exprs <- c(exprs(mainz7g)[aurkaPaf,], exprs(transbig7g)[aurkaPaf,], exprs(unt7g)[aurkaPaf,], exprs(upp7g)[aurkaPaf,], exprs(vdx7g)[aurkaPaf,], exprs(nki7g)[aurkaPagi,])
aurka.exprs.length <- c(length(exprs(mainz7g)[aurkaPaf,]), length(exprs(transbig7g)[aurkaPaf,]), length(exprs(unt7g)[aurkaPaf,]), length(exprs(upp7g)[aurkaPaf,]), length(exprs(vdx7g)[aurkaPaf,]), length(exprs(nki7g)[aurkaPagi,]))

pos <- 0
mygroup <- NULL
for(i in aurka.exprs.length){
  qq <- aurka.exprs[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.33, 0.66), na.rm=TRUE)
  qq[aurka.exprs[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[aurka.exprs[(pos+1):(pos+i)] >= myq[1] & aurka.exprs[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[aurka.exprs[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup <- c(mygroup,qq)
  pos <- pos + i
}

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("upp", nrow(pData(upp7g))), rep("vdx", nrow(pData(vdx7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Probability of survival for PLAU", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))

##dev.off()

dev.copy2eps(file="KMPlot-PLAU.eps")




###################################################
### STAT1 KM curve
##pdf("KMPlot-STAT1.pdf")

aurkaGs <- "stat1"
aurkaGid <- 6772
aurkaPaf <- "209969_s_at"
aurkaPagi <- "NM_007315"

surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(vdx7g)[ ,"t.dmfs"], pData(nki7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(vdx7g)[ ,"e.dmfs"], pData(nki7g)[ ,"e.dmfs"])            
aurka.exprs <- c(exprs(mainz7g)[aurkaPaf,], exprs(transbig7g)[aurkaPaf,], exprs(unt7g)[aurkaPaf,], exprs(upp7g)[aurkaPaf,], exprs(vdx7g)[aurkaPaf,], exprs(nki7g)[aurkaPagi,])
aurka.exprs.length <- c(length(exprs(mainz7g)[aurkaPaf,]), length(exprs(transbig7g)[aurkaPaf,]), length(exprs(unt7g)[aurkaPaf,]), length(exprs(upp7g)[aurkaPaf,]), length(exprs(vdx7g)[aurkaPaf,]), length(exprs(nki7g)[aurkaPagi,]))

pos <- 0
mygroup <- NULL
for(i in aurka.exprs.length){
  qq <- aurka.exprs[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.33, 0.66), na.rm=TRUE)
  qq[aurka.exprs[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[aurka.exprs[(pos+1):(pos+i)] >= myq[1] & aurka.exprs[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[aurka.exprs[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup <- c(mygroup,qq)
  pos <- pos + i
}

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("upp", nrow(pData(upp7g))), rep("vdx", nrow(pData(vdx7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Probability of survival for STAT1", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))

##dev.off()

dev.copy2eps(file="KMPlot-STAT1.eps")



###################################################
### CASP3 KM curve
##pdf("KMPlot-CASP3.pdf")

aurkaGs <- "casp3"
aurkaGid <- 836
aurkaPaf <- "202763_at"
aurkaPagi <- "NM_004346"

surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(vdx7g)[ ,"t.dmfs"], pData(nki7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(vdx7g)[ ,"e.dmfs"], pData(nki7g)[ ,"e.dmfs"])            
aurka.exprs <- c(exprs(mainz7g)[aurkaPaf,], exprs(transbig7g)[aurkaPaf,], exprs(unt7g)[aurkaPaf,], exprs(upp7g)[aurkaPaf,], exprs(vdx7g)[aurkaPaf,], exprs(nki7g)[aurkaPagi,])
aurka.exprs.length <- c(length(exprs(mainz7g)[aurkaPaf,]), length(exprs(transbig7g)[aurkaPaf,]), length(exprs(unt7g)[aurkaPaf,]), length(exprs(upp7g)[aurkaPaf,]), length(exprs(vdx7g)[aurkaPaf,]), length(exprs(nki7g)[aurkaPagi,]))

pos <- 0
mygroup <- NULL
for(i in aurka.exprs.length){
  qq <- aurka.exprs[(pos+1):(pos+i)]
  myq <- quantile(qq, probs=c(0.33, 0.66), na.rm=TRUE)
  qq[aurka.exprs[(pos+1):(pos+i)] < myq[1]] <- 1
  qq[aurka.exprs[(pos+1):(pos+i)] >= myq[1] & aurka.exprs[(pos+1):(pos+i)] < myq[2]] <- 2
  qq[aurka.exprs[(pos+1):(pos+i)] > myq[2]] <- 3
  qq <- factor(x=qq, levels=1:3)
  mygroup <- c(mygroup,qq)
  pos <- pos + i
}

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("upp", nrow(pData(upp7g))), rep("vdx", nrow(pData(vdx7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Probability of survival for CASP3", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))

##dev.off()

dev.copy2eps(file="KMPlot-CASP3.eps")






allData <- rbind(cindexall.mainz.small[,"cindex"], cindexall.transbig.small[,"cindex"], cindexall.upp.small[,"cindex"], cindexall.unt.small[,"cindex"], cindexall.vdx.small[,"cindex"], cindexall.nki.small[,"cindex"])
boxplot(allData)
lines(allData)


plot(cindexall.mainz.small[,"cindex"],type="p",col="red", axes=FALSE, xlim=c(0.5,7.5), ylim=c(0.3,0.7), xlab="Gene Symbol", ylab="Concordance Index", main="Concordance Indices from all Datasets plotted for each Gene", lwd=3)
lines(cindexall.transbig.small[,"cindex"],type="p",col="blue", lwd=3)
lines(cindexall.upp.small[,"cindex"],type="p",col="lightblue", lwd=3)
lines(cindexall.unt.small[,"cindex"],type="p",col="green", lwd=3)
lines(cindexall.vdx.small[,"cindex"],type="p",col="black", lwd=3)
lines(cindexall.nki.small[,"cindex"],type="p",col="yellow", lwd=3)
box()
axis(1, at=1:7, lab=c(tolower(fData(mainz7g)[,"Gene.symbol"])))
axis(2, las=1, at=0.1*3:7)
abline(v=1:length(cindexall.mainz.small), lty=3)
abline(h=c(0.45,0.55), lty=2, col="blue")
abline(h=0.5, lty=3, col="darkred")
legend("bottomright", c("MAINZ","TRANSBIG","UPP","UNT","VDX","NKI"), col=c("red","blue","lightblue","green","black","yellow"), pch=1, pt.lwd=3, bg="white")
boxplot(allData, add=TRUE, axes=FALSE, ylim=FALSE, xlim=FALSE, border="grey")





dev.copy2eps(file="hist-boxplot-cindex-all-datasets.eps")





####################################
### TDROCCs
####################################
tc <- 10 * 365

####################################
## AURKA MAINZ
#surv.data <- censor.time(surv.time=pData(mainz7g)[,"t.dmfs"] / 365, surv.event=pData(mainz7g)[,"e.dmfs"], time.cens=tc / 365)


surv.data.mainz <- censor.time(surv.time=c(pData(mainz7g)[ ,"t.dmfs"]) / 365, surv.event=c(pData(mainz7g)[ ,"e.dmfs"]), time.cens=tc / 365)
surv.data.transbig <- censor.time(surv.time=c(pData(transbig7g)[ ,"t.dmfs"]) / 365, surv.event=c(pData(transbig7g)[ ,"e.dmfs"]), time.cens=tc / 365)
surv.data.upp <- censor.time(surv.time=c(pData(upp7g)[ ,"t.rfs"]) / 365, surv.event=c(pData(upp7g)[ ,"e.rfs"]), time.cens=tc / 365)
surv.data.unt <- censor.time(surv.time=c(pData(unt7g)[ ,"t.dmfs"]) / 365, surv.event=c(pData(unt7g)[ ,"e.dmfs"]), time.cens=tc / 365)
surv.data.vdx <- censor.time(surv.time=c(pData(vdx7g)[ ,"t.dmfs"]) / 365, surv.event=c(pData(vdx7g)[ ,"e.dmfs"]), time.cens=tc / 365)
surv.data.nki <- censor.time(surv.time=c(pData(nki7g)[ ,"t.dmfs"]) / 365, surv.event=c(pData(nki7g)[ ,"e.dmfs"]), time.cens=tc / 365)


stime.mainz <- surv.data.mainz$surv.time.cens
sevent.mainz <- surv.data.mainz$surv.event.cens

stime.transbig <- surv.data.transbig$surv.time.cens
sevent.transbig <- surv.data.transbig$surv.event.cens

stime.upp <- surv.data.upp$surv.time.cens
sevent.upp <- surv.data.upp$surv.event.cens

stime.unt <- surv.data.unt$surv.time.cens
sevent.unt <- surv.data.unt$surv.event.cens

stime.vdx <- surv.data.vdx$surv.time.cens
sevent.vdx <- surv.data.vdx$surv.event.cens

stime.nki <- surv.data.nki$surv.time.cens
sevent.nki <- surv.data.nki$surv.event.cens


exprsData1 <- c(rescale(exprs(mainz7g)[3,] , q=0.05, na.rm=TRUE))
exprsData2 <- c(rescale(exprs(transbig7g)[3,] , q=0.05, na.rm=TRUE))
exprsData3 <- c(rescale(exprs(upp7g)[3,] , q=0.05, na.rm=TRUE))
exprsData4 <- c(rescale(exprs(unt7g)[3,] , q=0.05, na.rm=TRUE))
exprsData5 <- c(rescale(exprs(vdx7g)[3,] , q=0.05, na.rm=TRUE))
exprsData6 <- c(rescale(exprs(nki7g)[3,] , q=0.05, na.rm=TRUE))

timeTest <- unique(sort(stime[sevent == 1]))


################################
## 3 YEARS
tdroc1 <- tdrocc(x=exprsData1, surv.time=stime.mainz, surv.event=sevent.mainz, time=3, na.rm=TRUE, verbose=FALSE)
tdroc2 <- tdrocc(x=exprsData2, surv.time=stime.transbig, surv.event=sevent.transbig, time=3, na.rm=TRUE, verbose=FALSE)
tdroc3 <- tdrocc(x=exprsData3, surv.time=stime.upp, surv.event=sevent.upp, time=3, na.rm=TRUE, verbose=FALSE)
tdroc4 <- tdrocc(x=exprsData4, surv.time=stime.unt, surv.event=sevent.unt, time=3, na.rm=TRUE, verbose=FALSE)
tdroc5 <- tdrocc(x=exprsData5, surv.time=stime.vdx, surv.event=sevent.vdx, time=3, na.rm=TRUE, verbose=FALSE)
tdroc6 <- tdrocc(x=exprsData6, surv.time=stime.nki, surv.event=sevent.nki, time=3, na.rm=TRUE, verbose=FALSE)

##plot the time-dependent ROC curve
plot(x=1-tdroc1$spec, y=tdroc1$sens, type="l", xlab="1 - specificity", ylab="sensitivity", xlim=c(0, 1), ylim=c(0, 1), col="red")
lines(x=c(0,1), y=c(0,1), lty=3, col="red")

lines(x=1-tdroc2$spec, y=tdroc2$sens, type="l", col="royalblue")
lines(x=1-tdroc3$spec, y=tdroc3$sens, type="l", col="darkorange")
lines(x=1-tdroc4$spec, y=tdroc4$sens, type="l", col="green")
lines(x=1-tdroc5$spec, y=tdroc5$sens, type="l", col="black")
lines(x=1-tdroc6$spec, y=tdroc6$sens, type="l", col="yellow")
legend("bottomright", c("MAINZ","TRANSBIG","UPP","UNT","VDX","NKI"), col=c("red","royalblue","darkorange","green","black","yellow"), lty=1, bg="white")

dev.copy2eps(file="tdrocc-AURKA-3years.eps")



################################
## 5 YEARS
tdroc1 <- tdrocc(x=exprsData1, surv.time=stime.mainz, surv.event=sevent.mainz, time=5, na.rm=TRUE, verbose=FALSE)
tdroc2 <- tdrocc(x=exprsData2, surv.time=stime.transbig, surv.event=sevent.transbig, time=5, na.rm=TRUE, verbose=FALSE)
tdroc3 <- tdrocc(x=exprsData3, surv.time=stime.upp, surv.event=sevent.upp, time=5, na.rm=TRUE, verbose=FALSE)
tdroc4 <- tdrocc(x=exprsData4, surv.time=stime.unt, surv.event=sevent.unt, time=5, na.rm=TRUE, verbose=FALSE)
tdroc5 <- tdrocc(x=exprsData5, surv.time=stime.vdx, surv.event=sevent.vdx, time=5, na.rm=TRUE, verbose=FALSE)
tdroc6 <- tdrocc(x=exprsData6, surv.time=stime.nki, surv.event=sevent.nki, time=5, na.rm=TRUE, verbose=FALSE)

##plot the time-dependent ROC curve
plot(x=1-tdroc1$spec, y=tdroc1$sens, type="l", xlab="1 - specificity", ylab="sensitivity", xlim=c(0, 1), ylim=c(0, 1), col="red")
lines(x=c(0,1), y=c(0,1), lty=3, col="red")

lines(x=1-tdroc2$spec, y=tdroc2$sens, type="l", col="royalblue")
lines(x=1-tdroc3$spec, y=tdroc3$sens, type="l", col="darkorange")
lines(x=1-tdroc4$spec, y=tdroc4$sens, type="l", col="green")
lines(x=1-tdroc5$spec, y=tdroc5$sens, type="l", col="black")
lines(x=1-tdroc6$spec, y=tdroc6$sens, type="l", col="yellow")
legend("bottomright", c("MAINZ","TRANSBIG","UPP","UNT","VDX","NKI"), col=c("red","royalblue","darkorange","green","black","yellow"), lty=1, bg="white")

dev.copy2eps(file="tdrocc-AURKA-6years.eps")



################################
## 10 YEARS
tdroc1 <- tdrocc(x=exprsData1, surv.time=stime.mainz, surv.event=sevent.mainz, time=10, na.rm=TRUE, verbose=FALSE)
tdroc2 <- tdrocc(x=exprsData2, surv.time=stime.transbig, surv.event=sevent.transbig, time=10, na.rm=TRUE, verbose=FALSE)
tdroc3 <- tdrocc(x=exprsData3, surv.time=stime.upp, surv.event=sevent.upp, time=10, na.rm=TRUE, verbose=FALSE)
tdroc4 <- tdrocc(x=exprsData4, surv.time=stime.unt, surv.event=sevent.unt, time=10, na.rm=TRUE, verbose=FALSE)
tdroc5 <- tdrocc(x=exprsData5, surv.time=stime.vdx, surv.event=sevent.vdx, time=10, na.rm=TRUE, verbose=FALSE)
tdroc6 <- tdrocc(x=exprsData6, surv.time=stime.nki, surv.event=sevent.nki, time=10, na.rm=TRUE, verbose=FALSE)

##plot the time-dependent ROC curve
plot(x=1-tdroc1$spec, y=tdroc1$sens, type="l", xlab="1 - specificity", ylab="sensitivity", xlim=c(0, 1), ylim=c(0, 1), col="red")
lines(x=c(0,1), y=c(0,1), lty=3, col="red")

lines(x=1-tdroc2$spec, y=tdroc2$sens, type="l", col="royalblue")
lines(x=1-tdroc3$spec, y=tdroc3$sens, type="l", col="darkorange")
lines(x=1-tdroc4$spec, y=tdroc4$sens, type="l", col="green")
lines(x=1-tdroc5$spec, y=tdroc5$sens, type="l", col="black")
lines(x=1-tdroc6$spec, y=tdroc6$sens, type="l", col="yellow")
legend("bottomright", c("MAINZ","TRANSBIG","UPP","UNT","VDX","NKI"), col=c("red","royalblue","darkorange","green","black","yellow"), lty=1, bg="white")

dev.copy2eps(file="tdrocc-AURKA-10years.eps")



####################################
## VEGFA MAINZ
surv.data <- censor.time(surv.time=pData(mainz7g)[,"t.dmfs"] / 365, surv.event=pData(mainz7g)[,"e.dmfs"], time.cens=tc / 365)
stime <- surv.data$surv.time.cens
sevent <- surv.data$surv.event.cens

tdroc <- tdrocc(x=rescale(exprs(mainz7g)[3,] , q=0.05, na.rm=TRUE), surv.time=stime, surv.event=sevent, time=1, na.rm=TRUE, verbose=FALSE)

##plot the time-dependent ROC curve
plot(x=1-tdroc$spec, y=tdroc$sens, type="l", xlab="1 - specificity", ylab="sensitivity", xlim=c(0, 1), ylim=c(0, 1))
lines(x=c(0,1), y=c(0,1), lty=3, col="red")

##dev.off()







####################################
## Brier score wrt the time

#library(survcomp)
library(Biobase)
tc <- 10 * 365  

gsList <- c("ESR1", "ERBB2", "AURKA", "PLAU", "VEGF", "STAT1", "CASP3")
probesAffy <- c("205225_at", "216836_s_at", "208079_s_at", "211668_s_at", "211527_x_at", "209969_s_at", "202763_at")



#surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"])             
#surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"])
#mygene <- t(cbind(exprs(mainz7g), exprs(transbig7g)))


surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"])
mygene <- t(cbind(exprs(mainz7g)))

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd.ts <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], mygene)
dd.ts <- dd.ts[complete.cases(dd.ts), ,drop=FALSE]

sbrier.score2proba.ts.all <- NULL

## KM model
sf.tr <- survfit(Surv(time, event) ~ 1, data=dd.ts)
lsf.tr <- rep(list(sf.tr), nrow(dd.ts))
names(lsf.tr) <- dimnames(dd.ts)[[1]]

lsf.trts <- rep(list(sf.tr), nrow(dd.ts))
names(lsf.trts) <- dimnames(dd.ts)[[1]]
tt <- sort(dd.ts$time)
btime <- tt[tt >= 0 & tt <= max(tt, na.rm=TRUE)]
utime <- unique(sort(dd.ts$time[dd.ts$event == 1]))
bsc <- rep(NA, length(btime))
for(i in 1:length(utime)) {
	bsc[is.na(bsc) & btime <= utime[i]] <- sbrier(obj=Surv(dd.ts$time, dd.ts$event), pred=lsf.trts, btime=utime[i])
}
bsc[is.na(bsc)] <- bsc[ min(which(is.na(bsc)))-1] 
diffs <- c(btime[1], btime[2:length(btime)] - btime[1:(length(btime) - 1)])
bsc.int <- sum(diffs * bsc)/max(btime)
res.ts <- list("time"=btime, "bsc"=bsc, "bsc.integrated"=bsc.int)

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "KM"

## for ESR1
ddtt <- cbind(dd.ts[ , c("time", "event", "X205225_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "ESR1"

## for ERBB2
ddtt <- cbind(dd.ts[ , c("time", "event", "X216836_s_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "ERBB2"

## for AURKA
ddtt <- cbind(dd.ts[ , c("time", "event", "X208079_s_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "AURKA"

## for PLAU
ddtt <- cbind(dd.ts[ , c("time", "event", "X211668_s_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "PLAU"

## for VEGF
ddtt <- cbind(dd.ts[ , c("time", "event", "X211527_x_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "VEGF"

## for STAT1
ddtt <- cbind(dd.ts[ , c("time", "event", "X209969_s_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "STAT1"

## for CASP3
ddtt <- cbind(dd.ts[ , c("time", "event", "X202763_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "CASP3"

myperf <-  sbrier.score2proba.ts.all

mycol <- c("#000000FF", rainbow(length(myperf))[-1])
mylty <- 1:(length(myperf))
mylwd <- c(3, rep(2, length(myperf)-1))
for(i in 1:length(myperf)) {
	if(i == 1) {
		plot(myperf[[i]]$time, myperf[[i]]$bsc, typ="l", xlab="Time (years)", ylab="Brier score", col=mycol[i], lty=mylty[i], ylim=c(0,0.22), lwd=mylwd[i], main="Brier Scores for MAINZ")
	} else {  lines(myperf[[i]]$time, myperf[[i]]$bsc, col=mycol[i], lty=mylty[i], lwd=mylwd[i]) }
	mysbrier.int <- unlist(lapply(myperf, function(x) { return(x$bsc.integrated)}))
}
smartlegend(x="left", y="top", legend=sprintf("%s, IBSC = %.3g", names(myperf), mysbrier.int), col=mycol, lty=mylty, lwd=mylwd)

dev.copy2eps(file="BSC-mainz.eps")




####################################################
## for TRANSBIG
####################################################

surv.time.all <- c(pData(transbig7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(transbig7g)[ ,"e.dmfs"])
mygene <- t(cbind(exprs(transbig7g)))

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd.ts <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], mygene)
dd.ts <- dd.ts[complete.cases(dd.ts), ,drop=FALSE]

sbrier.score2proba.ts.all <- NULL

## KM model
sf.tr <- survfit(Surv(time, event) ~ 1, data=dd.ts)
lsf.tr <- rep(list(sf.tr), nrow(dd.ts))
names(lsf.tr) <- dimnames(dd.ts)[[1]]

lsf.trts <- rep(list(sf.tr), nrow(dd.ts))
names(lsf.trts) <- dimnames(dd.ts)[[1]]
tt <- sort(dd.ts$time)
btime <- tt[tt >= 0 & tt <= max(tt, na.rm=TRUE)]
utime <- unique(sort(dd.ts$time[dd.ts$event == 1]))
bsc <- rep(NA, length(btime))
for(i in 1:length(utime)) {
	bsc[is.na(bsc) & btime <= utime[i]] <- sbrier(obj=Surv(dd.ts$time, dd.ts$event), pred=lsf.trts, btime=utime[i])
}
bsc[is.na(bsc)] <- bsc[ min(which(is.na(bsc)))-1] 
diffs <- c(btime[1], btime[2:length(btime)] - btime[1:(length(btime) - 1)])
bsc.int <- sum(diffs * bsc)/max(btime)
res.ts <- list("time"=btime, "bsc"=bsc, "bsc.integrated"=bsc.int)

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "KM"

## for ESR1
ddtt <- cbind(dd.ts[ , c("time", "event", "X205225_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "ESR1"

## for ERBB2
ddtt <- cbind(dd.ts[ , c("time", "event", "X216836_s_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "ERBB2"

## for AURKA
ddtt <- cbind(dd.ts[ , c("time", "event", "X208079_s_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "AURKA"

## for PLAU
ddtt <- cbind(dd.ts[ , c("time", "event", "X211668_s_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "PLAU"

## for VEGF
ddtt <- cbind(dd.ts[ , c("time", "event", "X211527_x_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "VEGF"

## for STAT1
ddtt <- cbind(dd.ts[ , c("time", "event", "X209969_s_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "STAT1"

## for CASP3
ddtt <- cbind(dd.ts[ , c("time", "event", "X202763_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "CASP3"

myperf <-  sbrier.score2proba.ts.all

mycol <- c("#000000FF", rainbow(length(myperf))[-1])
mylty <- 1:(length(myperf))
mylwd <- c(3, rep(2, length(myperf)-1))
for(i in 1:length(myperf)) {
	if(i == 1) {
		plot(myperf[[i]]$time, myperf[[i]]$bsc, typ="l", xlab="Time (years)", ylab="Brier score", col=mycol[i], lty=mylty[i], ylim=c(0,0.3), lwd=mylwd[i], main="Brier Scores for TRANSBIG")
	} else {  lines(myperf[[i]]$time, myperf[[i]]$bsc, col=mycol[i], lty=mylty[i], lwd=mylwd[i]) }
	mysbrier.int <- unlist(lapply(myperf, function(x) { return(x$bsc.integrated)}))
}
smartlegend(x="left", y="top", legend=sprintf("%s, IBSC = %.3g", names(myperf), mysbrier.int), col=mycol, lty=mylty, lwd=mylwd)

dev.copy2eps(file="BSC-transbig.eps")



####################################################
## for UPP
####################################################

surv.time.all <- c(pData(upp7g)[ ,"t.rfs"])             
surv.event.all <- c(pData(upp7g)[ ,"e.rfs"])
mygene <- t(cbind(exprs(upp7g)))

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd.ts <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], mygene)
dd.ts <- dd.ts[complete.cases(dd.ts), ,drop=FALSE]

sbrier.score2proba.ts.all <- NULL

## KM model
sf.tr <- survfit(Surv(time, event) ~ 1, data=dd.ts)
lsf.tr <- rep(list(sf.tr), nrow(dd.ts))
names(lsf.tr) <- dimnames(dd.ts)[[1]]

lsf.trts <- rep(list(sf.tr), nrow(dd.ts))
names(lsf.trts) <- dimnames(dd.ts)[[1]]
tt <- sort(dd.ts$time)
btime <- tt[tt >= 0 & tt <= max(tt, na.rm=TRUE)]
utime <- unique(sort(dd.ts$time[dd.ts$event == 1]))
bsc <- rep(NA, length(btime))
for(i in 1:length(utime)) {
	bsc[is.na(bsc) & btime <= utime[i]] <- sbrier(obj=Surv(dd.ts$time, dd.ts$event), pred=lsf.trts, btime=utime[i])
}
bsc[is.na(bsc)] <- bsc[ min(which(is.na(bsc)))-1] 
diffs <- c(btime[1], btime[2:length(btime)] - btime[1:(length(btime) - 1)])
bsc.int <- sum(diffs * bsc)/max(btime)
res.ts <- list("time"=btime, "bsc"=bsc, "bsc.integrated"=bsc.int)

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "KM"

## for ESR1
ddtt <- cbind(dd.ts[ , c("time", "event", "X205225_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "ESR1"

## for ERBB2
ddtt <- cbind(dd.ts[ , c("time", "event", "X216836_s_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "ERBB2"

## for AURKA
ddtt <- cbind(dd.ts[ , c("time", "event", "X208079_s_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "AURKA"

## for PLAU
ddtt <- cbind(dd.ts[ , c("time", "event", "X211668_s_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "PLAU"

## for VEGF
ddtt <- cbind(dd.ts[ , c("time", "event", "X211527_x_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "VEGF"

## for STAT1
ddtt <- cbind(dd.ts[ , c("time", "event", "X209969_s_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "STAT1"

## for CASP3
ddtt <- cbind(dd.ts[ , c("time", "event", "X202763_at")])
colnames(ddtt) <- c("time", "event", "score")
res.ts <- sbrier.score2proba(data.tr=ddtt, data.ts=ddtt, method="cox")

sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- "CASP3"

myperf <-  sbrier.score2proba.ts.all

mycol <- c("#000000FF", rainbow(length(myperf))[-1])
mylty <- 1:(length(myperf))
mylwd <- c(3, rep(2, length(myperf)-1))
for(i in 1:length(myperf)) {
	if(i == 1) {
		plot(myperf[[i]]$time, myperf[[i]]$bsc, typ="l", xlab="Time (years)", ylab="Brier score", col=mycol[i], lty=mylty[i], ylim=c(0,0.2), lwd=mylwd[i], main="Brier Score for UPP")
	} else {  lines(myperf[[i]]$time, myperf[[i]]$bsc, col=mycol[i], lty=mylty[i], lwd=mylwd[i]) }
	mysbrier.int <- unlist(lapply(myperf, function(x) { return(x$bsc.integrated)}))
}
smartlegend(x="left", y="top", legend=sprintf("%s, IBSC = %.3g", names(myperf), mysbrier.int), col=mycol, lty=mylty, lwd=mylwd)

dev.copy2eps(file="BSC-upp.eps")













#######################################
### VENN DIAGRAM
#######################################

#######################
## top low/high risk genes MAINZ
#######################
qm <- quantile(cindexall.mainz[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE)
##cindexall.mainz[cindexall.mainz[,"cindex"] > qm[[1]],]
low.risk.probes.mainz <- as.data.frame(cindexall.mainz[cindexall.mainz[,"cindex"] > qm[[2]],])
low.risk.probes.mainz <- low.risk.probes.mainz[order(low.risk.probes.mainz$cindex, decreasing=TRUE),]
high.risk.probes.mainz <- as.data.frame(cindexall.mainz[cindexall.mainz[,"cindex"] < qm[[1]],])
high.risk.probes.mainz <- high.risk.probes.mainz[order(high.risk.probes.mainz$cindex),]

low.risk.probes.mainz$gs <- fData(mainz)[fData(mainz)[rownames(low.risk.probes.mainz),"probe"],"Gene.symbol"]
high.risk.probes.mainz$gs <- fData(mainz)[fData(mainz)[rownames(high.risk.probes.mainz),"probe"],"Gene.symbol"]

low.risk.probes.mainz$gid <- fData(mainz)[fData(mainz)[rownames(low.risk.probes.mainz),"probe"],"EntrezGene.ID"]
high.risk.probes.mainz$gid <- fData(mainz)[fData(mainz)[rownames(high.risk.probes.mainz),"probe"],"EntrezGene.ID"]


#######################
## top low/high risk genes TRANSBIG
#######################
qm <- quantile(cindexall.transbig[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE)
##cindexall.transbig[cindexall.transbig[,"cindex"] > qm[[1]],]
low.risk.probes.transbig <- as.data.frame(cindexall.transbig[cindexall.transbig[,"cindex"] > qm[[2]],])
low.risk.probes.transbig <- low.risk.probes.transbig[order(low.risk.probes.transbig$cindex, decreasing=TRUE),]
high.risk.probes.transbig <- as.data.frame(cindexall.transbig[cindexall.transbig[,"cindex"] < qm[[1]],])
high.risk.probes.transbig <- high.risk.probes.transbig[order(high.risk.probes.transbig$cindex),]

low.risk.probes.transbig$gs <- fData(transbig)[fData(transbig)[rownames(low.risk.probes.transbig),"probe"],"Gene.symbol"]
high.risk.probes.transbig$gs <- fData(transbig)[fData(transbig)[rownames(high.risk.probes.transbig),"probe"],"Gene.symbol"]

low.risk.probes.transbig$gid <- fData(transbig)[fData(transbig)[rownames(low.risk.probes.transbig),"probe"],"EntrezGene.ID"]
high.risk.probes.transbig$gid <- fData(transbig)[fData(transbig)[rownames(high.risk.probes.transbig),"probe"],"EntrezGene.ID"]

#######################
## top low/high risk genes UPP
#######################
qm <- quantile(cindexall.upp[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE)
##cindexall.upp[cindexall.upp[,"cindex"] > qm[[1]],]
low.risk.probes.upp <- as.data.frame(cindexall.upp[cindexall.upp[,"cindex"] > qm[[2]],])
low.risk.probes.upp <- low.risk.probes.upp[order(low.risk.probes.upp$cindex, decreasing=TRUE),]
high.risk.probes.upp <- as.data.frame(cindexall.upp[cindexall.upp[,"cindex"] < qm[[1]],])
high.risk.probes.upp <- high.risk.probes.upp[order(high.risk.probes.upp$cindex),]

low.risk.probes.upp$gs <- fData(upp)[fData(upp)[rownames(low.risk.probes.upp),"probe"],"Gene.symbol"]
high.risk.probes.upp$gs <- fData(upp)[fData(upp)[rownames(high.risk.probes.upp),"probe"],"Gene.symbol"]

low.risk.probes.upp$gid <- fData(upp)[fData(upp)[rownames(low.risk.probes.upp),"probe"],"EntrezGene.ID"]
high.risk.probes.upp$gid <- fData(upp)[fData(upp)[rownames(high.risk.probes.upp),"probe"],"EntrezGene.ID"]

#######################
## top low/high risk genes UNT
#######################
qm <- quantile(cindexall.unt[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE)
##cindexall.unt[cindexall.unt[,"cindex"] > qm[[1]],]
low.risk.probes.unt <- as.data.frame(cindexall.unt[cindexall.unt[,"cindex"] > qm[[2]],])
low.risk.probes.unt <- low.risk.probes.unt[order(low.risk.probes.unt$cindex, decreasing=TRUE),]
high.risk.probes.unt <- as.data.frame(cindexall.unt[cindexall.unt[,"cindex"] < qm[[1]],])
high.risk.probes.unt <- high.risk.probes.unt[order(high.risk.probes.unt$cindex),]
low.risk.probes.unt$gs <- fData(unt)[fData(unt)[rownames(low.risk.probes.unt),"probe"],"Gene.symbol"]
high.risk.probes.unt$gs <- fData(unt)[fData(unt)[rownames(high.risk.probes.unt),"probe"],"Gene.symbol"]

low.risk.probes.unt$gid <- fData(unt)[fData(unt)[rownames(low.risk.probes.unt),"probe"],"EntrezGene.ID"]
high.risk.probes.unt$gid <- fData(unt)[fData(unt)[rownames(high.risk.probes.unt),"probe"],"EntrezGene.ID"]


#######################
## top low/high risk genes VDX
#######################
qm <- quantile(cindexall.vdx[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE)
##cindexall.vdx[cindexall.vdx[,"cindex"] > qm[[1]],]
low.risk.probes.vdx <- as.data.frame(cindexall.vdx[cindexall.vdx[,"cindex"] > qm[[2]],])
low.risk.probes.vdx <- low.risk.probes.vdx[order(low.risk.probes.vdx$cindex, decreasing=TRUE),]
high.risk.probes.vdx <- as.data.frame(cindexall.vdx[cindexall.vdx[,"cindex"] < qm[[1]],])
high.risk.probes.vdx <- high.risk.probes.vdx[order(high.risk.probes.vdx$cindex),]
low.risk.probes.vdx$gs <- fData(vdx)[fData(vdx)[rownames(low.risk.probes.vdx),"probe"],"Gene.symbol"]
high.risk.probes.vdx$gs <- fData(vdx)[fData(vdx)[rownames(high.risk.probes.vdx),"probe"],"Gene.symbol"]

low.risk.probes.vdx$gid <- fData(vdx)[fData(vdx)[rownames(low.risk.probes.vdx),"probe"],"EntrezGene.ID"]
high.risk.probes.vdx$gid <- fData(vdx)[fData(vdx)[rownames(high.risk.probes.vdx),"probe"],"EntrezGene.ID"]

#######################
## top low/high risk genes NKI
#######################
qm <- quantile(cindexall.nki[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE)
##cindexall.nki[cindexall.nki[,"cindex"] > qm[[1]],]
low.risk.probes.nki <- as.data.frame(cindexall.nki[cindexall.nki[,"cindex"] > qm[[2]],])
low.risk.probes.nki <- low.risk.probes.nki[order(low.risk.probes.nki$cindex, decreasing=TRUE),]
high.risk.probes.nki <- as.data.frame(cindexall.nki[cindexall.nki[,"cindex"] < qm[[1]],])
high.risk.probes.nki <- high.risk.probes.nki[order(high.risk.probes.nki$cindex),]

low.risk.probes.nki$gs <- fData(nki)[fData(nki)[rownames(low.risk.probes.nki),"probe"],"NCBI.gene.symbol"]
high.risk.probes.nki$gs <- fData(nki)[fData(nki)[rownames(high.risk.probes.nki),"probe"],"NCBI.gene.symbol"]

low.risk.probes.nki$gid <- fData(nki)[fData(nki)[rownames(low.risk.probes.nki),"probe"],"EntrezGene.ID"]
high.risk.probes.nki$gid <- fData(nki)[fData(nki)[rownames(high.risk.probes.nki),"probe"],"EntrezGene.ID"]

low.risk.probes.nki$gs2 <- fData(nki)[fData(nki)[rownames(low.risk.probes.nki),"probe"],"HUGO.gene.symbol"]
high.risk.probes.nki$gs2 <- fData(nki)[fData(nki)[rownames(high.risk.probes.nki),"probe"],"HUGO.gene.symbol"]


#######################
## VENN DIAGRAM sample data
#######################
gLowMainz <- as.numeric(unique(na.omit(low.risk.probes.mainz[1:1000,"gid"])))
gLowTransbig <- as.numeric(unique(na.omit(low.risk.probes.transbig[1:1000,"gid"])))
gLowUpp <- as.numeric(unique(na.omit(low.risk.probes.upp[1:1000,"gid"])))
gLowUnt <- as.numeric(unique(na.omit(low.risk.probes.unt[1:1000,"gid"])))
gLowVdx <- as.numeric(unique(na.omit(low.risk.probes.vdx[1:1000,"gid"])))
gLowNki <- as.numeric(unique(na.omit(low.risk.probes.nki[1:1000,"gid"])))

gHighMainz <- as.numeric(unique(na.omit(high.risk.probes.mainz[1:1000,"gid"])))
gHighTransbig <- as.numeric(unique(na.omit(high.risk.probes.transbig[1:1000,"gid"])))
gHighUpp <- as.numeric(unique(na.omit(high.risk.probes.upp[1:1000,"gid"])))
gHighUnt <- as.numeric(unique(na.omit(high.risk.probes.unt[1:1000,"gid"])))
gHighVdx <- as.numeric(unique(na.omit(high.risk.probes.vdx[1:1000,"gid"])))
gHighNki <- as.numeric(unique(na.omit(high.risk.probes.nki[1:1000,"gid"])))


g <- cbind(unique(gLowMainz[1:500]),gLowTransbig[1:500],gLowUpp[1:500],gLowUnt[1:500],unique(gHighMainz[1:500]),gHighTransbig[1:500],gHighUpp[1:500],gHighUnt[1:500])

colnames(g) <- c("MainzLow","TransbigLow", "UppLow", "UntLow", "MainzHigh", "TransbigHigh", "UppHigh", "UntHigh")

write.table(g, file="vennData.csv", sep=" ")

## get the intersect of all datasets
geneIdAll <- intersect(gHighMainz[1:500],intersect(gHighTransbig[1:500],intersect(gHighUpp[1:500],intersect(gHighUnt[1:500],intersect(gHighVdx[1:500],gHighNki[1:500])))))
geneIdMinusNki <- intersect(gHighMainz[1:500],intersect(gHighTransbig[1:500],intersect(gHighUpp[1:500],intersect(gHighUnt[1:500],gHighNki[1:500]))))
geneIdMinusVdx <- intersect(gHighMainz[1:500],intersect(gHighTransbig[1:500],intersect(gHighUpp[1:500],intersect(gHighUnt[1:500],gHighVdx[1:500]))))
geneIdMinusMainz <- intersect(gHighNki[1:500],intersect(gHighTransbig[1:500],intersect(gHighUpp[1:500],intersect(gHighUnt[1:500],gHighVdx[1:500]))))
geneIdMinusTransbig <- intersect(gHighMainz[1:500],intersect(gHighNki[1:500],intersect(gHighUpp[1:500],intersect(gHighUnt[1:500],gHighVdx[1:500]))))
geneIdMinusUpp <- intersect(gHighMainz[1:500],intersect(gHighTransbig[1:500],intersect(gHighNki[1:500],intersect(gHighUnt[1:500],gHighVdx[1:500]))))
geneIdMinusUnt <- intersect(gHighMainz[1:500],intersect(gHighTransbig[1:500],intersect(gHighUpp[1:500],intersect(gHighNki[1:500],gHighVdx[1:500]))))




getAnnot <- function(x) {
  tt <- fData(mainz)[ ,colnames(fData(mainz))=='EntrezGene.ID']==x
  tt[is.na(tt)]<-F
  fData(mainz)[tt,c("EntrezGene.ID", "probe", "Gene.symbol", "Gene.title", "GO.Function", "GO.Process", "GO.Component")]
}

rr <- NULL
> for(i in geneID){
+ rr <- rbind(rr,getAnnot(i))
+ }




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






#######################################
### Heatmap for Cindices
#######################################

cindex.heatmap.small <- cbind("mainz"=cindexall.mainz.small[,1], "transbig"=cindexall.transbig.small[,1], "upp"=cindexall.upp.small[,1], "unt"=cindexall.unt.small[,1], "vdx"=cindexall.vdx.small[,1], "nki"=cindexall.nki.small[,1], "overall"=ccindex.small[,1])

rownames(cindex.heatmap.small) <- gsList

heatmap(1-cindex.heatmap.small, Rowv=NA, Colv=NA, col=cm.colors(8), scale="none", margins=c(5,10), main="Concordance Indices as Heatmap, >0.5 = red")
dev.copy2eps(file="heatmap.cindex.small-low.eps")

heatmap.2(cindex.heatmap.small, Rowv=NA, Colv=NA, col=heat.colors(256), scale="none", margins=c(5,10), main="Concordance Indices as Heatmap, <0.5 = red")
dev.copy2eps(file="heatmap.cindex.small-high.eps") 

heatplot(cindex.heatmap.small, dend="none", scale="none", lowcol="blue", highcol="red", zlim=c(0,1))








#######################################
### subtype identification
#######################################
scmgene.mainz <- subtype.cluster(module.ESR1=scmgene$mod$ESR1, module.ERBB2=scmgene$mod$ERBB2, module.AURKA=scmgene$mod$AURKA, data=t(exprs(mainz)), annot=fData(mainz), do.mapping=FALSE, do.scale=TRUE, plot=TRUE, verbose=TRUE)
#table(scmgene.mainz$subtype2)
scmgene.mainz.p <- subtype.cluster.predict(sbt.model=scmgene.mainz$model, data=t(exprs(mainz)), annot=fData(mainz), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)



scmgene.mainz <- subtype.cluster.predict(sbt.model=scmgene, data=t(exprs(mainz)), annot=fData(mainz), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
scmgene.transbig <- subtype.cluster.predict(sbt.model=scmgene, data=t(exprs(transbig)), annot=fData(transbig), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
scmgene.upp <- subtype.cluster.predict(sbt.model=scmgene, data=t(exprs(upp)), annot=fData(upp), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
scmgene.unt <- subtype.cluster.predict(sbt.model=scmgene, data=t(exprs(unt)), annot=fData(unt), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
scmgene.vdx <- subtype.cluster.predict(sbt.model=scmgene, data=t(exprs(vdx)), annot=fData(vdx), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)
scmgene.nki <- subtype.cluster.predict(sbt.model=scmgene, data=t(exprs(nki)), annot=fData(nki), do.mapping=TRUE, do.scale=TRUE, do.prediction.strength=FALSE, do.BIC=FALSE, plot=TRUE, verbose=TRUE)


table(scmgene.mainz$subtype2)

###################################################
### survival curves AURKA using subtypes, all datasets
###################################################

aurkaGs <- "aurka"
aurkaGid <- 6790
aurkaPaf <- "208079_s_at"
aurkaPagi <- "NM_003600"

surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(vdx7g)[ ,"t.dmfs"], pData(nki7g)[ ,"t.dmfs"])  
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(vdx7g)[ ,"e.dmfs"], pData(nki7g)[ ,"e.dmfs"]) 
aurka.exprs <- c(exprs(mainz7g)[aurkaPaf,], exprs(transbig7g)[aurkaPaf,], exprs(unt7g)[aurkaPaf,], exprs(upp7g)[aurkaPaf,], exprs(vdx7g)[aurkaPaf,], exprs(nki7g)[aurkaPagi,])

aurka.exprs.group <- c(as.character(scmgene.mainz$subtype2), as.character(scmgene.transbig$subtype2), as.character(scmgene.upp$subtype2), as.character(scmgene.unt$subtype2), as.character(scmgene.vdx$subtype2), as.character(scmgene.nki$subtype2))

mygroup <- as.factor(aurka.exprs.group)

surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=as.numeric(mygroup))
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("upp", nrow(pData(upp7g))), rep("vdx", nrow(pData(vdx7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))

km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Prob. of survival for subgroups of patients according to subtypes", sub.title="", leg.text=c(levels(mygroup)), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred", "orange"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))



dev.copy2eps(file="KMPlot-AURKA-subtypes.eps")



###################################################
### survival curves STAT1, quantiles, all datasets
###################################################
pdf("KMPlot-STAT1-quantiles.pdf")

stat1Gs <- "stat1"
stat1Gid <- 6772
stat1Paf <- "209969_s_at"
stat1Pagi <- "NM_007315"

surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(vdx7g)[ ,"t.dmfs"], pData(nki7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(vdx7g)[ ,"e.dmfs"], pData(nki7g)[ ,"e.dmfs"])            
stat1.exprs <- c(exprs(mainz7g)[stat1Paf,], exprs(transbig7g)[stat1Paf,], exprs(unt7g)[stat1Paf,], exprs(upp7g)[stat1Paf,], exprs(vdx7g)[stat1Paf,], exprs(nki7g)[stat1Pagi,])

stat1.exprs.length <- c(length(exprs(mainz7g)[stat1Paf,]), length(exprs(transbig7g)[stat1Paf,]), length(exprs(unt7g)[stat1Paf,]), length(exprs(upp7g)[stat1Paf,]), length(exprs(vdx7g)[stat1Paf,]), length(exprs(nki7g)[stat1Pagi,]))

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
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("upp", nrow(pData(upp7g))), rep("vdx", nrow(pData(vdx7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Probability of survival for STAT1", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))

dev.off()

dev.copy2eps(file="KMPlot-STAT1-subtypes.eps")




###################################################
### survival curves STAT1 using subtypes, all datasets
###################################################
pdf("KMPlot-STAT1-subtype-pmh.pdf")

stat1Gs <- "stat1"
stat1Gid <- 6772
stat1Paf <- "209969_s_at"
stat1Pagi <- "NM_007315"

stat1.exprs.group <- c(as.character(scmgene.mainz$subtype2), as.character(scmgene.transbig$subtype2), as.character(scmgene.upp$subtype2), as.character(scmgene.unt$subtype2), as.character(scmgene.vdx$subtype2), as.character(scmgene.nki$subtype2))
mm <- stat1.exprs.group=="ER-/HER2-"
mm[is.na(mm)] <- F

pmh <- stat1.exprs.group=="ER+/HER2- High Prolif"
pmh[is.na(pmh)] <- F

pml <- stat1.exprs.group=="ER+/HER2- Low Prolif"
pml[is.na(pml)] <- F

her2p <- stat1.exprs.group=="HER2+"
her2p[is.na(her2p)] <- F

surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(vdx7g)[ ,"t.dmfs"], pData(nki7g)[ ,"t.dmfs"])
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(vdx7g)[ ,"e.dmfs"], pData(nki7g)[ ,"e.dmfs"]) 
stat1.exprs <- c(exprs(mainz7g)[stat1Paf,], exprs(transbig7g)[stat1Paf,], exprs(unt7g)[stat1Paf,], exprs(upp7g)[stat1Paf,], exprs(vdx7g)[stat1Paf,], exprs(nki7g)[stat1Pagi,])

stat1.exprs.length <- c(length(exprs(mainz7g)[stat1Paf,]), length(exprs(transbig7g)[stat1Paf,]), length(exprs(unt7g)[stat1Paf,]), length(exprs(upp7g)[stat1Paf,]), length(exprs(vdx7g)[stat1Paf,]), length(exprs(nki7g)[stat1Pagi,]))

#mygroup <- as.factor(stat1.exprs.group)                      

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
dd <- dd[her2p,]

km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Probability of survival for STAT1, HER2+", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))

dev.off()

dev.copy2eps(file="KMPlot-STAT1-subtypes.eps")



###################################################
### survival curves STAT1 using subtypes, all datasets
###################################################
pdf("KMPlot-STAT1-subtypes.pdf")

stat1Gs <- "stat1"
stat1Gid <- 6772
stat1Paf <- "209969_s_at"
stat1Pagi <- "NM_007315"

surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"], pData(unt7g)[ ,"t.dmfs"], pData(upp7g)[ ,"t.rfs"], pData(vdx7g)[ ,"t.dmfs"], pData(nki7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"], pData(unt7g)[ ,"e.dmfs"], pData(upp7g)[ ,"e.rfs"], pData(vdx7g)[ ,"e.dmfs"], pData(nki7g)[ ,"e.dmfs"])            
stat1.exprs <- c(exprs(mainz7g)[stat1Paf,], exprs(transbig7g)[stat1Paf,], exprs(unt7g)[stat1Paf,], exprs(upp7g)[stat1Paf,], exprs(vdx7g)[stat1Paf,], exprs(nki7g)[stat1Pagi,])

stat1.exprs.length <- c(length(exprs(mainz7g)[stat1Paf,]), length(exprs(transbig7g)[stat1Paf,]), length(exprs(unt7g)[stat1Paf,]), length(exprs(upp7g)[stat1Paf,]), length(exprs(vdx7g)[stat1Paf,]), length(exprs(nki7g)[stat1Pagi,]))

stat1.exprs.group <- c(as.character(scmgene.mainz$subtype2), as.character(scmgene.transbig$subtype2), as.character(scmgene.upp$subtype2), as.character(scmgene.unt$subtype2), as.character(scmgene.vdx$subtype2), as.character(scmgene.nki$subtype2))

#mygroup <- as.factor(stat1.exprs.group)      

mm <- stat1.exprs.group=="ER-/HER2-"
pmh <- stat1.exprs.group=="ER+/HER2- High Prolif"
pml <- stat1.exprs.group=="ER+/HER2- Low Prolif"
her2p <- stat1.exprs.group=="HER2+"

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
gg <- factor(c(rep("mainz", nrow(pData(mainz7g))), rep("transbig", nrow(pData(transbig7g))), rep("unt", nrow(pData(unt7g))), rep("upp", nrow(pData(upp7g))), rep("vdx", nrow(pData(vdx7g))), rep("nki", nrow(pData(nki7g)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ mygroup), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="Probability of survival for STAT1", sub.title="", leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))

dev.off()

dev.copy2eps(file="KMPlot-STAT1-subtypes.eps")











