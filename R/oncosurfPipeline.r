#######################
##
##List of genes: ESR1, ERBB2, AURKA, PLAU, VEGF, STAT1, CASP3
##List of geneids: 2099, 2064, 6790, 5328, 7422, 6772, 836
##List of probes:
##VDX/TRANSBIG/UPP/UNT/MAINSZ: 205225_at, 216836_s_at, 208079_s_at, 211668_s_at, 211527_x_at, 209969_s_at, 202763_at
##NKI: NM_000125, NM_004448, NM_003600, NM_002658, NM_003376, NM_007315, NM_004346
#######################

#######################
## get packages and dependcies
#######################
#install.packages("../../git/survcomp_1.1.10.tar.gz")
library(survcomp)
library(genefu)
library(Biobase)
library(Hmisc)


#######################
## load eSets
#######################
load("breast-cancer_exp-packages/breastCancerNKI/data/nki.rda")
load("breast-cancer_exp-packages/breastCancerTRANSBIG/data/transbig.rda")
load("breast-cancer_exp-packages/breastCancerUPP/data/upp.rda")
load("breast-cancer_exp-packages/breastCancerUNT/data/unt.rda")
load("breast-cancer_exp-packages/breastCancerVDX/data/vdx.rda")
load("breast-cancer_exp-packages/breastCancerMAINZ/data/mainz.rda")
                       
#######################
## set variables
#######################

##Oncotype DX
##gs.list <- c("ACTB","GAPDH","GUS","RPLPO","TFRC","GRB7","ER","HER2","PGR","BCL2","SCUBE2","Survivin","KI67","MYBL2","CCNB1","STK15","CTSL2","MMP11","GSTM1","BAG1","CD68")
## gid.list <- c(7037,60,5241,4320,4605,573,2597,968,891,2944,1515,2886,596,57758)
## geneid.list <- c("geneid.7037","geneid.60","geneid.5241","geneid.4320","geneid.4605","geneid.573","geneid.2597","geneid.891","geneid.2944","geneid.1515","geneid.2886","geneid.596","geneid.57758")
## gs.list <- c("TFRC","PGR","MMP11","MYBL2","BAG1","GAPDH","CCNB1","GSTM1","CTSL2","GRB7","BCL2","SCUBE2")
## gs.list <- c("TFRC","ACTB","PGR","MMP11","MYBL2","BAG1","GAPDH","<NA>","CCNB1","GSTM1","CTSL2","GRB7","BCL2","SCUBE2")
## geneid.968 has no gene symbol in the annotation data (for mainz)

gs.list <- c("ESR1", "ERBB2", "AURKA", "PLAU", "VEGF", "STAT1", "CASP3")
gid.list <- c(2099, 2064, 6790, 5328, 7422, 6772, 836)
geneid.list <- c("geneid.2099","geneid.2064","geneid.6790","geneid.5328","geneid.7422","geneid.6772","geneid.836")
p.nki.list <- c("NM_000125", "NM_004448", "NM_003600", "NM_002658", "NM_003376", "NM_007315", "NM_004346")
p.affy.list <- c("205225_at", "216836_s_at", "208079_s_at", "211668_s_at", "211527_x_at", "209969_s_at", "202763_at")
##affy.data <- c("mainz", "transbig", "unt", "upp", "vdx")
##agilent.data <- c("nki")
##dataset.list <- c("mainz","transbig","upp","unt","vdx","nki","all")
dataset.list <- c("MAINZ","TRANSBIG","UPP","UNT","VDX","NKI","ALL")
myspace <- " "
mybigspace <- "    "
tc <- 10 * 365


pdf("all-figures_with-mapping.pdf")

#######################
## create eSets only for the 7 genes
#######################
## mainz
## mapping for mainz
gid <- sort(unique(as.numeric(as.character(featureData(mainz)@data[ ,"EntrezGene.ID"]))))
names(gid) <- paste("geneid", gid, sep=".")
gid.orig <- as.numeric(as.character(featureData(mainz)@data[ ,"EntrezGene.ID"]))
names(gid.orig) <- dimnames(featureData(mainz)@data)[[1]]
data.orig <- t(exprs(mainz))
rr <- geneid.map(geneid1=gid.orig, data1=data.orig, geneid2=gid)
## create a new (reduced) mainz eSet
## annot
tt <- featureData(mainz)@data[names(rr$geneid1), , drop=FALSE]
dimnames(tt)[[1]] <- paste("geneid", rr$geneid1, sep=".")
metadata<-data.frame(labelDescription=colnames(tt), row.names=colnames(tt))  
featureD <- new("AnnotatedDataFrame", data=tt, varMetadata=metadata)
## gene expression
tt <- rr$data1
dimnames(tt)[[2]] <- paste("geneid", rr$geneid1, sep=".")
arrayT <- "hgu133a.geneid"
mainz.geneid <-new("ExpressionSet", 
          exprs=t(tt), 
          phenoData=phenoData(mainz),
          featureData=featureD, 
          annotation=arrayT)
          
geneid.fData <- fData(mainz.geneid)[geneid.list,]
geneid.eData <- exprs(mainz.geneid)[geneid.list,]
metadata<-data.frame(labelDescription=colnames(geneid.fData), row.names=colnames(geneid.fData))  
featureD <- new("AnnotatedDataFrame", data=geneid.fData, varMetadata=metadata)
mainz.small <-new("ExpressionSet",
        exprs=geneid.eData, 
        phenoData=phenoData(mainz.geneid),
        featureData=featureD, 
        annotation=annotation(mainz.geneid),
        experimentData=experimentData(mainz.geneid))
        
        
        
## transbig
## mapping for transbig
gid <- sort(unique(as.numeric(as.character(featureData(transbig)@data[ ,"EntrezGene.ID"]))))
names(gid) <- paste("geneid", gid, sep=".")
gid.orig <- as.numeric(as.character(featureData(transbig)@data[ ,"EntrezGene.ID"]))
names(gid.orig) <- dimnames(featureData(transbig)@data)[[1]]
data.orig <- t(exprs(transbig))
rr <- geneid.map(geneid1=gid.orig, data1=data.orig, geneid2=gid)
## create a new (reduced) transbig eSet
## annot
tt <- featureData(transbig)@data[names(rr$geneid1), , drop=FALSE]
dimnames(tt)[[1]] <- paste("geneid", rr$geneid1, sep=".")
metadata<-data.frame(labelDescription=colnames(tt), row.names=colnames(tt))  
featureD <- new("AnnotatedDataFrame", data=tt, varMetadata=metadata)
## gene expression
tt <- rr$data1
dimnames(tt)[[2]] <- paste("geneid", rr$geneid1, sep=".")
arrayT <- "hgu133a.geneid"
transbig.geneid <-new("ExpressionSet", 
          exprs=t(tt), 
          phenoData=phenoData(transbig),
          featureData=featureD, 
          annotation=arrayT)
          
geneid.fData <- fData(transbig.geneid)[geneid.list,]
geneid.eData <- exprs(transbig.geneid)[geneid.list,]
metadata<-data.frame(labelDescription=colnames(geneid.fData), row.names=colnames(geneid.fData))  
featureD <- new("AnnotatedDataFrame", data=geneid.fData, varMetadata=metadata)
transbig.small <-new("ExpressionSet",
        exprs=geneid.eData, 
        phenoData=phenoData(transbig.geneid),
        featureData=featureD, 
        annotation=annotation(transbig.geneid),
        experimentData=experimentData(transbig.geneid))


## vdx
## mapping for vdx
gid <- sort(unique(as.numeric(as.character(featureData(vdx)@data[ ,"EntrezGene.ID"]))))
names(gid) <- paste("geneid", gid, sep=".")
gid.orig <- as.numeric(as.character(featureData(vdx)@data[ ,"EntrezGene.ID"]))
names(gid.orig) <- dimnames(featureData(vdx)@data)[[1]]
data.orig <- t(exprs(vdx))
rr <- geneid.map(geneid1=gid.orig, data1=data.orig, geneid2=gid)
## create a new (reduced) vdx eSet
## annot
tt <- featureData(vdx)@data[names(rr$geneid1), , drop=FALSE]
dimnames(tt)[[1]] <- paste("geneid", rr$geneid1, sep=".")
metadata<-data.frame(labelDescription=colnames(tt), row.names=colnames(tt))  
featureD <- new("AnnotatedDataFrame", data=tt, varMetadata=metadata)
## gene expression
tt <- rr$data1
dimnames(tt)[[2]] <- paste("geneid", rr$geneid1, sep=".")
arrayT <- "hgu133a.geneid"
vdx.geneid <-new("ExpressionSet", 
          exprs=t(tt), 
          phenoData=phenoData(vdx),
          featureData=featureD, 
          annotation=arrayT)

geneid.fData <- fData(vdx.geneid)[geneid.list,]
geneid.eData <- exprs(vdx.geneid)[geneid.list,]
metadata<-data.frame(labelDescription=colnames(geneid.fData), row.names=colnames(geneid.fData))  
featureD <- new("AnnotatedDataFrame", data=geneid.fData, varMetadata=metadata)
vdx.small <-new("ExpressionSet",
        exprs=geneid.eData, 
        phenoData=phenoData(vdx.geneid),
        featureData=featureD, 
        annotation=annotation(vdx.geneid),
        experimentData=experimentData(vdx.geneid))



## upp
## mapping for upp
gid <- sort(unique(as.numeric(as.character(featureData(upp)@data[ ,"EntrezGene.ID"]))))
names(gid) <- paste("geneid", gid, sep=".")
gid.orig <- as.numeric(as.character(featureData(upp)@data[ ,"EntrezGene.ID"]))
names(gid.orig) <- dimnames(featureData(upp)@data)[[1]]
data.orig <- t(exprs(upp))
rr <- geneid.map(geneid1=gid.orig, data1=data.orig, geneid2=gid)
## create a new (reduced) upp eSet
## annot
tt <- featureData(upp)@data[names(rr$geneid1), , drop=FALSE]
dimnames(tt)[[1]] <- paste("geneid", rr$geneid1, sep=".")
metadata<-data.frame(labelDescription=colnames(tt), row.names=colnames(tt))  
featureD <- new("AnnotatedDataFrame", data=tt, varMetadata=metadata)
## gene expression
tt <- rr$data1
dimnames(tt)[[2]] <- paste("geneid", rr$geneid1, sep=".")
arrayT <- "hgu133a.geneid"
upp.geneid <-new("ExpressionSet", 
          exprs=t(tt), 
          phenoData=phenoData(upp),
          featureData=featureD, 
          annotation=arrayT)

geneid.fData <- fData(upp.geneid)[geneid.list,]
geneid.eData <- exprs(upp.geneid)[geneid.list,]
metadata<-data.frame(labelDescription=colnames(geneid.fData), row.names=colnames(geneid.fData))  
featureD <- new("AnnotatedDataFrame", data=geneid.fData, varMetadata=metadata)
upp.small <-new("ExpressionSet",
        exprs=geneid.eData, 
        phenoData=phenoData(upp.geneid),
        featureData=featureD, 
        annotation=annotation(upp.geneid),
        experimentData=experimentData(upp.geneid))



## unt
## mapping for unt
gid <- sort(unique(as.numeric(as.character(featureData(unt)@data[ ,"EntrezGene.ID"]))))
names(gid) <- paste("geneid", gid, sep=".")
gid.orig <- as.numeric(as.character(featureData(unt)@data[ ,"EntrezGene.ID"]))
names(gid.orig) <- dimnames(featureData(unt)@data)[[1]]
data.orig <- t(exprs(unt))
rr <- geneid.map(geneid1=gid.orig, data1=data.orig, geneid2=gid)
## create a new (reduced) unt eSet
## annot
tt <- featureData(unt)@data[names(rr$geneid1), , drop=FALSE]
dimnames(tt)[[1]] <- paste("geneid", rr$geneid1, sep=".")
metadata<-data.frame(labelDescription=colnames(tt), row.names=colnames(tt))  
featureD <- new("AnnotatedDataFrame", data=tt, varMetadata=metadata)
## gene expression
tt <- rr$data1
dimnames(tt)[[2]] <- paste("geneid", rr$geneid1, sep=".")
arrayT <- "hgu133a.geneid"
unt.geneid <-new("ExpressionSet", 
          exprs=t(tt), 
          phenoData=phenoData(unt),
          featureData=featureD, 
          annotation=arrayT)

geneid.fData <- fData(unt.geneid)[geneid.list,]
geneid.eData <- exprs(unt.geneid)[geneid.list,]
metadata<-data.frame(labelDescription=colnames(geneid.fData), row.names=colnames(geneid.fData))  
featureD <- new("AnnotatedDataFrame", data=geneid.fData, varMetadata=metadata)
unt.small <-new("ExpressionSet",
        exprs=geneid.eData, 
        phenoData=phenoData(unt.geneid),
        featureData=featureD, 
        annotation=annotation(unt.geneid),
        experimentData=experimentData(unt.geneid))



## nki
## mapping for nki
gid <- sort(unique(as.numeric(as.character(featureData(nki)@data[ ,"EntrezGene.ID"]))))
names(gid) <- paste("geneid", gid, sep=".")
gid.orig <- as.numeric(as.character(featureData(nki)@data[ ,"EntrezGene.ID"]))
names(gid.orig) <- dimnames(featureData(nki)@data)[[1]]
data.orig <- t(exprs(nki))
rr <- geneid.map(geneid1=gid.orig, data1=data.orig, geneid2=gid)
## create a new (reduced) nki eSet
## annot
tt <- featureData(nki)@data[names(rr$geneid1), , drop=FALSE]
dimnames(tt)[[1]] <- paste("geneid", rr$geneid1, sep=".")
metadata<-data.frame(labelDescription=colnames(tt), row.names=colnames(tt))  
featureD <- new("AnnotatedDataFrame", data=tt, varMetadata=metadata)
## gene expression
tt <- rr$data1
dimnames(tt)[[2]] <- paste("geneid", rr$geneid1, sep=".")
arrayT <- "hgu133a.geneid"
nki.geneid <-new("ExpressionSet", 
          exprs=t(tt), 
          phenoData=phenoData(nki),
          featureData=featureD, 
          annotation=arrayT)

geneid.fData <- fData(nki.geneid)[geneid.list,]
geneid.eData <- exprs(nki.geneid)[geneid.list,]
metadata<-data.frame(labelDescription=colnames(geneid.fData), row.names=colnames(geneid.fData))  
featureD <- new("AnnotatedDataFrame", data=geneid.fData, varMetadata=metadata)
nki.small <-new("ExpressionSet",
        exprs=geneid.eData, 
        phenoData=phenoData(nki.geneid),
        featureData=featureD, 
        annotation=annotation(nki.geneid),
        experimentData=experimentData(nki.geneid))


#save(list=c("mainz","transbig","upp","unt","vdx","nki"), compress=TRUE, file=paste("sampleData.RData",sep=""))
#mainz.small <- mainz
#transbig.small <- transbig
#nki.small <- nki
#unt.small <- unt
#upp.small <- upp
#vdx.small <- vdx



#######################
## concordance.index computation
#######################
cindexall.mainz.small <- t(apply(X=exprs(mainz.small), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainz.small)[ ,"t.dmfs"], z=pData(mainz.small)[ ,"e.dmfs"]))

cindexall.transbig.small <- t(apply(X=exprs(transbig.small), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbig.small)[ ,"t.dmfs"], z=pData(transbig.small)[ ,"e.dmfs"]))

cindexall.upp.small <- t(apply(X=exprs(upp.small), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(upp.small)[ ,"t.rfs"], z=pData(upp.small)[ ,"e.rfs"]))

cindexall.unt.small <- t(apply(X=exprs(unt.small), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(unt.small)[ ,"t.dmfs"], z=pData(unt.small)[ ,"e.dmfs"]))

cindexall.vdx.small <- t(apply(X=exprs(vdx.small), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdx.small)[ ,"t.dmfs"], z=pData(vdx.small)[ ,"e.dmfs"])) 

cindexall.nki.small <- t(apply(X=exprs(nki.small), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nki.small)[ ,"t.dmfs"], z=pData(nki.small)[ ,"e.dmfs"]))



#######################
## D.index computation
#######################
dindexall.mainz.small <- t(apply(X=exprs(mainzSample), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainzSample)[ ,"t.dmfs"], z=pData(mainzSample)[ ,"e.dmfs"]))

dindexall.transbig.small <- t(apply(X=exprs(transbigSample), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbigSample)[ ,"t.dmfs"], z=pData(transbigSample)[ ,"e.dmfs"]))

dindexall.upp.small <- t(apply(X=exprs(uppSample), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(uppSample)[ ,"t.rfs"], z=pData(uppSample)[ ,"e.rfs"]))

dindexall.unt.small <- t(apply(X=exprs(untSample), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(untSample)[ ,"t.dmfs"], z=pData(untSample)[ ,"e.dmfs"]))

dindexall.vdx.small <- t(apply(X=exprs(vdxSample), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdxSample)[ ,"t.dmfs"], z=pData(vdxSample)[ ,"e.dmfs"]))

dindexall.nki.small <- t(apply(X=exprs(nkiSample), MARGIN=1, function(x, y, z) { tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nkiSample)[ ,"t.dmfs"], z=pData(nkiSample)[ ,"e.dmfs"]))



#######################
## hazard.ratio computation
#######################
hratio.mainz.small <- t(apply(X=exprs(mainzSample), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainzSample)[ ,"t.dmfs"], z=pData(mainzSample)[ ,"e.dmfs"]))

hratio.transbig.small <- t(apply(X=exprs(transbigSample), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbigSample)[ ,"t.dmfs"], z=pData(transbigSample)[ ,"e.dmfs"]))

hratio.upp.small <- t(apply(X=exprs(uppSample), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(uppSample)[ ,"t.rfs"], z=pData(uppSample)[ ,"e.rfs"]))

hratio.unt.small <- t(apply(X=exprs(untSample), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(untSample)[ ,"t.dmfs"], z=pData(untSample)[ ,"e.dmfs"]))

hratio.vdx.small <- t(apply(X=exprs(vdxSample), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdxSample)[ ,"t.dmfs"], z=pData(vdxSample)[ ,"e.dmfs"]))

hratio.nki.small <- t(apply(X=exprs(nkiSample), MARGIN=1, function(x, y, z) { tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE); return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nkiSample)[ ,"t.dmfs"], z=pData(nkiSample)[ ,"e.dmfs"]))
        


#######################
## CINDEX: combining the cindices for each gene from all datasets
#######################        
tt <- as.data.frame(NULL)
for(i in 1:length(gs.list)){
  tt <- rbind(
    tt,combine.est(x=cbind(   cindexall.mainz.small[i,"cindex"],
                              cindexall.transbig.small[i,"cindex"],
                              cindexall.upp.small[i,"cindex"],
                              cindexall.unt.small[i,"cindex"],
                              cindexall.vdx.small[i,"cindex"],
                              cindexall.nki.small[i,"cindex"]),
                   x.se=cbind(cindexall.mainz.small[i,"cindex"],
                              cindexall.transbig.small[i,"cindex.se"],
                              cindexall.upp.small[i,"cindex.se"],
                              cindexall.unt.small[i,"cindex.se"],
                              cindexall.vdx.small[i,"cindex.se"],
                              cindexall.nki.small[i,"cindex.se"]),)
              )
}
tt$lower <- tt$estimate + qnorm(0.025, lower.tail=TRUE) * tt$se
tt$upper <- tt$estimate + qnorm(0.025, lower.tail=FALSE) * tt$se
rownames(tt) <- gs.list
colnames(tt) <- c("cindex","cindex.se","lower","upper")
ccindex <- tt

        
        
#######################
## DINDEX: combining the dindices for each gene from all datasets
#######################        
tt <- as.data.frame(NULL)
for(i in 1:length(gs.list)){
  tt <- rbind(
    tt,combine.est(x=cbind(   dindexall.mainz.small[i,"dindex"],
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
rownames(tt) <- gs.list
colnames(tt) <- c("dindex","dindex.se","lower","upper")
cdindex <- tt



#######################
## HRATIO: combining the hazard ratios for each gene from all datasets
#######################        
tt <- as.data.frame(NULL)
for(i in 1:length(gs.list)){
  tt <- rbind(
    tt,combine.est(x=cbind(   hratio.mainz.small[i,"hratio"],
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
rownames(tt) <- gs.list
colnames(tt) <- c("hratio","hratio.se","lower","upper")
chratio <- tt
        
        

#######################
## CINDEX: forestplot for each gene combining all datasets
#######################
##pdf("forestplot-cindex-all-genes.pdf")                  
labeltext <- cbind(c("Gene Symbol",gs.list),c(rep(myspace,length(gs.list))))
bs <- rep(0.5, nrow(labeltext))                              
r.mean <- c(NA,ccindex$cindex)
r.lower <- c(NA,ccindex$cindex + qnorm(0.025, lower.tail=TRUE) * ccindex$cindex.se)
r.upper <- c(NA,ccindex$cindex + qnorm(0.025, lower.tail=FALSE) * ccindex$cindex.se)

forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.3,0.7,0.05), xlab=paste("concordance index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"), box.size=bs, clip=c(0.3,1))
title(paste("cindex forestplot, each gene, ccindex"))
##dev.off()



#######################
## DINDEX: forestplot for each gene combining all datasets
#######################
##pdf("forestplot-dindex-all-genes.pdf")                  
labeltext <- cbind(c("Gene Symbol",gs.list),c(rep(mybigspace,length(gs.list))),c(rep(mybigspace,length(gs.list))))
bs <- rep(0.5, nrow(labeltext))                              
tt <- log2(cdindex)

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,1.5,0.1), xlab=paste("log2 D.index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"), box.size=bs, clip=c(-0.5,1.5))
title(paste("log2 dindex forestplot, each gene, cdindex"))
##dev.off()



#######################
## HRATIO: forestplot for each gene combining all datasets
#######################
##pdf("forestplot-hratio-all-genes.pdf")                  
labeltext <- cbind(c("Gene Symbol",gs.list),c(rep(mybigspace,length(gs.list))),c(rep(mybigspace,length(gs.list))))
bs <- rep(0.5, nrow(labeltext))                              
tt <- log2(chratio)

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$hratio), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,1,0.1), xlab=paste("log2 hazard.ratio", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"), box.size=bs, clip=c(-0.5,1))
title(paste("log2 hratio forestplot, each gene, chratio"))
##dev.off()


##pdf("AURKA-and-VEGF-meta-analysis.pdf")



#######################
## CINDEX: AURKA and different datasets as forestplot
#######################
##pdf("AURKA-all-datasets-and-combined-cindex-dindex-hratio.pdf")
tt <- rbind(cindexall.mainz.small[3,],
            cindexall.transbig.small[3,],
            cindexall.upp.small[3,],
            cindexall.unt.small[3,],
            cindexall.vdx.small[3,],
            cindexall.nki.small[3,],
            as.numeric(ccindex[3,]))

rownames(tt) <- dataset.list
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",dataset.list),c(mybigspace,rep(mybigspace,length(dataset.list))))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.8,0.05), xlab=paste("AURKA cindex", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"), box.size=bs, clip=c(0.5,1), is.summary=(c(rep(FALSE,7),TRUE)))
title(paste("cindex forestplot, AURKA"))



#######################
## CINDEX: VEGF and different datasets as forestplot
#######################
##pdf("VEGF-all-datasets-and-combined-cindex-dindex-hratio.pdf")
tt <- rbind(cindexall.mainz.small[5,],
            cindexall.transbig.small[5,],
            cindexall.upp.small[5,],
            cindexall.unt.small[5,],
            cindexall.vdx.small[5,],
            cindexall.nki.small[5,],
            as.numeric(ccindex[5,]))

rownames(tt) <- dataset.list
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",dataset.list),c(mybigspace,rep(mybigspace,length(dataset.list))))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.75,0.05), xlab=paste("VEGF cindex", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"), box.size=bs, clip=c(0.4,0.75), is.summary=(c(rep(FALSE,7),TRUE)))
title(paste("cindex forestplot, VEGF"))


#######################
## CINDEX: AURKA and VEGF together with different datasets as forestplot
#######################
##pdf("AURKA-VEGF-all-datasets-and-combined-cindex-dindex-hratio.pdf")
tt <- rbind(cindexall.mainz.small[3,],
            cindexall.mainz.small[5,],
            NA,
            cindexall.transbig.small[3,],            
            cindexall.transbig.small[5,],
            NA,
            cindexall.upp.small[3,],            
            cindexall.upp.small[5,],
            NA,
            cindexall.unt.small[3,],            
            cindexall.unt.small[5,],
            NA,
            cindexall.vdx.small[3,],            
            cindexall.vdx.small[5,],
            NA,
            cindexall.nki.small[3,],
            cindexall.nki.small[5,],
            NA,
            as.numeric(ccindex[3,]),            
            as.numeric(ccindex[5,]))

rownames(tt) <- c("MAINZa","MAINZv","a","TRANSBIGa","TRANSBIGv","b","UPPa","UPPv","c","UNTa","UNTv","d","VDXa","VDXv","e","NKIa","NKIv","f","ALLa","ALLv")
##dataset.list
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset","MAINZ",NA,NA,"TRANSBIG",NA,NA,"UPP",NA,NA,"UNT",NA,NA,"VDX",NA,NA,"NKI",NA,NA,"ALL",NA),
                   c("Gene",rep(c("aurka","vegf",NA),length(dataset.list)-1),c("aurka","vegf")),
                   c(rep(mybigspace,21)))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.8,0.05), xlab=paste("AURKA and VEGF cindex", myspace, sep=""), col=meta.colors(line=c(rep(c(NA,"darkblue","seagreen"),7)),zero="firebrick",box=c(rep(c(NA,"royalblue","forestgreen"),7))), box.size=bs, clip=c(0.4,1), is.summary=(c(rep(FALSE,19),TRUE,TRUE)))
title(paste("cindex forestplot, AURKA and VEGF"))    

                         

#######################
## DINDEX: AURKA and different datasets as forestplot
#######################
tt <- rbind(dindexall.mainz.small[3,],
            dindexall.transbig.small[3,],
            dindexall.upp.small[3,],
            dindexall.unt.small[3,],
            dindexall.vdx.small[3,],
            dindexall.nki.small[3,],
            as.numeric(cdindex[3,]))

rownames(tt) <- dataset.list
tt <- as.data.frame(tt)
tt <- log2(tt)
labeltext <- cbind(c("Dataset",dataset.list),c(mybigspace,rep(mybigspace,length(dataset.list))))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,2,0.5), xlab=paste("AURKA log2 D.index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"), box.size=bs, clip=c(-0.25,2), is.summary=(c(rep(FALSE,7),TRUE)))
title(paste("log2 D.index forestplot, AURKA"))


#######################
## DINDEX: VEGF and different datasets as forestplot
#######################
tt <- rbind(dindexall.mainz.small[5,],
            dindexall.transbig.small[5,],
            dindexall.upp.small[5,],
            dindexall.unt.small[5,],
            dindexall.vdx.small[5,],
            dindexall.nki.small[5,],
            as.numeric(cdindex[5,]))

rownames(tt) <- dataset.list
tt <- as.data.frame(tt)
tt <- log2(tt)
labeltext <- cbind(c("Dataset",dataset.list),c(mybigspace,rep(mybigspace,length(dataset.list))))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-1,1.5,0.5), xlab=paste("VEGF log2 D.index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="firebrick"), box.size=bs, clip=c(-1,1.5), is.summary=(c(rep(FALSE,7),TRUE)))
title(paste("log2 D.index forestplot, VEGF"))


#######################
## DINDEX: AURKA and VEGF together with different datasets as forestplot
#######################
##pdf("AURKA-VEGF-all-datasets-and-combined-cindex-dindex-hratio.pdf")
tt <- rbind(dindexall.mainz.small[3,],
            dindexall.mainz.small[5,],
            NA,
            dindexall.transbig.small[3,],            
            dindexall.transbig.small[5,],
            NA,
            dindexall.upp.small[3,],            
            dindexall.upp.small[5,],
            NA,
            dindexall.unt.small[3,],            
            dindexall.unt.small[5,],
            NA,
            dindexall.vdx.small[3,],            
            dindexall.vdx.small[5,],
            NA,
            dindexall.nki.small[3,],
            dindexall.nki.small[5,],
            NA,
            as.numeric(cdindex[3,]),            
            as.numeric(cdindex[5,]))

rownames(tt) <- c("MAINZa","MAINZv","a","TRANSBIGa","TRANSBIGv","b","UPPa","UPPv","c","UNTa","UNTv","d","VDXa","VDXv","e","NKIa","NKIv","f","ALLa","ALLv")
##dataset.list
tt <- as.data.frame(tt)
tt <- log2(tt)
labeltext <- cbind(c("Dataset","MAINZ",NA,NA,"TRANSBIG",NA,NA,"UPP",NA,NA,"UNT",NA,NA,"VDX",NA,NA,"NKI",NA,NA,"ALL",NA),
                   c("Gene",rep(c("aurka","vegf",NA),length(dataset.list)-1),c("aurka","vegf")),
                   c(rep(mybigspace,21)))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-1,2,0.5), xlab=paste("AURKA and VEGF log2 D.index", myspace, sep=""), col=meta.colors(line=c(rep(c(NA,"darkblue","seagreen"),7)),zero="firebrick",box=c(rep(c(NA,"royalblue","forestgreen"),7))), box.size=bs, clip=c(-1,2), is.summary=(c(rep(FALSE,19),TRUE,TRUE)))
title(paste("log2 D.index forestplot, AURKA and VEGF"))    

                         




#######################
## HRATIO: AURKA and different datasets as forestplot
#######################
tt <- rbind(hratio.mainz.small[3,],
            hratio.transbig.small[3,],
            hratio.upp.small[3,],
            hratio.unt.small[3,],
            hratio.vdx.small[3,],
            hratio.nki.small[3,],
            as.numeric(chratio[3,]))

rownames(tt) <- dataset.list
tt <- as.data.frame(tt)
tt <- log2(tt)
labeltext <- cbind(c("Dataset",dataset.list),c(rep(myspace,length(dataset.list))))
bs <- rep(0.5, nrow(labeltext))
   
forestplot.surv(labeltext=labeltext, mean=c(NA,tt$hratio), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,3.5,0.5), xlab=paste("AURKA log2 hazard.ratio", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.5,3.5),is.summary=(c(rep(FALSE,7),TRUE)))
title(paste("log2 hratio forestplot, AURKA"))
##dev.off()


          
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



#######################
## DINDEX: every SINGLE gene and the different datasets as forestplot
#######################
##pdf("dindex_forestplot-for-each-gene-showing-all-datasets.pdf")
for(i in 1:length(dataset.list)) {
  myspace <- " "
  tt <- rbind(dindexall.mainz.small[i,],
              dindexall.transbig.small[i,],
              dindexall.upp.small[i,],
              dindexall.unt.small[i,],
              dindexall.vdx.small[i,],
              dindexall.nki.small[i,],
              as.numeric(cdindex[i,]))

  rownames(tt) <- dataset.list
  tt <- as.data.frame(tt)
  tt <- log2(tt)
  labeltext <- cbind(c("Dataset",dataset.list),c(rep(mybigspace,length(dataset.list))))
  bs <- rep(0.5, nrow(labeltext))
   
  forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-1,2,0.25), xlab=paste(gs.list[i], myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-1,2),is.summary=(c(rep(FALSE,7),TRUE)))
  title(paste("log2 dindex forestplot, ", gs.list[i]))
}
##dev.off()


#######################
## HRATIO: every SINGLE gene and the different datasets as forestplot
#######################
##pdf("hratio_forestplot-for-each-gene-showing-all-datasets.pdf")
for(i in 1:length(dataset.list)) {
  myspace <- " "
  tt <- rbind(hratio.mainz.small[i,],
              hratio.transbig.small[i,],
              hratio.upp.small[i,],
              hratio.unt.small[i,],
              hratio.vdx.small[i,],
              hratio.nki.small[i,],
              as.numeric(chratio[i,]))

  rownames(tt) <- dataset.list
  tt <- as.data.frame(tt)
  tt <- log2(tt)
  labeltext <- cbind(c("Dataset",dataset.list),c(rep(myspace,length(dataset.list))))
  bs <- rep(0.5, nrow(labeltext))
   
  forestplot.surv(labeltext=labeltext, mean=c(NA,tt$hratio), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-1,3.5,0.5), xlab=paste(gs.list[i], myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.75,3.5),is.summary=(c(rep(FALSE,7),TRUE)))
  title(paste("log2 hratio forestplot, ", gs.list[i]))
}
##dev.off()



#######################
## concatenate datasets and make survival curve
#######################     
surv.data <- censor.time(surv.time=c(pData(mainz.small)[ ,"t.dmfs"], pData(transbig.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(mainz.small)[ ,"e.dmfs"], pData(transbig.small)[ ,"e.rfs"]), time.cens=tc / 365)
gg <- factor(c(rep("mainz", nrow(pData(mainz.small))), rep("transbig", nrow(pData(transbig.small)))), levels=c("mainz", "transbig"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkorange", "darkviolet"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
title(paste("km survival curve, MAINZ and TRANSBIG"))



#######################
## draw survival curve for all datasets and one survival curve for each datasets
#######################
##pdf("survival-curve-test_all-small-datasets.pdf")
surv.data <- censor.time(surv.time=c(pData(mainz.small)[ ,"t.dmfs"], pData(transbig.small)[ ,"t.dmfs"], pData(unt.small)[ ,"t.dmfs"], pData(vdx.small)[ ,"t.dmfs"], pData(upp.small)[ ,"t.rfs"], pData(nki.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(mainz.small)[ ,"e.dmfs"], pData(transbig.small)[ ,"e.dmfs"], pData(unt.small)[ ,"e.dmfs"], pData(vdx.small)[ ,"e.dmfs"], pData(upp.small)[ ,"e.rfs"], pData(nki.small)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("mainz", nrow(pData(mainz.small))), rep("transbig", nrow(pData(transbig.small))), rep("unt", nrow(pData(unt.small))), rep("vdx", nrow(pData(vdx.small))), rep("upp", nrow(pData(upp.small))), rep("nki", nrow(pData(nki.small)))), levels=c("mainz", "transbig", "unt", "vdx", "upp", "nki"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS/RFS(upp)", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkorange", "red", "darkblue", "darkgreen", "black", "brown"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
title(paste("km survival curve, all datasets"))



#######################
## survival curve for each dataset
#######################
surv.data <- censor.time(surv.time=c(pData(mainz.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(mainz.small)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("mainz", nrow(pData(mainz.small)))), levels=c("mainz"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkorange"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
title(paste("km survival curve, MAINZ"))

surv.data <- censor.time(surv.time=c(pData(transbig.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(transbig.small)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("transbig", nrow(pData(transbig.small)))), levels=c("transbig"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("red"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
title(paste("km survival curve, TRANSBIG"))

surv.data <- censor.time(surv.time=c(pData(upp.small)[ ,"t.rfs"]) / 365, surv.event=c(pData(upp.small)[ ,"e.rfs"]), time.cens=tc / 365)
gg <- factor(c(rep("upp", nrow(pData(upp.small)))), levels=c("upp"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of RFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("brown"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE) 
title(paste("km survival curve, UPP"))

surv.data <- censor.time(surv.time=c(pData(unt.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(unt.small)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("unt", nrow(pData(unt.small)))), levels=c("unt"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("black"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
title(paste("km survival curve, UNT"))

surv.data <- censor.time(surv.time=c(pData(vdx.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(vdx.small)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("vdx", nrow(pData(vdx.small)))), levels=c("vdx"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkblue"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)
title(paste("km survival curve, VDX"))

surv.data <- censor.time(surv.time=c(pData(nki.small)[ ,"t.dmfs"]) / 365, surv.event=c(pData(nki.small)[ ,"e.dmfs"]), time.cens=tc / 365)
gg <- factor(c(rep("nki", nrow(pData(nki.small)))), levels=c("nki"))
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "group"=gg)
km.coxph.plot(formula.s=formula(Surv(time, event) ~ group), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of DMFS", main.title="", sub.title=NULL, leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .lty=rep(1, length(levels(gg))), show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, .col=c("darkgreen"), leg.text=paste(levels(gg), myspace, sep=""), verbose=FALSE)      
title(paste("km survival curve, NKI"))
##dev.off()                              



#######################
## censor.time,  km.coxph.plot
#######################

#gene symbol: AURKA
#entrezgene id: 6790
#VDX/TRANSBIG/UPP/UNT/MAINSZ probe id: 208079_s_at
#NKI probe id: NM_003600

aurkaGs <- "AURKA"
aurkaGid <- 6790
aurkaPaf <- "208079_s_at"
aurkaPagi <- "NM_003600"


surv.time.all <- c(pData(mainz.small)[ ,"t.dmfs"], pData(transbig.small)[ ,"t.dmfs"], pData(unt.small)[ ,"t.dmfs"], pData(vdx.small)[ ,"t.dmfs"], pData(upp.small)[ ,"t.rfs"], pData(nki.small)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz.small)[ ,"e.dmfs"], pData(transbig.small)[ ,"e.dmfs"], pData(unt.small)[ ,"e.dmfs"], pData(vdx.small)[ ,"e.dmfs"], pData(upp.small)[ ,"e.rfs"], pData(nki.small)[ ,"e.dmfs"])            
aurka.exprs <- c(exprs(mainz)[aurkaPaf,], exprs(transbig)[aurkaPaf,], exprs(unt)[aurkaPaf,], exprs(vdx)[aurkaPaf,], exprs(upp)[aurkaPaf,], exprs(nki)[aurkaPagi,])


#quantile for cutoff
#34% in the wang2005geneexpression publication
#41% in the foekens2006multicenter publication
#28% in the desmedt2007strong publication
##myprobs <- c(0.33, 0.66)
##mycutoff <- quantile(aurka.exprs, probs=myprobs, na.rm=TRUE)
##mycutoff <- cut2(aurka.exprs, myprobs)
##cut2(x=aurka.exprs, cuts=quantile(aurka.exprs, myprobs), g=3)

qq <- aurka.exprs
qq1 <- aurka.exprs <= quantile(qq, probs=0.33, na.rm=TRUE)
qq2 <- aurka.exprs > quantile(qq, probs=0.33, na.rm=TRUE)
qq3 <- aurka.exprs > quantile(qq, probs=0.66, na.rm=TRUE)
qq[qq2] <- 2
qq[qq3] <- 3
qq[qq1] <- 1
mygroup <- qq

##mygroup <- as.numeric(aurka.exprs >= mycutoff)

##pdf("km-plot-test.pdf")
surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
gg <- factor(c(rep("mainz", nrow(pData(mainz.small))), rep("transbig", nrow(pData(transbig.small))), rep("unt", nrow(pData(unt.small))), rep("vdx", nrow(pData(vdx.small))), rep("upp", nrow(pData(upp.small))), rep("nki", nrow(pData(nki.small)))), levels=c("mainz", "transbig", "unt", "vdx", "upp", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Survival", main.title="", sub.title=NULL, leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .col=c("darkblue", "darkred", "orange"), .lty=1, show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE)
title(paste("km survival curve, AURKA"))
##dev.off()




vegfGs <- "VEGF"
vegfGid <- 7422
vegfPaf <- "211527_x_at"
vegfPagi <- "NM_003376"





#######################
## FIGURES FOR APPLICATION NOTE
datasetList <- c("MAINZ","TRANSBIG","UPP","UNT","VDX","NKI","","Overall")

#######################
## CINDEX: AURKA and different datasets as forestplot
#######################
tt <- rbind(cindexall.mainz.small[3,],
            cindexall.transbig.small[3,],
            cindexall.upp.small[3,],
            cindexall.unt.small[3,],
            cindexall.vdx.small[3,],
            cindexall.nki.small[3,],
            NA,
            as.numeric(ccindex[3,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
labeltext <- cbind(c("Dataset",datasetList),c(rep(myspace,length(gsList)+2)))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.8,0.05), xlab=paste("Concordance Index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(0.5,1), is.summary=(c(rep(FALSE,8),TRUE)))
##title(paste("cindex forestplot, AURKA"))
dev.copy2eps(file="figure1a.eps")



#######################
## DINDEX: AURKA and different datasets as forestplot
#######################
tt <- rbind(dindexall.mainz.small[3,],
            dindexall.transbig.small[3,],
            dindexall.upp.small[3,],
            dindexall.unt.small[3,],
            dindexall.vdx.small[3,],
            dindexall.nki.small[3,],
            NA,
            as.numeric(cdindex[3,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
tt <- log2(tt)
labeltext <- cbind(c("Dataset",datasetList))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,2,0.5), xlab=paste("log2 D Index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.5,2), is.summary=(c(rep(FALSE,8),TRUE)))
##title(paste("log2 dindex forestplot, AURKA"))
dev.copy2eps(file="figure1b.eps")



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
            as.numeric(chratio[3,]))

rownames(tt) <- datasetList
tt <- as.data.frame(tt)
tt <- log2(tt)
labeltext <- cbind(c("Dataset",datasetList))
bs <- rep(0.5, nrow(labeltext))
  
forestplot.surv(labeltext=labeltext, mean=c(NA,tt$hratio), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,3.5,0.5), xlab=paste("log2 Hazard Ratio", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.75,3.5),is.summary=(c(rep(FALSE,8),TRUE)))
##title(paste("log2 hratio forestplot, AURKA"))
dev.copy2eps(file="figure1c.eps")



#######################
## censor.time,  km.coxph.plot
#######################

#quantile for cutoff
#34% in the wang2005geneexpression publication
#41% in the foekens2006multicenter publication
#28% in the desmedt2007strong publication
##myprobs <- c(0.33, 0.66)
##mycutoff <- quantile(aurka.exprs, probs=myprobs, na.rm=TRUE)
##mycutoff <- cut2(aurka.exprs, myprobs)
##cut2(x=aurka.exprs, cuts=quantile(aurka.exprs, myprobs), g=3)


#gene symbol: AURKA
#entrezgene id: 6790
#VDX/TRANSBIG/UPP/UNT/MAINSZ probe id: 208079_s_at
#NKI probe id: NM_003600
aurkaGs <- "AURKA"
aurkaGid <- 6790
aurkaPaf <- "208079_s_at"
aurkaPagi <- "NM_003600"


##oncotype test
#aurkaGs <- "MYBL2"
#aurkaGid <- "geneid.4320"
#aurkaPaf <- "203876_s_at"
#aurkaPagi <- "NM_005940"



surv.time.all <- c(pData(mainzSample)[ ,"t.dmfs"], pData(transbigSample)[ ,"t.dmfs"], pData(untSample)[ ,"t.dmfs"], pData(uppSample)[ ,"t.rfs"], pData(vdxSample)[ ,"t.dmfs"], pData(nkiSample)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainzSample)[ ,"e.dmfs"], pData(transbigSample)[ ,"e.dmfs"], pData(untSample)[ ,"e.dmfs"], pData(uppSample)[ ,"e.rfs"], pData(vdxSample)[ ,"e.dmfs"], pData(nkiSample)[ ,"e.dmfs"])            
aurka.exprs <- c(exprs(mainzSample)[aurkaPaf,], exprs(transbigSample)[aurkaPaf,], exprs(untSample)[aurkaPaf,], exprs(uppSample)[aurkaPaf,], exprs(vdxSample)[aurkaPaf,], exprs(nkiSample)[aurkaPagi,])
aurka.exprs.length <- c(length(exprs(mainzSample)[aurkaPaf,]), length(exprs(transbigSample)[aurkaPaf,]), length(exprs(untSample)[aurkaPaf,]), length(exprs(uppSample)[aurkaPaf,]), length(exprs(vdxSample)[aurkaPaf,]), length(exprs(nkiSample)[aurkaPagi,]))

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

##mygroup <- as.numeric(aurka.exprs >= mycutoff)

##pdf("km-plot-test.pdf")
surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
gg <- factor(c(rep("mainz", nrow(pData(mainzSample))), rep("transbig", nrow(pData(transbigSample))), rep("unt", nrow(pData(untSample))), rep("upp", nrow(pData(uppSample))), rep("vdx", nrow(pData(vdxSample))), rep("nki", nrow(pData(nkiSample)))), levels=c("mainz", "transbig", "unt", "upp", "vdx", "nki"))
km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", ylim=c(0.4,1.0), x.label="Time (years)", y.label="Probability of DMFS/RFS", main.title="", sub.title=NULL, leg.text=c("Low   ", "Intermediate   ", "High   "), leg.pos="bottomright", leg.inset=0.05, o.text=NULL, v.line=NULL, h.line=NULL, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=TRUE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE)
##title(paste("survival curve for MYBL2 (4320, 203876_s_at)"))
dev.copy2eps(file="figure1d.eps")

dev.off()
