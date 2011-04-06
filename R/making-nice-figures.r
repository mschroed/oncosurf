#######################
## get packages and dependcies
#######################
#install.packages("../../git/survcomp_1.1.8.tar.gz")
library(survcomp)
library(genefu)
library(Biobase)
library(Hmisc)


#######################
## load eSets
#######################
library("breastCancerNKI")
library("breastCancerTRANSBIG")
library("breastCancerUPP")
library("breastCancerUNT")
library("breastCancerVDX")
library("breastCancerMAINZ")

data(mainz,transbig,upp,unt,vdx,nki)


#######################
## concordance.index computation for complete datasets
#######################
system.time(cindexall.mainz <- t(apply(X=exprs(mainz), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainz)[ ,"t.dmfs"], z=pData(mainz)[ ,"e.dmfs"])))

system.time(cindexall.transbig <- t(apply(X=exprs(transbig), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbig)[ ,"t.dmfs"], z=pData(transbig)[ ,"e.dmfs"])))

system.time(cindexall.vdx <- t(apply(X=exprs(vdx), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdx)[ ,"t.dmfs"], z=pData(vdx)[ ,"e.dmfs"])))

system.time(cindexall.upp <- t(apply(X=exprs(upp), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(upp)[ ,"t.rfs"], z=pData(upp)[ ,"e.rfs"])))

system.time(cindexall.unt <- t(apply(X=exprs(unt), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(unt)[ ,"t.dmfs"], z=pData(unt)[ ,"e.dmfs"])))

system.time(cindexall.nki <- t(apply(X=exprs(nki), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nki)[ ,"t.dmfs"], z=pData(nki)[ ,"e.dmfs"])))


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


#######################
## top low/high risk genes NKI
#######################
qm <- quantile(cindexall.nki[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE)
##cindexall.nki[cindexall.nki[,"cindex"] > qm[[1]],]
low.risk.probes.nki <- as.data.frame(cindexall.nki[cindexall.nki[,"cindex"] > qm[[2]],])
low.risk.probes.nki <- low.risk.probes.nki[order(low.risk.probes.nki$cindex, decreasing=TRUE),]
high.risk.probes.nki <- as.data.frame(cindexall.nki[cindexall.nki[,"cindex"] < qm[[1]],])
high.risk.probes.nki <- high.risk.probes.nki[order(high.risk.probes.nki$cindex),]

## TODO
low.risk.probes.nki$gs <- fData(nki)[fData(nki)[rownames(low.risk.probes.nki),"probe"],"Gene.symbol"]
high.risk.probes.nki$gs <- fData(nki)[fData(nki)[rownames(high.risk.probes.nki),"probe"],"Gene.symbol"]



#######################
## VENN DIAGRAM sample data
#######################
gLowMainz <- as.numeric(low.risk.probes.mainz[1:1000,"gid"])
gLowTransbig <- as.numeric(low.risk.probes.transbig[1:1000,"gid"])
gLowUpp <- as.numeric(low.risk.probes.upp[1:1000,"gid"])
g <- cbind(gLowMainz,gLowTransbig,gLowUpp)
c <- vennCounts(g, include="both")
vennDiagram(c)





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









#######################
## some plots
#######################
pdf("cindex-histogram-all-datasets.pdf")
## all datasets
hist(c(cindexall.mainz[,"cindex"], cindexall.transbig[,"cindex"], cindexall.upp[,"cindex"], cindexall.unt[,"cindex"], cindexall.vdx[,"cindex"], cindexall.nki[,"cindex"]), breaks="FD", main="Concordance Indices from 6 Breast Cancer Datasets", col=c("darkblue"), xlab="Concordance Index")
abline(v=quantile(c(cindexall.mainz[,"cindex"], cindexall.transbig[,"cindex"], cindexall.upp[,"cindex"], cindexall.unt[,"cindex"], cindexall.vdx[,"cindex"], cindexall.nki[,"cindex"]), probs=c(0.05,0.95),na.rm=TRUE))


den.mainz <- density(cindexall.mainz[,"cindex"])
den.transbig <- density(cindexall.transbig[,"cindex"])
den.upp <- density(cindexall.upp[,"cindex"])
den.unt <- density(cindexall.unt[,"cindex"])
den.vdx <- density(cindexall.vdx[,"cindex"])
den.nki <- density(cindexall.nki[,"cindex"])

plot(den.vdx, col="red", xlab="Concordance Index", main="Concordance Indices from 6 Breast Cancer Datasets")
lines(den.transbig, col="darkgreen")
lines(den.upp, col="darkblue")
lines(den.unt, col="brown")
lines(den.mainz, col="darkorange")
lines(den.nki, col="black")
legend("topright", c("MAINZ","TRANSBIG","UPP","UNT","VDX","NKI"), col=c("darkorange","darkgreen","darkblue","brown","red","black"), pch=-1, lty=1)
abline(v=0.5, col="red", lty=3)


## single datasets
hist(cindexall.mainz[,"cindex"], breaks = "FD", main="Concordance Indices from MAINZ", col=c("darkblue"), xlab="Concordance Index")
abline(v=quantile(cindexall.mainz[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE))

hist(cindexall.transbig[,"cindex"], breaks = "FD", main="Concordance Indices from TRANSBIG", col=c("darkgreen"), xlab="Concordance Index")
abline(v=quantile(cindexall.transbig[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE))

hist(cindexall.upp[,"cindex"], breaks = "FD", main="Concordance Indices from UPP", col=c("red"), xlab="Concordance Index")
abline(v=quantile(cindexall.upp[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE))

hist(cindexall.unt[,"cindex"], breaks = "FD", main="Concordance Indices from UNT", col=c("brown"), xlab="Concordance Index")
abline(v=quantile(cindexall.unt[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE))

hist(cindexall.vdx[,"cindex"], breaks = "FD", main="Concordance Indices from VDX", col=c("darkorange"), xlab="Concordance Index")
abline(v=quantile(cindexall.vdx[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE))

hist(cindexall.nki[,"cindex"], breaks = "FD", main="Concordance Indices from NKI", col=c("black"), xlab="Concordance Index")
abline(v=quantile(cindexall.nki[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE))

##plot(c(quantile(cindexall.mainz[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE), quantile(cindexall.transbig[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE), quantile(cindexall.upp[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE), quantile(cindexall.unt[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE), quantile(cindexall.vdx[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE), quantile(cindexall.nki[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE)))

boxplot(quantile(cindexall.mainz[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE), quantile(cindexall.transbig[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE), quantile(cindexall.upp[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE), quantile(cindexall.unt[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE), quantile(cindexall.vdx[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE), quantile(cindexall.nki[,"cindex"], probs=c(0.05,0.95),na.rm=TRUE), main="Cindex Top/Bottom 5% Quantiles", xlab="MAINZ       TRANSBIG         UPP             UNT             VDX             NKI", ylab="Concordance Index")

dev.off()









#######################
## concordance.index computation for SAMPLE DATASETS
#######################
system.time(cindexSample.mainz <- t(apply(X=exprs(mainzSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(mainzSample)[ ,"t.dmfs"], z=pData(mainzSample)[ ,"e.dmfs"])))

system.time(cindexSample.transbig <- t(apply(X=exprs(transbigSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(transbigSample)[ ,"t.dmfs"], z=pData(transbigSample)[ ,"e.dmfs"])))

system.time(cindexSample.vdx <- t(apply(X=exprs(vdxSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(vdxSample)[ ,"t.dmfs"], z=pData(vdxSample)[ ,"e.dmfs"])))

system.time(cindexSample.upp <- t(apply(X=exprs(uppSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(uppSample)[ ,"t.rfs"], z=pData(uppSample)[ ,"e.rfs"])))

system.time(cindexSample.unt <- t(apply(X=exprs(untSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(untSample)[ ,"t.dmfs"], z=pData(untSample)[ ,"e.dmfs"])))

system.time(cindexSample.nki <- t(apply(X=exprs(nkiSample), MARGIN=1, function(x, y, z) { tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE); return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "lower"=tt$lower, "upper"=tt$upper)); }, y=pData(nkiSample)[ ,"t.dmfs"], z=pData(nkiSample)[ ,"e.dmfs"])))




allData <- rbind(cindexSample.mainz[,"cindex"], cindexSample.transbig[,"cindex"], cindexSample.upp[,"cindex"], cindexSample.unt[,"cindex"], cindexSample.vdx[,"cindex"], cindexSample.nki[,"cindex"])
boxplot(allData)


pdf("making-nice-figures.pdf")


plot(cindexSample.mainz[,"cindex"],type="b",col="red", axes=FALSE, xlim=c(0.5,7.5), ylim=c(0.3,0.7), xlab="Gene Symbol", ylab="Concordance Index", main="Concordance Indices from all Datasets plotted for each Gene", lwd=2)
lines(cindexSample.transbig[,"cindex"],type="b",col="blue", lwd=2)
lines(cindexSample.upp[,"cindex"],type="b",col="lightblue", lwd=2)
lines(cindexSample.unt[,"cindex"],type="b",col="green", lwd=2)
lines(cindexSample.vdx[,"cindex"],type="b",col="black", lwd=2)
lines(cindexSample.nki[,"cindex"],type="b",col="yellow", lwd=2)
box()
axis(1, at=1:7, lab=c(fData(mainzSample)[,"Gene.symbol"]))
axis(2, las=1, at=0.1*3:7)
abline(v=1:length(cindexSample.mainz), lty=3)
abline(h=c(0.45,0.55), lty=2, col="blue")
abline(h=0.5, lty=3, col="darkred")
legend("bottomright", c("mainz","transbig","upp","unt","vdx","nki"), col=c("red","blue","lightblue","green","black","yellow"), lty=1, lwd=2, bg="white")
boxplot(allData, add=TRUE, axes=FALSE, ylim=FALSE, xlim=FALSE, border="grey")


plot(cindexSample.mainz[,"cindex"],type="p",col="red", axes=FALSE, xlim=c(0.5,7.5), ylim=c(0.3,0.7), xlab="Gene Symbol", ylab="Concordance Index", main="Concordance Indices from all Datasets plotted for each Gene", lwd=5)
lines(cindexSample.transbig[,"cindex"],type="p",col="blue", lwd=5)
lines(cindexSample.upp[,"cindex"],type="p",col="lightblue", lwd=5)
lines(cindexSample.unt[,"cindex"],type="p",col="green", lwd=5)
lines(cindexSample.vdx[,"cindex"],type="p",col="black", lwd=5)
lines(cindexSample.nki[,"cindex"],type="p",col="yellow", lwd=5)
box()
axis(1, at=1:7, lab=c(fData(mainzSample)[,"Gene.symbol"]))
axis(2, las=1, at=0.1*3:7)
abline(v=1:length(cindexSample.mainz), lty=3)
abline(h=c(0.45,0.55), lty=2, col="blue")
abline(h=0.5, lty=3, col="darkred")
legend("bottomright", c("mainz","transbig","upp","unt","vdx","nki"), col=c("red","blue","lightblue","green","black","yellow"), pch=1, pt.lwd=5, bg="white")
boxplot(allData, add=TRUE, axes=FALSE, ylim=FALSE, xlim=FALSE, border="grey")

dev.off()









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
labeltext <- cbind(c("Dataset",dataset.list),c(rep(mybigspace,length(dataset.list))))
bs <- rep(0.5, nrow(labeltext))                              

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.8,0.05), xlab=paste("AURKA cindex", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(0.5,1), is.summary=(c(rep(FALSE,7),TRUE)))
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

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.75,0.05), xlab=paste("VEGF cindex", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(0.4,0.75), is.summary=(c(rep(FALSE,7),TRUE)))
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

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$cindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0.5, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(0.4,0.8,0.05), xlab=paste("AURKA and VEGF cindex", myspace, sep=""), col=meta.colors(line=c(rep(c(NA,"darkblue","darkred"),7)),zero="darkred",box=c(rep(c(NA,"royalblue","red"),7))), box.size=bs, clip=c(0.4,1), is.summary=(c(rep(FALSE,19),TRUE,TRUE)))
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

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,2,0.5), xlab=paste("AURKA log2 D.index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.25,2), is.summary=(c(rep(FALSE,7),TRUE)))
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

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,1.5,0.5), xlab=paste("VEGF log2 D.index", myspace, sep=""), col=meta.colors(box="royalblue",line="darkblue",zero="darkred"), box.size=bs, clip=c(-0.5,1.5), is.summary=(c(rep(FALSE,7),TRUE)))
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

forestplot.surv(labeltext=labeltext, mean=c(NA,tt$dindex), lower=c(NA,tt$lower), upper=c(NA,tt$upper), zero=0, align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(-0.5,2,0.5), xlab=paste("AURKA and VEGF log2 D.index", myspace, sep=""), col=meta.colors(line=c(rep(c(NA,"darkblue","darkred"),7)),zero="darkred",box=c(rep(c(NA,"royalblue","red"),7))), box.size=bs, clip=c(-0.5,2), is.summary=(c(rep(FALSE,19),TRUE,TRUE)))
title(paste("log2 D.index forestplot, AURKA and VEGF"))    

                         
                         
                         
                         
                         
##########################
## expression plot sorted by gene expression value
##########################                        
plot(sort(exprs(mainzSample)["208079_s_at",]))
plot(sort(c(exprs(mainzSample)["208079_s_at",],exprs(transbigSample)["208079_s_at",],exprs(uppSample)["208079_s_at",],exprs(untSample)["208079_s_at",],exprs(vdxSample)["208079_s_at",],exprs(nkiSample)["NM_003600",])))

plot(sort(c(exprs(mainzSample)["208079_s_at",],exprs(transbigSample)["208079_s_at",],exprs(uppSample)["208079_s_at",],exprs(untSample)["208079_s_at",],exprs(vdxSample)["208079_s_at",],exprs(nkiSample)["NM_003600",])))
 