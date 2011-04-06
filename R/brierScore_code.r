####################################
## Brier score wrt the time

#library(survcomp)
library(Biobase)
tc <- 10 * 365  

gsList <- c("ESR1", "ERBB2", "AURKA", "PLAU", "VEGF", "STAT1", "CASP3")
probesAffy <- c("205225_at", "216836_s_at", "208079_s_at", "211668_s_at", "211527_x_at", "209969_s_at", "202763_at")



surv.time.all <- c(pData(mainz7g)[ ,"t.dmfs"], pData(transbig7g)[ ,"t.dmfs"])             
surv.event.all <- c(pData(mainz7g)[ ,"e.dmfs"], pData(transbig7g)[ ,"e.dmfs"])
mygene <- t(cbind(exprs(mainz7g), exprs(transbig7g)))
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
		plot(myperf[[i]]$time, myperf[[i]]$bsc, typ="l", xlab="Time (years)", ylab="Brier score", col=mycol[i], lty=mylty[i], ylim=c(0,0.4), lwd=mylwd[i], main="Brier Score")
	} else {  lines(myperf[[i]]$time, myperf[[i]]$bsc, col=mycol[i], lty=mylty[i], lwd=mylwd[i]) }
	mysbrier.int <- unlist(lapply(myperf, function(x) { return(x$bsc.integrated)}))
}
smartlegend(x="left", y="top", legend=sprintf(names(myperf), mysbrier.int), col=mycol, lty=mylty, lwd=mylwd)

