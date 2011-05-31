require("rjson")
require("RCurl")
require("Biobase")

.onLoad <- function(lib, pkg) {
    ##options(GXA_API_URL = "http://lime.ebi.ac.uk:14068/gxa/api/v2")
    options(GXA_API_URL = "http://localhost:9999/gxa/api/v2")
}

.onLoad()

getGXAUrl <- function(section) {
    paste(getOption("GXA_API_URL"), section, sep = "/")
}

getFromServer <- function(domain, request) {
    response <- fromJSON(getURL(url = getGXAUrl(domain), post = TRUE, postfields=toJSON(request)))
    if (!is.null(response$errorText)) {
        stop(response$errorText)
    }
    response
}

listExperimentsWithFactor <- function(factorName) {
    query <- list("hasFactor" = list(list("name" = factorName)))
    filter <- list("experiments" = list("accession" = "FIELD"))
    request <- list(query = query, filter = filter)
    response <- getFromServer("experiments", request)
    sapply(response$experiments, function(experiment) { experiment$accession })
}

listAssays <- function(experimentAccession){
    query <- list("experimentAccession" = experimentAccession)
    filter <- list("*" = "ALL")
    request <- list(query = query, filter = filter)
    assays <- getFromServer("assays", request)
    names(assays) <- sapply(assays, function(a) { a$accession })
    assays
}

loadExpressionSetsForAssayAccessions <- function(experimentAccession, assayAccessions, geneNames = NULL, geneIdentifiers = NULL) {
    if (is.null(geneNames) && is.null(geneIdentifiers)) {
        stop("You must specify geneNames OR geneIdentifiers"); 
    }
    if (!is.null(geneNames) && !is.null(geneIdentifiers)) {
        stop("You cannot specify both geneNames AND geneIdentifiers"); 
    }

    # Obtaining assay data from server, filtering for accessions
    assays <- listAssays(experimentAccession)
    if(assayAccessions != "*"){assays <- assays[assayAccessions]}
    accessionsNotFound <- assayAccessions[is.na(names(assays))]
    if (length(accessionsNotFound) > 0) {
        print("Assays not found for accessions:"); 
        print(accessionsNotFound)
        assays <- assays[!is.na(names(assays))]
    }
    loadExpressionSetsForAssays(experimentAccession, assays, geneNames = geneNames, geneIdentifiers = geneIdentifiers)
}


loadExpressionSetsForAssays <- function(experimentAccession, assays, geneNames = NULL, geneIdentifiers = NULL) {
    arrayDesignAccessions <- sapply(assays, function(a) { a$arrayDesignAccession })
    sapply(unique(arrayDesignAccessions), function(ad) {
        loadExpressionSetForArrayDesignAndAssays(
            experimentAccession, ad, assays[arrayDesignAccessions == ad],
            geneNames = geneNames,
            geneIdentifiers = geneIdentifiers
        )
    })
}

loadExpressionSetForArrayDesignAndAssays <- function(experimentAccession, arrayDesignAccession, assays, geneNames = NULL, geneIdentifiers = NULL) {
    # Requesting expression data from server
    if (!is.null(geneNames)) {
        query <- list("experimentAccession" = experimentAccession, "assayAccessions" = names(assays), "geneNames" = geneNames)
    } else {
        query <- list("experimentAccession" = experimentAccession, "assayAccessions" = names(assays), "geneIdentifiers" = geneIdentifiers)
    }
    filter <- list("*" = "ALL")
    request <- list(query = query, filter = filter)
    response <- getFromServer("data", request)

    # Creating assays data (aData)
    as <- response[[1]]$assayAccessions
    de <- sapply(response[[1]]$genes, function(x) { x$designElementAccession } )
    bdc <- t(sapply(response[[1]]$genes, function(x) { x$expressionLevels } ))
    colnames(bdc) <- as
    rownames(bdc) <- de
    aData <- assayDataNew(storage.mode = "lockedEnvironment", exprs = bdc)
    sampleNames(aData) <- as
    featureNames(aData) <- de

    # Creating pheno data (pData)
    ef <- unique(c(sapply(assays, function(a) { names(a$properties) })))
    efOrder <- match(as,names(assays))
    efv <- data.frame(t(sapply(assays[efOrder], function(a) { sapply(a$properties[ef], function(p) p[[1]]$value)})))
    pData = new("AnnotatedDataFrame", data = efv)
    sampleNames(pData) <- as
    
    if (!is.null(geneNames)) {
        gn <- sapply(response[[1]]$genes, function(x) { x$geneName } )
    } else {
        gn <- sapply(response[[1]]$genes, function(x) { x$geneIdentifier } )
    }
    gn <- sapply(response[[1]]$genes, function(x) { x$geneName } )
    fDataFrame <- data.frame(de = de, deacc = de, gn = gn)
    fData = new("AnnotatedDataFrame", data = fDataFrame)
    featureNames(fData) = de
    
    # Creating experiment data (eData)
    eData = new("MIAME",
        other = list(
            accession = experimentAccession,
            experimentid = experimentAccession,
            arraydesignid = arrayDesignAccession,
            qt = 0, # TODO: what is this? (quantitation type)
            arraydesignaccnum = arrayDesignAccession,
            arraydesignname = arrayDesignAccession # TODO: have to use description instead
        )
    )

    # Returning an ExpressionSet
    new("ExpressionSet", assayData = aData, phenoData = pData, featureData = fData, experimentData = eData)
}

listGenesByNames <- function(geneNames) {
    query <- list("geneNames" = geneNames)
    filter <- list("*" = "ALL")
    request <- list(query = query, filter = filter)
    response <- getFromServer("genes", request)
    response
}

listGenesByIdentifiers <- function(geneIdentifiers) {
    query <- list("geneIdentifiers" = geneIdentifiers)
    filter <- list("*" = "ALL")
    request <- list(query = query, filter = filter)
    response <- getFromServer("genes", request)
    response
}


## precompute statistic
precomputeResults <- function(){
  timeValue <- proc.time()
  timeHistory <- NULL
  timeHistory <- as.data.frame(t(c("start",as.numeric(proc.time() - proc.time()))))
  library(genefu)
  library(survcomp)
  
  experimentAccession <- c("E-MAIN-1234", "E-TRAN-1234", "E-UPPP-1234", "E-UNTT-1234", "E-VDXX-1234")
  subtypesOrder <- c("ER-/HER2-","HER2+","ER+/HER2-") 
  subtypesOrderAll <- c("ER-/HER2-","HER2+","ER+/HER2-","All Subtypes")
  
  gxaData <- as.list(NULL)
  gxaESets <- as.list(NULL)
  
  timeHistory <- rbind(timeHistory, t(c("start loop",as.numeric(proc.time() - timeValue))))
  
  for(expAcc in 1:length(experimentAccession)){
    assays <- listAssays(experimentAccession[expAcc])
    assayAccessions <- "*"  ## c("UPP_103B41", "UPP_104B91", "UPP_280C43")
    geneNames <- "*" ## c("AURKA", "ERBB2", "ESR1", "STAT1", "CASP3", "VEGFA", "PLAU")
    E <- loadExpressionSetsForAssayAccessions(experimentAccession[expAcc], assayAccessions, geneNames)
    E <- E[[1]]
    E
    timeHistory <- rbind(timeHistory, t(c("loaded eSet",as.numeric(proc.time() - timeValue))))
    ## gxaData <- array(NA, dim=c(length(unionGeneList), 2, length(data.all), length(subtypesOrder)), dimnames=list(unionGeneListN, c("cindex", "cindex.se"), names(data.all), subtypesOrder))
    
    ##sampleData <- as.data.frame(t(sapply(pData(E)$properties,function(x)sapply(x,function(y)y[[1]]$value))))
    sampleData <- pData(E)
    
    ###################################################
    ### get subtypes
    ###################################################
    subtypeData <- NULL
    domap <- FALSE
    subtypeData <- subtype.cluster.predict(sbt.model=scmgene,
      data=t(exprs(E)), annot=fData(E), do.mapping=domap, do.scale=TRUE,
      do.prediction.strength=FALSE, do.BIC=FALSE, plot=FALSE, verbose=TRUE)
    
    timeHistory <- rbind(timeHistory, t(c("computed subtypes",as.numeric(proc.time() - timeValue))))
    
    ###################################################
    ### data
    ###################################################  
    expData <- array(NA, dim=c(nrow(exprs(E)), 3, 1, 3, length(subtypesOrder)+1), dimnames=list(rownames(exprs(E)), c("est", "se", "p.value"), experimentAccession[expAcc], c("cindices","dindices","hratios"), subtypesOrderAll))
  
    exprsDataTmp <- exprs(E)
    if(expAcc != 3){
      timeDataTmp <- as.numeric(as.character(sampleData[,"t.dmfs"]))
      eventDataTmp <- as.numeric(as.character(sampleData[,"e.dmfs"]))
    } else {
      timeDataTmp <- as.numeric(as.character(sampleData[,"t.rfs"]))
      eventDataTmp <- as.numeric(as.character(sampleData[,"e.rfs"]))
    }
    sbtprob <- subtypeData$subtype.proba
    
    timeHistory <- rbind(timeHistory, t(c("data",as.numeric(proc.time() - timeValue))))
    
    ###################################################
    ### compute concordance indices
    ###################################################
    for(j in 1:length(subtypesOrder)){
      subtypeWeights <- as.numeric(sbtprob[,j])    
      cindexTmp <- t(apply(X=exprsDataTmp, MARGIN=1, FUN=function(w, x, y, z) {
        tt <- concordance.index(x=x, surv.time=y, surv.event=z, weights=w, method="noether", na.rm=TRUE);
        return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "p.value"=tt$p.value)); },
        y=timeDataTmp, z=eventDataTmp, w=subtypeWeights))
      expData[rownames(exprsDataTmp), , , "cindices", subtypesOrder[j]] <- cindexTmp
    }
    cindexTmp <- t(apply(X=exprsDataTmp, MARGIN=1, FUN=function(x, y, z) {
        tt <- concordance.index(x=x, surv.time=y, surv.event=z, method="noether", na.rm=TRUE);
        return(c("cindex"=tt$c.index, "cindex.se"=tt$se, "p.value"=tt$p.value)); },
        y=timeDataTmp, z=eventDataTmp))
    expData[rownames(exprsDataTmp), , , "cindices", "All Subtypes"] <- cindexTmp
    
    timeHistory <- rbind(timeHistory, t(c("cindex",as.numeric(proc.time() - timeValue))))
    
    ###################################################
    ### compute D indices
    ###################################################
    for(j in 1:length(subtypesOrder)){
      subtypeWeights <- as.numeric(sbtprob[,j])    
      dindexTmp <- t(apply(X=exprsDataTmp, MARGIN=1, FUN=function(w, x, y, z) {
        tt <- D.index(x=x, surv.time=y, surv.event=z, weights=w, na.rm=TRUE);
        return(c("dindex"=tt$d.index, "dindex.se"=tt$se, "p.value"=tt$p.value)); },
        y=timeDataTmp, z=eventDataTmp, w=subtypeWeights))
      expData[rownames(exprsDataTmp), , , "dindices", subtypesOrder[j]] <- dindexTmp
    }
    dindexTmp <- t(apply(X=exprsDataTmp, MARGIN=1, FUN=function(x, y, z) {
        tt <- D.index(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
        return(c("dindex"=tt$coef, "dindex.se"=tt$se, "p.value"=tt$p.value)); },
        y=timeDataTmp, z=eventDataTmp))
    expData[rownames(exprsDataTmp), , , "dindices", "All Subtypes"] <- dindexTmp
    
    timeHistory <- rbind(timeHistory, t(c("dindex",as.numeric(proc.time() - timeValue))))
    
    ###################################################
    ### rescale function
    ###################################################
    rescale <- function(x, na.rm=FALSE, q=0.05){
      ma <- quantile(x, probs=1-(q/2), na.rm=na.rm)
      mi <- quantile(x, probs=q/2, na.rm=na.rm)
      x <- (x - mi) / (ma - mi)
      return((x - 0.5) * 4)
    }
    
    ###################################################
    ### compute hazard ratio
    ###################################################
    for(j in 1:length(subtypesOrder)){
      subtypeWeights <- as.numeric(sbtprob[,j])                            
      hratioTmp <- t(apply(X=rescale(exprsDataTmp , q=0.05, na.rm=TRUE), MARGIN=1, FUN=function(w, x, y, z) {
        tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, weights=w, na.rm=TRUE);
        return(c("hratio"=tt$hazard.ratio, "hratio.se"=tt$se, "p.value"=tt$p.value)); },
        y=timeDataTmp, z=eventDataTmp, w=subtypeWeights))
      expData[rownames(exprsDataTmp), , , "hratios", subtypesOrder[j]] <- hratioTmp
    }
    hratioTmp <- t(apply(X=rescale(exprsDataTmp , q=0.05, na.rm=TRUE), MARGIN=1, FUN=function(x, y, z) {
      tt <- hazard.ratio(x=x, surv.time=y, surv.event=z, na.rm=TRUE);
      return(c("hratio"=tt$coef, "hratio.se"=tt$se, "p.value"=tt$p.value)); },
      y=timeDataTmp, z=eventDataTmp))
    expData[rownames(exprsDataTmp), , , "hratios", "All Subtypes"] <- hratioTmp
    
    timeHistory <- rbind(timeHistory, t(c("hratio",as.numeric(proc.time() - timeValue))))
    
    ###################################################
    ### store data
    ###################################################
    gxaData <- c(gxaData, list(list(expData, subtypeData)))
  }
  names(gxaData) <- experimentAccession
  save(gxaData,file="gxaData_bu.RData")
  save(timeHistory,file="timeHistory.RData")
  return(timeHistory)
}


###################################################
### draw boxplots for each gene in each dataset
###################################################
gxa.boxplot <- function(E = NULL, probe = NULL, subtypeData = NULL, subtypesOrderAll = NULL){
  if(is.null(E) || is.null(probe) || is.null(subtypeData) || is.null(subtypesOrderAll)){
    stop("gxa.boxplot: You must specifiy eSet, probe, subtype data and subtype order!")
  }
  ##pdf("boxplot-test.pdf")
  selectProbe <- fData(E)$deacc == probe
  gs <- as.character(fData(E)$gn[selectProbe])
  patientsSubType <- dd <- gg <- as.list(NULL)
  ddMin <- ddMax <- NULL
  exprsData <- exprs(E)
  patientsSubType <- as.character(subtypeData$subtype)
  for(j in 1:length(subtypesOrderAll)){
    if(j != length(subtypesOrderAll)){
      gg[[j]] <- patientsSubType == subtypesOrder[j]
      dd[[j]] <- exprsData[probe,gg[[j]]]
    } else {
      dd[[j]] <- exprsData[probe,]
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
    main=paste(gs," ", probe), xlab="", yaxt="n", xaxt="n")
  labels <- subtypesOrderAll
  text(1:4, par("usr")[3] - 0.05, srt = 35, adj = 1,
    labels = labels, xpd = TRUE)
  ##dev.off()
}


###################################################
### getting the survival type
###################################################
getSurvivalType <- function(pDataEset = NULL){
  if(is.null(pDataEset)){
    stop("getSurvivalType: You must specify cov!")
  }
  cova <- colnames(pDataEset)
  cova <- cova[grep(pattern="^t\\.",x=cova)]
  covaNA <- NULL
  for(i in 1:length(cova)){
    covaTmp <- as.character(pDataEset[,cova[i]])
    covaTmp <- covaTmp[covaTmp != "NA"]
    covaNA <- rbind(covaNA, c(length(covaTmp),cova[i]))
  }
  sType <- covaNA[covaNA[,1] == max(covaNA[,1]),2]
  sType <- c(sType,sprintf("%s.%s","e",strsplit(sType,"\\.")[[1]][2]))
  return(sType)
}


###################################################
### integrated Brier Score depending on time
###################################################
gxa.sbrier.plot <- function(E = NULL, probe = NULL, subtypeData = NULL, subtypesOrderAll = NULL) {
  if(is.null(E) || is.null(probe) || is.null(subtypeData) || is.null(subtypesOrderAll)){
    stop("gxa.sbrier.plot: You must specify eSet, probe, subtype data and subtype order!")
  }
  ##pdf("sbier-test.pdf")  
  sType <- getSurvivalType(pData(E))   
  surv.time.all <- as.numeric(as.character(pData(E)[ ,sType[1]]))
  surv.event.all <- as.numeric(as.character(pData(E)[ ,sType[2]]))
  mygene <- cbind(exprs(E)[rownames(exprs(E))==probe,])
  
  tc <- 10 * 365
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
  
  ## for probe and subtypes
  ddtt <- cbind(dd.ts[ , c("time", "event", "mygene")])
  colnames(ddtt) <- c("time", "event", "score")
  
  patientsSubType <- as.character(subtypeData$subtype)
  for(j in 1:length(subtypesOrderAll)){
    if(j != length(subtypesOrderAll)){
      gg <- patientsSubType == subtypesOrder[j]
    } else {
      gg <- rep(TRUE,nrow(ddtt))
    }
    res.ts <- sbrier.score2proba(data.tr=ddtt[gg,], data.ts=ddtt[gg,], method="cox")    
    sbrier.score2proba.ts.all <- c(sbrier.score2proba.ts.all, list(res.ts))
    names(sbrier.score2proba.ts.all)[length(sbrier.score2proba.ts.all)] <- subtypesOrderAll[j]
  }
  myperf <-  sbrier.score2proba.ts.all

  ylimMax <- 0.25
  for(k in 1:length(myperf)){
    maxVal <- max(myperf[[k]]$bsc)
    if(ylimMax < maxVal){
      ylimMax <- maxVal
    }
  }
  
  mycol <- c("#000000FF", rainbow(length(myperf))[-1])
  mylty <- 1:(length(myperf))
  mylwd <- c(3, rep(2, length(myperf)-1))
  for(i in 1:length(myperf)) {
  	if(i == 1) {
  		plot(myperf[[i]]$time, myperf[[i]]$bsc, typ="l", xlab="Time (years)", ylab="Brier score", col=mycol[i], lty=mylty[i], ylim=c(0,ylimMax), lwd=mylwd[i], main=sprintf("Brier Scores for %s", probe))
  	} else {  lines(myperf[[i]]$time, myperf[[i]]$bsc, col=mycol[i], lty=mylty[i], lwd=mylwd[i]) }
  	mysbrier.int <- unlist(lapply(myperf, function(x) { return(x$bsc.integrated)}))
  }
  smartlegend(x="left", y="top", legend=sprintf("%s, IBSC = %.3g", names(myperf), mysbrier.int), col=mycol, lty=mylty, lwd=mylwd)
  ##dev.off()
}


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
### time dependent receiver operating characteristic curve
###################################################
gxa.tdrocc.plot <- function(E = NULL, probe = NULL, subtypeData = NULL, subtypesOrderAll = NULL, timeData = NULL) {
  if(is.null(E) || is.null(probe) || is.null(subtypeData) || is.null(subtypesOrderAll) || is.null(timeData)){
    stop("gxa.tdrocc.plot: You must specify eSet, probe, subtype data, subtype order and time!")
  }
  ##pdf("tdrocc-test.pdf")
  tc <- 10 * 365 
  sType <- getSurvivalType(pData(E))   
  surv.time.all <- as.numeric(as.character(pData(E)[ ,sType[1]]))
  surv.event.all <- as.numeric(as.character(pData(E)[ ,sType[2]]))
  surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
  stime <- surv.data$surv.time.cens
  sevent <- surv.data$surv.event.cens
  mygene <- exprs(E)[rownames(exprs(E))==probe,]
  exprsData <- as.numeric(rescale(mygene, q=0.05, na.rm=TRUE))
  
  tdroc.all <- NULL
  
  patientsSubType <- as.character(subtypeData$subtype)
  for(i in 1:length(subtypesOrderAll)){
    if(i != length(subtypesOrderAll)){
      gg <- patientsSubType == subtypesOrder[i]
    } else {
      gg <- rep(TRUE,length(exprsData))
    }    
    tdroc <- tdrocc(x=exprsData[gg], surv.time=stime[gg], surv.event=sevent[gg], time=timeData, na.rm=TRUE, verbose=FALSE)
    tdroc.all <- c(tdroc.all, list(tdroc))
  }
  names(tdroc.all) <- subtypesOrderAll
  
  ## plot the time-dependent ROC curve
  colorOrder <- mycol <- c("#000000FF", rainbow(length(tdroc.all))[-1])
  plot(x=1-tdroc.all[[1]]$spec, y=tdroc.all[[1]]$sens, type="l", xlab="1 - specificity", ylab="sensitivity", xlim=c(0, 1), ylim=c(0, 1), col=colorOrder[1], main=sprintf("TDROCC for %s", probe))
  lines(x=c(0,1), y=c(0,1), lty=3, col="red")
  for(j in 2:length(subtypesOrderAll)){
    lines(x=1-tdroc.all[[j]]$spec, y=tdroc.all[[j]]$sens, type="l", col=colorOrder[j])
  }
  legend("bottomright", names(tdroc.all), col=colorOrder, lty=1, bg="white")
  ##dev.off()
}


###################################################
### forestplot
###################################################
gxa.forest.plot <- function(estName = NULL, estData = NULL, probe = NULL, subtypeData = NULL, subtypesOrderAll = NULL) {
  if(is.null(estName) || is.null(estData) || is.null(probe) || is.null(subtypeData) || is.null(subtypesOrderAll)){
    stop("gxa.forest.plot: You must specify estName, estData, probe, subtype data and subtype order!")
  }
  ##pdf(sprintf("forest-test-%s.pdf",estName))
  myest <- estData[rownames(estData[,,,estName,1]) == probe,,,estName,]
  tt <- as.data.frame(t(rbind(apply(myest,2,function(x)x[1]))))
  colnames(tt) <- "est"
  tt$se <- t(rbind(apply(myest,2,function(x)x[2])))
  labeltext <- cbind(c("Subtype",subtypesOrderAll),
                     c(rep("   ",length(subtypesOrderAll)+1)))
  bs <- rep(0.5, nrow(labeltext))
  if(estName == "cindices"){
    r.mean <- c(NA,tt$est)
    r.lower <- c(NA, tt$est + qnorm(0.025, lower.tail=TRUE) * tt$se)
    r.upper <- c(NA, tt$est + qnorm(0.025, lower.tail=FALSE) * tt$se)
    xticksLower <- 0.3 ##round(min(r.lower,na.rm=TRUE),digits=2)
    xticksUpper <- 0.8 ##round(max(r.lower,na.rm=TRUE),digits=2)
    xticksStep <- 0.1
  } else {
    r.mean <- c(NA,log2(tt$est))
    r.lower <-  c(NA, log2(tt$est + qnorm(0.025, lower.tail=TRUE) * tt$se))
    r.upper <- c(NA, log2(tt$est + qnorm(0.025, lower.tail=FALSE) * tt$se))
    xticksLower <- -2
    xticksUpper <- 2
    xticksStep <- 0.25
  }
  forestplot.surv(labeltext=labeltext, mean=r.mean, lower=r.lower, upper=r.upper, zero=0.5,
      align=c("l"), graphwidth=unit(2, "inches"), x.ticks=seq(xticksLower,xticksUpper,xticksStep),
      col=meta.colors(line=c(rep(c(NA, rainbow(length(subtypesOrderAll))),length(subtypesOrderAll))), zero="firebrick",
      box=c(rep(c(NA,rainbow(length(subtypesOrderAll))),length(subtypesOrderAll)))), box.size=bs, clip=c(xticksLower,xticksUpper), is.summary=FALSE)
  title(sprintf("Forest plot for %s", probe))
  ##dev.off()
}


###################################################
### Kaplan Meier survival curve
###################################################
gxa.km.plot <- function(E = NULL, probe = NULL, subtypeData = NULL, subtypesOrderAll = NULL) {
  if(is.null(E) || is.null(probe) || is.null(subtypeData) || is.null(subtypesOrderAll)){
    stop("gxa.tdrocc.plot: You must specify eSet, probe, subtype data and subtype order!")
  }
  ##pdf("km-test.pdf")
  sType <- getSurvivalType(pData(E))   
  surv.time.all <- as.numeric(as.character(pData(E)[ ,sType[1]]))
  surv.event.all <- as.numeric(as.character(pData(E)[ ,sType[2]]))
  mygene <- exprs(E)[rownames(exprs(E))==probe,]

  mygroup <- rep(NA,length(mygene))
  patientsSubType <- as.character(subtypeData$subtype)
  for(i in 1:length(subtypesOrderAll)){
    if(i != length(subtypesOrderAll)){
      gg <- patientsSubType == subtypesOrder[i]
      mygroup[gg] <- i
    }
  }
  tc <- 10 * 365
  surv.data <- censor.time(surv.time=surv.time.all / 365, surv.event=surv.event.all, time.cens=tc / 365)
  dd <- data.frame("time"=surv.data[[1]], "event"=surv.data[[2]], "gg"=mygroup)
  km.coxph.plot(formula.s=formula(Surv(time, event) ~ gg), data.s=dd, sub.s="all", x.label="Time (years)", y.label="Probability of Survival", main.title=sprintf("Probability of Survival for %s",probe), sub.title="", leg.text=subtypesOrderAll[1:3], leg.pos="bottomright", leg.inset=0.05,  o.text=FALSE, v.line=FALSE, h.line=FALSE, .col=c("darkblue", "darkgreen", "darkred"), .lty=1, show.n.risk=FALSE, n.risk.step=1, n.risk.cex=0.85, verbose=FALSE, ylim=c(0.3,1))
  ##dev.off()
}


###################################################
### combining estimations (cindex, dindex, hratio)
###################################################
## TODO
gxa.comb.est <- function(estName = NULL, estData = NULL, probe = NULL, subtypeData = NULL, subtypesOrderAll = NULL) {
  if(is.null(estName) ||is.null(estData) || is.null(probe) || is.null(subtypeData) || is.null(subtypesOrderAll)){
    stop("gxa.tdrocc.plot: You must specify estName, estData, probe, subtype data and subtype order!")
  }
  
  myest <- estData[rownames(estData[,,,estName,1]) == probe,,,estName,]
  tt <- as.data.frame(t(rbind(apply(myest,2,function(x)x[1]))))
  colnames(tt) <- "est"
  tt$se <- t(rbind(apply(myest,2,function(x)x[2])))

  cindices.meta <- apply(X=cindices, MARGIN=c(1, 4), FUN=combcind)
  
  combcind <- function(x) {
  	rr <- unlist(combine.est(x=x["cindex", ], x.se=x["cindex.se", ], na.rm=TRUE))
  	## compute two-sided p-values
  	statt <- (rr[1] - 0.5)/rr[2]
  	rr <- c(rr, pnorm(statt, lower.tail = rr[1] < 0.5) * 2)
  	names(rr) <- c("cindex", "cindex.se", "cindex.p")
  	return(rr)
  }

}



load("gxaData.RData")
load("gxaESets.RData")

E <- gxaESets[["E-MAIN-1234"]]
probe <- "1053_at"
subtypesOrderAll <- subtypesOrderAll
subtypeData <- gxaData[["E-MAIN-1234"]][[2]]
timeData <- 5
estName <- "cindices"
estData <- gxaData[["E-MAIN-1234"]][[1]]
experimentNames <- c("E-MAIN-1234", "E-TRAN-1234", "E-UPPP-1234", "E-UNTT-1234", "E-VDXX-1234")

gxa.boxplot(E, probe, subtypeData, subtypesOrderAll)
gxa.sbrier.plot(E, probe, subtypeData, subtypesOrderAll)
gxa.tdrocc.plot(E, probe, subtypeData, subtypesOrderAll, timeData)
gxa.forest.plot("cindices", estData, probe, subtypeData, subtypesOrderAll)
gxa.km.plot(E, probe, subtypeData, subtypesOrderAll)

## gxa.forest.plot("dindices", gxaData[["E-MAIN-1234"]][[1]], "1053_at", gxaData[["E-MAIN-1234"]][[2]], subtypesOrderAll)
## gxa.forest.plot("hratios", gxaData[["E-MAIN-1234"]][[1]], "1053_at", gxaData[["E-MAIN-1234"]][[2]], subtypesOrderAll)