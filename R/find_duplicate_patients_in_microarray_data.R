library(survcomp)
library(Biobase)
data(breastCancerData)

cinfo <- colnames(pData(mainz7g))
##c("samplename", "dataset", "series", "id","filename","size", "age", "er", "grade", "pgr", "her2", "brca.mutation", "e.dmfs", "t.dmfs", "node", "t.rfs", "e.rfs", "treatment", "tissue", "t.os", "e.os")
data.all <- c("mainz7g"=mainz7g, "transbig7g"=transbig7g, "upp7g"=upp7g, "unt7g"=unt7g, "vdx7g"=vdx7g, "nki7g"=nki7g)


idtoremove.all <- NULL
duplres <- NULL
demo.all <- NULL

for(i in data.all) {
	demo.all <- rbind(demo.all, pData(i))
}

dn2 <- unique(demo.all$dataset)
duplid <- sort(unique(demo.all[duplicated(demo.all[ , "id"]), "id"]))
duplrest <- NULL
for(i in 1:length(duplid)) {
	tt <- NULL
	for(k in 1:length(dn2)) {
		myx <- sort(row.names(demo.all)[complete.cases(demo.all[ , c("id", "dataset")]) & demo.all[ , "id"] == duplid[i] & demo.all[ , "dataset"] == dn2[k]])
		if(length(myx) > 0) { tt <- c(tt, myx) }
	}
	duplrest <- c(duplrest, list(tt))
}
names(duplrest) <- duplid
duplres <- c(duplres, duplrest)

res <- matrix("", ncol=max(sapply(duplres, length)), nrow=length(duplres), dimnames=list(names(duplres), paste("duplicate", 1:max(sapply(duplres, length)), sep=".")))
for(i in 1:length(duplres)) { res[i, 1:length(duplres[[i]])] <- duplres[[i]] }

duplall <- c(res[,2:ncol(res)])
removeEmpty <- duplall != ""
duplall <- duplall[removeEmpty,]

newPL <- as.list(NULL)
patientsList <- patID <- NULL
count <- 1

for(i in data.all) {
	patientsList <- NULL
	for(j in 1:length(rownames(pData(i)))){
		tfMatrix <- duplall == rownames(pData(i)[j,])
		patientsList[j] <- is.na(table(tfMatrix)[2])[[1]]
	}
	newPL[[count]] <- patientsList
	count <- count + 1
}

names(newPL) <- names(data.all)