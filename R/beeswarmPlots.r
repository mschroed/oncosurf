##for mainz7g:
library(beeswarm)
library(survcomp)
library(Biobase)
data(breastCancerData)

par(mfrow=c(1,3))

breast2 <- NULL
breast2$ER <- as.factor(pData(mainz7g)[,"er"])
breast2$ESR1 <- as.numeric(exprs(mainz7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(mainz7g)[2,])
breast2$time_survival <- as.integer(pData(mainz7g)[,"t.dmfs"])
breast2$event_survival <- as.integer(pData(mainz7g)[,"e.dmfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(ER),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$ER), title = 'ER', pch = 16, col = 1:2)


breast2 <- NULL
breast2$grade <- as.factor(pData(mainz7g)[,"grade"])
breast2$ESR1 <- as.numeric(exprs(mainz7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(mainz7g)[2,])
breast2$time_survival <- as.integer(pData(mainz7g)[,"t.dmfs"])
breast2$event_survival <- as.integer(pData(mainz7g)[,"e.dmfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(grade),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$grade), title='Grade', pch=16, col=1:3)


breast2 <- NULL
breast2$age <- as.factor(pData(mainz7g)[,"age"])
breast2$ESR1 <- as.numeric(exprs(mainz7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(mainz7g)[2,])
breast2$time_survival <- as.integer(pData(mainz7g)[,"t.dmfs"])
breast2$event_survival <- as.integer(pData(mainz7g)[,"e.dmfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(age),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$age), title='Age', pch=16, col=1:length(unique(breast2$age)))


##for transbig7g
par(mfrow=c(1,3))

breast2 <- NULL
breast2$ER <- as.factor(pData(transbig7g)[,"er"])
breast2$ESR1 <- as.numeric(exprs(transbig7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(transbig7g)[2,])
breast2$time_survival <- as.integer(pData(transbig7g)[,"t.dmfs"])
breast2$event_survival <- as.integer(pData(transbig7g)[,"e.dmfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(ER),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$ER), title = 'ER', pch = 16, col = 1:2)


breast2 <- NULL
breast2$grade <- as.factor(pData(transbig7g)[,"grade"])
breast2$ESR1 <- as.numeric(exprs(transbig7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(transbig7g)[2,])
breast2$time_survival <- as.integer(pData(transbig7g)[,"t.dmfs"])
breast2$event_survival <- as.integer(pData(transbig7g)[,"e.dmfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(grade),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$grade), title='Grade', pch=16, col=1:3)


breast2 <- NULL
breast2$age <- as.factor(pData(transbig7g)[,"age"])
breast2$ESR1 <- as.numeric(exprs(transbig7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(transbig7g)[2,])
breast2$time_survival <- as.integer(pData(transbig7g)[,"t.dmfs"])
breast2$event_survival <- as.integer(pData(transbig7g)[,"e.dmfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(age),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$age), title='Age', pch=16, col=1:length(unique(breast2$age)))


## for upp7g

par(mfrow=c(1,3))

breast2 <- NULL
breast2$ER <- as.factor(pData(upp7g)[,"er"])
breast2$ESR1 <- as.numeric(exprs(upp7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(upp7g)[2,])
breast2$time_survival <- as.integer(pData(upp7g)[,"t.rfs"])
breast2$event_survival <- as.integer(pData(upp7g)[,"e.rfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(ER),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$ER), title = 'ER', pch = 16, col = 1:2)


breast2 <- NULL
breast2$grade <- as.factor(pData(upp7g)[,"grade"])
breast2$ESR1 <- as.numeric(exprs(upp7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(upp7g)[2,])
breast2$time_survival <- as.integer(pData(upp7g)[,"t.rfs"])
breast2$event_survival <- as.integer(pData(upp7g)[,"e.rfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(grade),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$grade), title='Grade', pch=16, col=1:3)


breast2 <- NULL
breast2$age <- as.factor(pData(upp7g)[,"age"])
breast2$ESR1 <- as.numeric(exprs(upp7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(upp7g)[2,])
breast2$time_survival <- as.integer(pData(upp7g)[,"t.rfs"])
breast2$event_survival <- as.integer(pData(upp7g)[,"e.rfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(age),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$age), title='Age', pch=16, col=1:length(unique(breast2$age)))



##unt7g
par(mfrow=c(1,3))

breast2 <- NULL
breast2$ER <- as.factor(pData(unt7g)[,"er"])
breast2$ESR1 <- as.numeric(exprs(unt7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(unt7g)[2,])
breast2$time_survival <- as.integer(pData(unt7g)[,"t.dmfs"])
breast2$event_survival <- as.integer(pData(unt7g)[,"e.dmfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(ER),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$ER), title = 'ER', pch = 16, col = 1:2)


breast2 <- NULL
breast2$grade <- as.factor(pData(unt7g)[,"grade"])
breast2$ESR1 <- as.numeric(exprs(unt7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(unt7g)[2,])
breast2$time_survival <- as.integer(pData(unt7g)[,"t.dmfs"])
breast2$event_survival <- as.integer(pData(unt7g)[,"e.dmfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(grade),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$grade), title='Grade', pch=16, col=1:3)


breast2 <- NULL
breast2$age <- as.factor(pData(unt7g)[,"age"])
breast2$ESR1 <- as.numeric(exprs(unt7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(unt7g)[2,])
breast2$time_survival <- as.integer(pData(unt7g)[,"t.dmfs"])
breast2$event_survival <- as.integer(pData(unt7g)[,"e.dmfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(age),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$age), title='Age', pch=16, col=1:length(unique(breast2$age)))


## for vdx7g

par(mfrow=c(1,3))

breast2 <- NULL
breast2$ER <- as.factor(pData(vdx7g)[,"er"])
breast2$ESR1 <- as.numeric(exprs(vdx7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(vdx7g)[2,])
breast2$time_survival <- as.integer(pData(vdx7g)[,"t.dmfs"])
breast2$event_survival <- as.integer(pData(vdx7g)[,"e.dmfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(ER),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$ER), title = 'ER', pch = 16, col = 1:2)


breast2 <- NULL
breast2$grade <- as.factor(pData(vdx7g)[,"grade"])
breast2$ESR1 <- as.numeric(exprs(vdx7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(vdx7g)[2,])
breast2$time_survival <- as.integer(pData(vdx7g)[,"t.dmfs"])
breast2$event_survival <- as.integer(pData(vdx7g)[,"e.dmfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(grade),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$grade), title='Grade', pch=16, col=1:3)


breast2 <- NULL
breast2$age <- as.factor(pData(vdx7g)[,"age"])
breast2$ESR1 <- as.numeric(exprs(vdx7g)[1,])
breast2$ERBB2 <- as.numeric(exprs(vdx7g)[2,])
breast2$time_survival <- as.integer(pData(vdx7g)[,"t.dmfs"])
breast2$event_survival <- as.integer(pData(vdx7g)[,"e.dmfs"])

beeswarm(time_survival ~ event_survival, data = breast2,method = 'smile',pch = 16, pwcol = as.numeric(age),xlab = '', ylab = 'Follow-up time (days)',labels = c('Censored', 'Metastasis'))
legend('topright', legend = levels(breast2$age), title='Age', pch=16, col=1:length(unique(breast2$age)))