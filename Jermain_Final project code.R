setwd("C:/Users/w10007346/Dropbox/Shark Dataset")
shark.length<-read.csv("Shark_Length.csv")
shark.cpue<-read.csv("Shark_cpue.csv")
require(dplyr)
require(tidyr)


################# Data manipulation 

shark.cpue.v2<-rename(shark.cpue, one=CARCHARHINUS.ACRONOTUS, 
                      two=CARCHARHINUS.BREVIPINNA, three=CARCHARHINUS.FALCIFORMIS,
                      four=CARCHARHINUS.LEUCAS, five=CARCHARHINUS.LIMBATUS, six=CARCHARHINUS.PLUMBEUS, 
                      seven=GALEOCERDO.CUVIER, eight=RHIZOPRIONODON.TERRAENOVAE)

shark.cpue.gath<-gather(shark.cpue.v2, "TAXON", "CPUE", one:eight)



six<-dplyr::filter(shark.cpue.gath, TAXON=="six")


plot(six$TURBBOTM,six$CPUE)

##################### Occupancy model data

six<-dplyr::filter(shark.cpue.gath, TAXON=="six")
six$CPUE[which(six$CPUE!=0)]<-1

six.3<-select(six, OXYBOTM, SALBOTM, TEMPBOTM, ENDDEPTH, TURBBOTM)
six.3$OXYBOTM2<-six.3$OXYBOTM*six.3$OXYBOTM
six.3$SALBOTM2<-six.3$SALBOTM*six.3$SALBOTM
six.3$TEMPBOTM2<-six.3$TEMPBOTM*six.3$TEMPBOTM
six.3$ENDDEPTH2<-six.3$ENDDEPTH*six.3$ENDDEPTH
six.3$TURBBOTM2<-six.3$TURBBOTM*six.3$TURBBOTM
head(six.3)

six.3<-scale(six.3)
six.3<-as.data.frame(six.3)
six.3$dombot4<-six$dombot4

six.3$CPUE<-six$CPUE
six.3<-na.omit(six.3)
View(six.3)
six.3[!complete.cases(six.3),]
length(six.3$OXYBOTM)


############# Dummy code substrate 

require(psych)

sub<-dummy.code(six.3$dombot4)
new.six<-data.frame(sub, six.3)
head(new.six)



write.csv(new.six, file="Pres_nodetect_2.csv", row.names=F)

############################################ Abundance data

six<-dplyr::filter(shark.cpue.gath, TAXON=="six")
six.4<-subset(six,six$CPUE>0)
six.3<-select(six.4, OXYBOTM, SALBOTM, TEMPBOTM, ENDDEPTH, TURBBOTM)
six.3$OXYBOTM2<-six.3$OXYBOTM*six.3$OXYBOTM
six.3$SALBOTM2<-six.3$SALBOTM*six.3$SALBOTM
six.3$TEMPBOTM2<-six.3$TEMPBOTM*six.3$TEMPBOTM
six.3$ENDDEPTH2<-six.3$ENDDEPTH*six.3$ENDDEPTH
six.3$TURBBOTM2<-six.3$TURBBOTM*six.3$TURBBOTM
head(six.3)

six.3<-scale(six.3)
six.3<-as.data.frame(six.3)
six.3$dombot4<-six.4$dombot4

six.3$CPUE<-six.4$CPUE
six.3<-na.omit(six.3)
View(six.3)
six.3[!complete.cases(six.3),]
length(six.3$OXYBOTM)


############# Dummy code substrate 

require(psych)

sub<-dummy.code(six.3$dombot4)
new.six<-data.frame(sub, six.3)
head(new.six)



write.csv(new.six, file="six.8.csv", row.names=F)


################################# Jags model for degree of presence model

library(data.table)
library(R2OpenBUGS)

library(rjags)
library(coda)

library(bayesplot)

setwd("C:/Users/Nathaniel/Dropbox/Bayes")
source("./post_summ_function.R")

data1 <- read.csv("six.8.csv", stringsAsFactors = F)
head(data1)
length(data1$OXYBOTM)

data = list(
  SALBOTM = data1$SALBOTM,
  SALBOTM2= data1$SALBOTM2,
  TURBBOTM2= data1$TURBBOTM2,
  TURBBOTM = data1$TURBBOTM,
  Gravel= data1$Gravel,
  Sand= data1$Sand,
  Rock = data1$Rock,
  CPUE = data1$CPUE
)


#Pooled model
mod = function(){
  #priors
  b0~dnorm(0,.001)
  b3~dnorm(0,.001)
  b4~dnorm(0,.001)
  b7~dnorm(0,.001)
  b8~dnorm(0,.001)
  b9~dnorm(0,.001)
  b10~dnorm(0,.001)
  b11~dnorm(0,.001)
  
  sigma~dunif(0,100)
  tau<-1/(sigma*sigma)
  
  
  
  #likelihood
  for(i in 1:321){
    mu[i]<-b0+b3*SALBOTM[i]+b4*SALBOTM2[i]+b7*Gravel[i]+b8*Sand[i]+b9*Rock[i]+
      b10*TURBBOTM[i]+b11*TURBBOTM2[i]
    
    CPUE[i]~dlnorm(mu[i], tau)
    CPUE.sim[i]~dlnorm(mu[i], tau)
  }
  mean.data<-mean(CPUE[])
  mean.sim<-mean(CPUE.sim[])
  p.mean<-step(mean.sim-mean.data)
  
  
  sd.data<-sd(CPUE[])
  sd.sim<-sd(CPUE.sim[])
  p.sd<-step(sd.sim-sd.data)
}

model.file="model.txt"
write.model(mod,model.file)


inits<-NULL

params = c("tau","b0", "b3", "b4", "b7", "b8", "b9", "b10", "b11","p.mean", "CPUE.sim", "p.sd")
ni = 100000
nb = 1000
nt = 1
nc = 3

jmod = jags.model(file = model.file, data = data, n.chains = nc, inits = inits, n.adapt = 1000)
update(jmod, n.iter=nb, by=1)
post = coda.samples(jmod, params, n.iter = ni, thin = nt)
mcmc <- as.data.frame(rbindlist(lapply(post, as.data.frame)))
plot(post)
summary(post)
gelman.diag(post)
post.summ(post, "b3")
post.summ(post, "b4")
post.summ(post, "b7")
post.summ(post, "b8")
post.summ(post, "b9")
post.summ(post, "b10")
post.summ(post, "b11")
post.summ(post, "p.mean")
post.summ(post, "p.sd")

# make quantiles   #############################
post.mat<-as.matrix(post)


ci.lower<-0
ci.upper<-0
pred.median<-0

for (i in 1:321){
  
  pred.median[i]<-median(post.mat[,i])
  ci.lower[i]<-quantile(post.mat[,i], c(0.025))
  ci.upper[i]<-quantile(post.mat[,i], c(0.975))
}

ci.upper

oneto<-seq(0,100)
one<-seq(0,100)

par(oma=c(.5,.5,0,0))
plot(x=data1$CPUE, y=pred.median, xlim=c(0,16), ylim=c(0,8), col="red", cex=1.5, xlab="Observed CPUE", 
     ylab="Predicted CPUE", cex.lab=1.5, cex.axis=1.5, las=1, pch=19, bg= "red")
lines(x=oneto, y=one, lwd=2)

for (i in 1:321){
  arrows(data1$CPUE[i], ci.lower[i], data1$CPUE[i], ci.upper[i], angle=FALSE ,code=3, col="blue" )
}
legend(11,8, legend ="Predicted Median", pch=19, col = "red", cex=1)
legend(11,6.5, legend = "Credible Interval", lty = 1, lwd = 2, col = "blue", cex=1)

### Histogram with fit #################################

post.summ(post, "p.mean")
CPUE.pred <- post.summ(post, "CPUE.sim[")

par(oma=c(.5,.5,.5,.5))
hist(data1$CPUE, freq = FALSE, ylim = c(0,1), main ="", xlab = "CPUE", breaks=20, 
     cex.lab=1.5, cex.axis=1.5, las=1, xlim=c(0,10))
lines(density(CPUE.pred["mean",]), type = "l", col = "red", lwd = 2)
legend("right", legend = "Predicted", lty = 1, lwd = 2, col = "red", cex=1.5)

#### DIC for this model


dic.mod1 <- dic.samples(jmod, 10000, "pD")




########################################################## Occupancy model in Jags


setwd("C:/Users/w10007346/Dropbox/Bayes")
data1 <- read.csv("Pres_nodetect_2.csv", stringsAsFactors = F)
head(data1)
length(data1$TURBBOTM)

data = list(
  SAL = data1$SALBOTM,
  SAL2= data1$SALBOTM2,
  TURB = data1$TURBBOTM,
  TURB2= data1$TURBBOTM2,
  TEMP = data1$TEMPBOTM,
  TEMP2= data1$TEMPBOTM2,
  Gravel= data1$Gravel,
  Sand= data1$Sand,
  Rock = data1$Rock,
  pres = data1$CPUE,
  OXY = data1$OXYBOTM,
  OXY2 = data1$OXYBOTM2
)



mod = function(){
  
  #priors
  b0~dnorm(0,.001)
  b1~dnorm(-.36,1)
  b2~dnorm(.37,1)
  b3~dnorm(-.01,1)
  b4~dnorm(.09,.1)
  b5~dnorm(-.16,1)
  b6~dnorm(-.10,1)
  b7~dnorm(-.15,1)
  b8~dnorm(0,.001)
  b9~dnorm(0,.001)
  b10~dnorm(0,.001)
  b11~dnorm(0,.001)
  
  
  
  
  
  #likelihood
  for ( i in 1:2230){
    logit(P[i])<-b0+b1*SAL[i]+ b2*SAL2[i]+b3*TURB[i]+b4*TURB2[i]+
      b5*Gravel[i]+b6*Sand[i]+b7*Rock[i]+b8*TEMP[i]+b9*TEMP2[i]+b10*OXY[i]+b11*OXY2[i]
    
    pres[i]~dbern(P[i])
    pres.sim[i]~dbern(P[i])
    
    
  }
  mean.data<-mean(pres[])
  mean.sim<-mean(pres.sim[])
  p.mean<-step(mean.sim-mean.data)
  
  sd.data<-sd(pres[])
  sd.sim<-sd(pres.sim[])
  p.sd<-step(sd.sim-sd.data)
  
}



model.file="model.txt"
write.model(mod,model.file)


inits<-NULL

params = c("b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "b11", "p.mean", "pres.sim", "p.sd")
ni = 10000
nb = 1000
nt = 1
nc = 3

jmod = jags.model(file = model.file, data = data, n.chains = nc, inits = inits, n.adapt = 1000)
update(jmod, n.iter=nb, by=1)
post = coda.samples(jmod, params, n.iter = ni, thin = nt)
mcmc <- as.data.frame(rbindlist(lapply(post, as.data.frame)))
plot(post)
summary(post)
gelman.diag(post)
post.summ(post, "b1")
post.summ(post, "b2")
post.summ(post, "b3")
post.summ(post, "b4")
post.summ(post, "b5")
post.summ(post, "b6")
post.summ(post, "b7")
post.summ(post, "b8")
post.summ(post, "b9")
post.summ(post, "b10")
post.summ(post, "b11")
post.summ(post, "p.mean")
post.summ(post, "p.sd")


### DIC for this model

dic.mod1 <- dic.samples(jmod, 5000, "pD")

### plot fit 

post.summ(post, "p.mean")
pres.pred <- post.summ(post, "pres.sim[")

par(oma=c(.5,.5,.5,.5))
hist(data1$CPUE, freq=FALSE, main ="", xlab = "Prob(Presence)", 
     cex.lab=1.5, cex.axis=1.5, las=1, ylim=c(0,10), xlim=c(0,1))
lines(density(pres.pred["mean",]), type = "l", col = "red", lwd = 2, xlim=c(0,1))
legend("right", legend = "Predicted", lty = 1, lwd = 2, col = "red", cex=1.5)

# make quantiles   #############################
post.mat<-as.matrix(post)


ci.lower<-0
ci.upper<-0
pred.median<-0

for (i in 1:2230){
  
  pred.median[i]<-median(post.mat[,i])
  ci.lower[i]<-quantile(post.mat[,i], c(0.025))
  ci.upper[i]<-quantile(post.mat[,i], c(0.975))
}

ci.upper

oneto<-seq(0,100)
one<-seq(0,100)

par(oma=c(.5,.5,.5,.5))
plot(x=data1$CPUE, y=pred.median, col="red", cex=1.5, xlab="Observed CPUE", 
     ylab="Predicted CPUE", cex.lab=1.5, cex.axis=1.5)
lines(x=oneto, y=one, lwd=2)

for (i in 1:321){
  arrows(data1$CPUE[i], ci.lower[i], data1$CPUE[i], ci.upper[i], angle=FALSE ,code=3, col="blue" )
}





