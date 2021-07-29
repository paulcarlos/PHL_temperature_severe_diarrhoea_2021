library(dlnm)
library(splines)
library(mixmeta)

# import data
dlist <- readRDS("ovr_data_wkly.rds")

# model specification
varfun <- "ns"
vardf <- 3
lagnk <- 3
tlag <- 4
dfseas <- 3
yr <- 12

# model formula
fmla <- as.formula(mort~cb+ns(wk,df=dfseas*yr)+holiday2+ns(dewp,df=3)) # for mortality
fmla <- as.formula(dhn2~cb+ns(wk,df=dfseas*yr)+holiday2+ns(dewp,df=3)) # severe dehydration


# blank coefficients and vcov matrices
dist <- names(dlist)
b.coef <- matrix(data=NA,nrow=length(dist),ncol=vardf,dimnames=list(dist))
b.vcov <- vector("list",length(dist)); names(b.vcov) <- dist



################################
############# BLUP ############# 
coef1 <- b.coef
vcov1 <- b.vcov

for (j in 1:length(dist)) {
  dat <- dlist[[dist[j]]]
  var1 <- dat$tmean
  
  # model and prediction
  argvar <- list(fun=varfun,knots=quantile(var1,c(0.33,0.67)))
  arglag <- list(fun=varfun,knots=logknots(tlag,nk=lagnk))
  cb <- crossbasis(var1,lag=tlag,argvar=argvar,arglag=arglag)
  mod <- glm(fmla,dat,family=quasipoisson)
  
  # reduction
  red <- crossreduce(cb,mod,cen=median(var1))
  coef1[j,] <- coef(red)
  vcov1[[j]] <- vcov(red)
  rm(dat,argvar,arglag,cb,mod,red,var1)
}
m1 <- mixmeta(formula=coef1~1,S=vcov1,control=list(showiter=TRUE))
summary(m1)
blup1 <- blup(m1,vcov=TRUE)



###########################################
############# NCR level plot ##############

# model and prediction
t1 <- unlist(lapply(dlist[dist],function(x)x$tmean))
tpct <- round(quantile(t1,c(0.01,0.05,0.5,0.95,0.99)),1)
argvar1 <- list(x=t1,fun=varfun,knots=quantile(t1,c(0.33,0.67)))
bvar1 <- do.call(onebasis,argvar1)
pred1 <- crosspred(bvar1,coef=m1$coefficients,vcov=m1$vcov,
                   model.link="log",by=0.1,cen=tpct["50%"])
cen1 <- as.numeric(names(which.min(pred1$allRRfit[as.character(seq(tpct["1%"],tpct["99%"],by=0.1))])))
pred1 <- crosspred(bvar1,coef=m1$coefficients,vcov=m1$vcov,
                   model.link="log",by=0.1,cen=cen1)

# plot
plot(pred1,ylab="relative risk",ylim=c(0,3),xlab="weekly mean temperature (°C)",xaxt="n",
     main="Mortality due to diarrhoea",
     lwd=4,col="red",ci.arg=list(density=30,angle=45,col="pink"))
axis(1,at=seq(floor(min(t1)),ceiling(max(t1)),by=0.5))
breaks <- c(min(t1)-1,seq(pred1$predvar[1],
                          pred1$predvar[length(pred1$predvar)],length=30),max(t1)+1)
hist <- hist(t1,breaks=breaks,plot=F)
hist$density <- hist$density/max(hist$density)*0.7
prop <- max(hist$density)/max(hist$counts)
counts <- pretty(hist$count,3)
plot(hist,ylim=c(0,max(hist$density)*3.5),axes=F,ann=F,col="lightgray",
     breaks=breaks,freq=F,add=T)
axis(4,at=counts*prop,labels=counts,cex.axis=0.7)
abline(v=cen1,col="red",lty=1)
abline(v=quantile(t1,c(0.05,0.95)),col=grey(0.3),lty=3)



#################################################
############# district level plots ##############

par(mfrow=c(2,2),mar=c(3,3,4,3),oma=c(3,3,1,1))
loc <- c("First District","Second District","Fourth District")
for (i in 1:length(dist)) {
  # data selection
  dat <- dlist[[dist[i]]]
  var1 <- dat$tmean
  vpct <- round(quantile(var1,c(0.01,0.05,0.5,0.95,0.99)),1)
  
  # centring at minimum risk temperature
  predvar <- quantile(var1,1:99/100,na.rm=T)
  argvar1 <- list(x=predvar,fun=varfun,knots=quantile(var1,c(0.33,0.67)),Bound=range(var1,na.rm=T))
  bvar <- do.call(onebasis,argvar1)
  minpercreg1 <- (1:99)[which.min((bvar%*%blup1[[i]]$blup))]
  cen <- round(quantile(var1,minpercreg1/100,na.rm=T),1)
  
  # model
  argvar <- list(x=var1,fun=varfun,knots=quantile(var1,c(0.33,0.67)))
  bvar <- do.call(onebasis,argvar)
  cp <- crosspred(bvar,coef=blup1[[i]]$blup,vcov=blup1[[i]]$vcov,
                  model.link="log",by=0.1,cen=cen)
  
  #plot
  plot(cp,"overall",ylim=c(0,3),lwd=3,col="red",ci.arg=list(density=30,angle=45,col="pink"),
       xaxt="n",xlab="", ylab="",main=loc[i])
  axis(1,at=seq(floor(min(var1)),ceiling(max(var1)),by=1))
  breaks <- c(min(var1)-1,seq(pred1$predvar[1],
                              pred1$predvar[length(pred1$predvar)],length=30),max(var1)+1)
  hist <- hist(var1,breaks=breaks,plot=F)
  hist$density <- hist$density/max(hist$density)*0.7
  prop <- max(hist$density)/max(hist$counts)
  counts <- pretty(hist$count,3)
  plot(hist,ylim=c(0,max(hist$density)*3.5),axes=F,ann=F,col="lightgray",
       breaks=breaks,freq=F,add=T)
  axis(4,at=counts*prop,labels=counts,cex.axis=0.7)
  abline(v=cen,col="red",lty=1,)
  rm(predvar,bvar,argvar1,minpercreg1,argvar,cen,cp,var1)
}
mtext(text="weekly mean temperature (°C)",side=1,line=1,outer=TRUE)
mtext(text="relative risk",side=2,line=1,outer=TRUE)


