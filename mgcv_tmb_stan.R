library(rdhs)
library(demogsurv)
library(ggplot2)
library(data.table)
library(MASS)
library(survival)

###DATA
countries <- dhs_countries()
load.dat<-function(x){
  cc <- countries[CountryName == x]$DHS_CountryCode
  surveys <- dhs_surveys(countryIds = cc, surveyYearStart=2005, surveyType = "DHS")
  ird <- dhs_datasets(fileType = "IR", fileFormat = "flat")[SurveyId %in% surveys$SurveyId]
  ird$path <- unlist(get_datasets(ird$FileName))
  ir <- list()
  for(survid in ird$SurveyId){
    print(survid)
    dat <- readRDS(ird[SurveyId == survid]$path)
    dat <- dat[grep("caseid|^v0|^v1|^b|^mm", names(dat))]
    ir[[survid]] <- dat
  }
  ir <- lapply(ir, haven::as_factor)
  ir <- Map(data.frame,
            SurveyId = surveys$SurveyId,
            CountryName = surveys$CountryName,
            SurveyYear = surveys$SurveyYear,
            ir)
  
  zw<-lapply(ir,reshape_sib_data)
  for(i in 1:length(zw)){
    zw[[i]]$death<-factor(zw[[i]]$mm2,c("dead","alive"))=="dead"
    zw[[i]]$tstop<-ifelse(zw[[i]]$death,zw[[i]]$mm8,zw[[i]]$v008)
    zw[[i]]$weight<-zw[[i]]$v005/1e6
  }
  
  aggr<-list()
  for(i in 1:length(zw))
    aggr[[names(zw)[i]]] <- demog_pyears(~mm1, zw[[i]],
                                         period = -15:0+as.numeric(surveys$SurveyYear[i]),
                                         agegr = seq(0,70),
                                         tips = 0:15,
                                         dob = "mm4",
                                         intv = "v008",
                                         event = "death",
                                         tstart = "mm4",
                                         tstop = "tstop",
                                         weights = "weight")$data
  
  return(aggr)
}

aggr<-load.dat("Zimbabwe")

aggr.mat<-c()
for(i in 1:length(aggr)){
  aggr.mat<-rbind(aggr.mat,aggr[[i]])
}

aggr.mat<-aggr.mat[!aggr.mat$mm1=="missing",]
aggr.mat.reduced<-aggr.mat[aggr.mat$pyears!=0,]
aggr.mat.reduced$agegr<-as.numeric(levels(aggr.mat.reduced$agegr))[aggr.mat.reduced$agegr]
aggr.mat.reduced$period<-as.numeric(levels(aggr.mat.reduced$period))[aggr.mat.reduced$period]
aggr.mat.reduced<-aggr.mat.reduced[order(aggr.mat.reduced$mm1,aggr.mat.reduced$period,aggr.mat.reduced$tips,aggr.mat.reduced$agegr),]
aggr.mat.reduced$tips<-as.factor(as.numeric(levels(aggr.mat.reduced$tips))[aggr.mat.reduced$tips])
aggr.mat.reduced$tips<-relevel(aggr.mat.reduced$tips,ref=3)


age.start<-10;age.end<-50
year.start<-1995;year.end<-2014

###MGCV

###2D spline with mgcv, UBRE/AIC and GCV tend to underfit
#fit<-mgcv::gam(event~te(agegr,period,bs="ps",k=15,by=as.factor(mm1))+as.factor(mm1)+as.factor(tips),offset=log(pyears),data=aggr.mat.reduced[aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,],family=poisson(link=log))
#fit.gcv.gamma1<-mgcv::gam(event~te(agegr,period,bs="ps",k=15,by=as.factor(mm1))+as.factor(mm1)+as.factor(tips),offset=log(pyears),data=aggr.mat.reduced[aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,],family=poisson(link=log),scale=-1)

###inflate EDF, set gamma=1.4
#fit.gcv<-mgcv::gam(event~te(agegr,period,bs="ps",k=15,by=as.factor(mm1))+as.factor(mm1)+as.factor(tips),offset=log(pyears),data=aggr.mat.reduced[aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,],family=poisson(link=log),scale=-1,gamma=1.4)

###ReML, round deaths to integers
system.time(fit.reml<-mgcv::gam(round(event)~te(agegr,period,bs="ps",k=15,by=as.factor(mm1))+as.factor(mm1)+as.factor(tips),offset=log(pyears),data=aggr.mat.reduced[aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,],family=poisson(link=log),method="REML"))

###TMB
no.basis=15
knots<-seq(0,1,length=no.basis-2)
dk<-knots[2]-knots[1]	
knots<-c(knots[1]-dk*(3:1),knots,knots[no.basis-2]+dk*(1:3))
age<-seq(0,1,length=age.end-age.start+1)
year<-seq(0,1,length=year.end-year.start+1)
cohort<-seq(0,1,length=year.end-age.start-year.start+age.end+1)

bspline<-function (x,k,i,m=2) {
  if (m==-1) {basis<-as.numeric(x<k[i+1] & x>=k[i])} else {
    z0<-(x-k[i])/(k[i+m+1]-k[i])
    z1<-(k[i+m+2]-x)/(k[i+m+2]-k[i+1])
    basis<-z0*bspline(x,k,i,m-1)+z1*bspline(x,k,i+1,m-1) }
  basis
}

A.age<-c()
for(j in 1:no.basis) {
  A.age<-cbind(A.age,bspline(age,knots,j))
}

A.year<-c()
for(j in 1:no.basis) {
  A.year<-cbind(A.year,bspline(year,knots,j))
}

te.spline<-A.year%x%A.age

P<-diff(diag(no.basis),differences=2)
P1<-diff(diag(no.basis))
S<-crossprod(P)
penal.age<-diag(no.basis)%x%S
penal.time<-S%x%diag(no.basis)
penal.age.time<-crossprod(P1%x%P1)

#Null penalties for Stan
u.null<-eigen(penal.age+penal.time)$vectors[,(no.basis^2-3):(no.basis^2)]
Q.null<-1e-3*tcrossprod(u.null)
u.null.agecrosstime<-eigen(penal.age+penal.time+penal.age.time)$vectors[,(no.basis^2-2):(no.basis^2)]
Q.null.agecrosstime<-1e-3*tcrossprod(u.null.agecrosstime)

library(TMB)
compile(".../2d_single_sex.cpp")
dyn.load(dynlib(".../2d_single_sex"))

data.mat.m<-aggr.mat.reduced[aggr.mat.reduced$mm1=="male" & aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,c(7,5,3,2,4)]
data.mat.m[,5]<-as.numeric(levels(data.mat.m[,5]))[data.mat.m[,5]]
DX.spline.m<-te.spline[(age.end-age.start+1)*(data.mat.m[,4]-year.start)+data.mat.m[,3]-age.start+1,]

data.mat.f<-aggr.mat.reduced[aggr.mat.reduced$mm1=="female" & aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,c(7,5,3,2,4)]
data.mat.f[,5]<-as.numeric(levels(data.mat.f[,5]))[data.mat.f[,5]]
DX.spline.f<-te.spline[(age.end-age.start+1)*(data.mat.f[,4]-year.start)+data.mat.f[,3]-age.start+1,]

tmb <- MakeADFun(data=list(dm=round(data.mat.m[,1]),Em=data.mat.m[,2],tpm=data.mat.m[,5],
                         df=round(data.mat.f[,1]),Ef=data.mat.f[,2],tpf=data.mat.f[,5],
                         Dm=DX.spline.m,Df=DX.spline.f,
                         penal_age=penal.age,penal_time=penal.time,null_penal=1e-3*diag(no.basis*no.basis)),
               parameters=list(log_lambda_age_m=3,log_lambda_time_m=3,log_lambda_age_f=3,log_lambda_time_f=3,spline_params_m=rep(0,no.basis*no.basis),spline_params_f=rep(0,no.basis*no.basis),tips_params=rep(0,14)),random=c("spline_params_m","spline_params_f"),
               DLL="2d_single_sex")
##optim BFGS default
system.time(op<- do.call("optim",tmb))
rep<-sdreport(tmb,getJointPrecision = TRUE)

###Stan
data.mat.m<-aggr.mat.reduced[aggr.mat.reduced$mm1=="male" & aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,c(7,5,3,2,4)]
data.mat.m[,5]<-as.numeric(levels(data.mat.m[,5]))[data.mat.m[,5]]
DX.spline.m<-te.spline[(age.end-age.start+1)*(data.mat.m[,4]-year.start)+data.mat.m[,3]-age.start+1,]
data.mat.f<-aggr.mat.reduced[aggr.mat.reduced$mm1=="female" & aggr.mat.reduced$agegr %in% age.start:age.end & aggr.mat.reduced$period %in% year.start:year.end,c(7,5,3,2,4)]
data.mat.f[,5]<-as.numeric(levels(data.mat.f[,5]))[data.mat.f[,5]]
DX.spline.f<-te.spline[(age.end-age.start+1)*(data.mat.f[,4]-year.start)+data.mat.f[,3]-age.start+1,]
max.tp<-max(c(data.mat.m[,5],data.mat.f[,5]))

mu.m<-data.mat.m[,1]+0.1
offset.m<-log(data.mat.m[,2])
QR<-qr(as.numeric(mu.m^0.5)*DX.spline.m)
QR.pivot<-QR$pivot
R<-qr.R(QR)
full.mat.m<-crossprod(R[,order(QR.pivot)])
XWz.m<-c(t(DX.spline.m)%*%(-0.1+mu.m*log(mu.m)-mu.m*offset.m))

mu.f<-data.mat.f[,1]+0.1
offset.f<-log(data.mat.f[,2])
QR<-qr(as.numeric(mu.f^0.5)*DX.spline.f)
QR.pivot<-QR$pivot
R<-qr.R(QR)
full.mat.f<-crossprod(R[,order(QR.pivot)])
XWz.f<-c(t(DX.spline.f)%*%(-0.1+mu.f*log(mu.f)-mu.f*offset.f))

stan.list<-list(k=no.basis,tips=max.tp,
                n_m=nrow(data.mat.m),n_f=nrow(data.mat.f),
                dx_m=round(data.mat.m[,1]),dx_f=round(data.mat.f[,1]),
                Ex_m=data.mat.m[,2],Ex_f=data.mat.f[,2],
                tp_m=data.mat.m[,5],tp_f=data.mat.f[,5],
                DX_m=DX.spline.m,DX_f=DX.spline.f,
                penal_age=penal.age,
                penal_time=penal.time,null_penal=Q.null,
                full_mat_m=full.mat.m,full_mat_f=full.mat.f,
                XWz_m=XWz.m,XWz_f=XWz.f)

stan.fit<- stan(file ='C:\\Users\\steve\\Desktop\\Imperial\\DHS 2 sexes uniform priors.stan', data = stan.list, iter = 3000, chains = 1)

####PLOTS
##simulate draws from estimated parameter distributions
param.sim<-t(mvrnorm(1000,fit.reml$coef,vcov(fit.reml))) ###mgcv, default Bayesian confidence intervals, conditional on estimated smoothing parameters
library(spam)
tmb.sim<-rmvnorm.prec(1000,mu=summary(rep)[,1],Q=rep$jointPrecision) ###tmb


stan.age<-function(x,DX,i){
  matrix(DX%*%x,age.end-age.start+1,year.end-year.start+1)[i-age.start+1,]
}
stan.year<-function(x,DX,i){
  matrix(DX%*%x,age.end-age.start+1,year.end-year.start+1)[,i-year.start+1]
}

plot.fun<-function(i,sex){
  if(sex=="Male") {
    stan.age.10.m<-t(apply(stan.params.m,1,stan.age,DX=te.spline,i=i))
    age.10.m<-t(apply(param.sim,2,stan.age,DX=Xm,i=i))
    tmb.age.10.m<-t(apply(tmb.sim[,19:(19+224)],1,stan.age,DX=te.spline,i=i))
  } else {
    stan.age.10.m<-t(apply(stan.params.f,1,stan.age,DX=te.spline,i=i))
    age.10.m<-t(apply(param.sim,2,stan.age,DX=Xf,i=i))
    tmb.age.10.m<-t(apply(tmb.sim[,(19+225):(19+225+224)],1,stan.age,DX=te.spline,i=i))
  }
  df<-data.frame(stan.mean=apply(stan.age.10.m,2,mean),
                 stan.025=apply(stan.age.10.m,2,quantile,probs=0.025),
                 stan.975=apply(stan.age.10.m,2,quantile,probs=0.975),
                 mgcv.mean=apply(age.10.m,2,mean),
                 mgcv.025=apply(age.10.m,2,quantile,probs=0.025),
                 mgcv.975=apply(age.10.m,2,quantile,probs=0.975),
                 tmb.mean=apply(tmb.age.10.m,2,mean),
                 tmb.025=apply(tmb.age.10.m,2,quantile,probs=0.025),
                 tmb.975=apply(tmb.age.10.m,2,quantile,probs=0.975)
  )
  
  gg.colors<-c("mgcv"="black","TMB"="blue","Stan"="red")
  
  ggplot(df)+
    geom_line(aes(x=1995:2014,y=mgcv.mean,color="mgcv"),size=1)+geom_ribbon(aes(x=1995:2014,ymin=mgcv.025,ymax=mgcv.975,fill="mgcv"),alpha=0.4)+
    geom_line(aes(x=1995:2014,y=stan.mean,color="TMB"),size=1)+geom_ribbon(aes(x=1995:2014,ymin=stan.025,ymax=stan.975,fill="Stan"),alpha=0.2)+
    geom_line(aes(x=1995:2014,y=tmb.mean,color="Stan"),size=1)+geom_ribbon(aes(x=1995:2014,ymin=tmb.025,ymax=tmb.975,fill="TMB"),alpha=0.3)+
    scale_color_manual(values=gg.colors,guide=FALSE)+scale_fill_manual(values=gg.colors,guide=FALSE)+
    ylab("Log Mortality")+xlab("Year")+ggtitle(paste(sex,"Age",i))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),legend.key.size=unit(0.5,"cm"),legend.key=element_rect(fill=NA,colour=NA,size=0.5),legend.position=c(0.1,0.2),plot.title = element_text(hjust = 0.5))
}

plot.fun.year<-function(i,sex){
  if(sex=="Male") {
    stan.age.10.m<-t(apply(stan.params.m,1,stan.year,DX=te.spline,i=i))
    age.10.m<-t(apply(param.sim,2,stan.year,DX=Xm,i=i))
    tmb.age.10.m<-t(apply(tmb.sim[,19:(19+224)],1,stan.year,DX=te.spline,i=i))
  } else {
    stan.age.10.m<-t(apply(stan.params.f,1,stan.year,DX=te.spline,i=i))
    age.10.m<-t(apply(param.sim,2,stan.year,DX=Xf,i=i))
    tmb.age.10.m<-t(apply(tmb.sim[,(19+225):(19+225+224)],1,stan.year,DX=te.spline,i=i))
  }
  df<-data.frame(stan.mean=apply(stan.age.10.m,2,mean),
                 stan.025=apply(stan.age.10.m,2,quantile,probs=0.025),
                 stan.975=apply(stan.age.10.m,2,quantile,probs=0.975),
                 mgcv.mean=apply(age.10.m,2,mean),
                 mgcv.025=apply(age.10.m,2,quantile,probs=0.025),
                 mgcv.975=apply(age.10.m,2,quantile,probs=0.975),
                 tmb.mean=apply(tmb.age.10.m,2,mean),
                 tmb.025=apply(tmb.age.10.m,2,quantile,probs=0.025),
                 tmb.975=apply(tmb.age.10.m,2,quantile,probs=0.975)
  )
  
  gg.colors<-c("mgcv"="black","TMB"="blue","Stan"="red")
  
  ggplot(df)+
    geom_line(aes(x=10:50,y=mgcv.mean,color="mgcv"),size=1)+geom_ribbon(aes(x=10:50,ymin=mgcv.025,ymax=mgcv.975,fill="mgcv"),alpha=0.4)+
    geom_line(aes(x=10:50,y=stan.mean,color="TMB"),size=1)+geom_ribbon(aes(x=10:50,ymin=stan.025,ymax=stan.975,fill="Stan"),alpha=0.2)+
    geom_line(aes(x=10:50,y=tmb.mean,color="Stan"),size=1)+geom_ribbon(aes(x=10:50,ymin=tmb.025,ymax=tmb.975,fill="TMB"),alpha=0.3)+
    scale_color_manual(values=gg.colors,guide=FALSE)+scale_fill_manual(values=gg.colors,guide=FALSE)+
    ylab("Log Mortality")+xlab("Year")+ggtitle(paste(sex,"Year",i))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),legend.key.size=unit(0.5,"cm"),legend.key=element_rect(fill=NA,colour=NA,size=0.5),legend.position=c(0.8,0.2),plot.title = element_text(hjust = 0.5))
}

q.45.15<-function(x,DX){
  1-exp(apply(-exp(matrix(DX%*%x,age.end-age.start+1,year.end-year.start+1)[6:36,]),2,sum))
}

plot.q4515<-function(sex){
  if(sex=="Male"){
    mgcv.q.4515<-t(apply(param.sim,2,q.45.15,DX=Xm))
    stan.q.4515<-t(apply(stan.params.m,1,q.45.15,DX=te.spline))
    tmb.q.4515<-t(apply(tmb.sim[,19:(19+224)],1,q.45.15,DX=te.spline))
  } else{
    mgcv.q.4515<-t(apply(param.sim,2,q.45.15,DX=Xf))
    stan.q.4515<-t(apply(stan.params.f,1,q.45.15,DX=te.spline))
    tmb.q.4515<-t(apply(tmb.sim[,(19+225):(19+225+224)],1,q.45.15,DX=te.spline))
  }
  
  df<-data.frame(stan.mean=apply(stan.q.4515,2,mean),
                 stan.025=apply(stan.q.4515,2,quantile,probs=0.025),
                 stan.975=apply(stan.q.4515,2,quantile,probs=0.975),
                 mgcv.mean=apply(mgcv.q.4515,2,mean),
                 mgcv.025=apply(mgcv.q.4515,2,quantile,probs=0.025),
                 mgcv.975=apply(mgcv.q.4515,2,quantile,probs=0.975),
                 tmb.mean=apply(tmb.q.4515,2,mean),
                 tmb.025=apply(tmb.q.4515,2,quantile,probs=0.025),
                 tmb.975=apply(tmb.q.4515,2,quantile,probs=0.975)
  )
  
  gg.colors<-c("mgcv"="black","TMB"="blue","Stan"="red")
  
  ggplot(df)+
    geom_line(aes(x=1995:2014,y=mgcv.mean,color="mgcv"),size=1)+geom_ribbon(aes(x=1995:2014,ymin=mgcv.025,ymax=mgcv.975,fill="mgcv"),alpha=0.4)+
    geom_line(aes(x=1995:2014,y=stan.mean,color="TMB"),size=1)+geom_ribbon(aes(x=1995:2014,ymin=stan.025,ymax=stan.975,fill="Stan"),alpha=0.2)+
    geom_line(aes(x=1995:2014,y=tmb.mean,color="Stan"),size=1)+geom_ribbon(aes(x=1995:2014,ymin=tmb.025,ymax=tmb.975,fill="TMB"),alpha=0.3)+
    scale_color_manual(values=gg.colors,guide=FALSE)+scale_fill_manual(values=gg.colors,guide=FALSE)+
    ylab(bquote(""[45]*q[15]))+xlab("Year")+ggtitle(paste(sex))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),legend.key.size=unit(0.5,"cm"),legend.key=element_rect(fill=NA,colour=NA,size=0.5),legend.position=c(0.4,0.2),plot.title = element_text(hjust = 0.5))
}

library(gridExtra)
###selected ages
grid.arrange(plot.fun(10,sex="Male")+scale_color_manual(values=gg.colors,guide_legend(title="")),
             plot.fun(20,sex="Male"),plot.fun(30,sex="Male"),plot.fun(40,sex="Male"),plot.fun(50,sex="Male"),nrow=2)

grid.arrange(plot.fun(10,sex="Female"),plot.fun(20,sex="Female"),plot.fun(30,sex="Female"),plot.fun(40,sex="Female"),plot.fun(50,sex="Female"),nrow=2)
###selected years
grid.arrange(plot.fun.year(1995,sex="Male")+scale_color_manual(values=gg.colors,guide_legend(title="")),
             plot.fun.year(2000,sex="Male"),plot.fun.year(2005,sex="Male"),plot.fun.year(2010,sex="Male"),plot.fun.year(2014,sex="Male"),nrow=2)

grid.arrange(plot.fun.year(1995,sex="Female"),plot.fun.year(2000,sex="Female"),plot.fun.year(2005,sex="Female"),plot.fun.year(2010,sex="Female"),plot.fun.year(2014,sex="Female"),nrow=2)

###45q15
grid.arrange(plot.q4515("Male")+scale_color_manual(values=gg.colors,guide_legend(title="")),plot.q4515("Female"),nrow=1)

