### CHANG_ETAL_SUPPLEMENTARY_INFORMATION_RCODE - Full R Code
### Quantification of causal links by CCM analysis
### Updated in Feb 3, 2022
rm(list = ls())
# Empirical dynamical modeling for CCM analysis was based on 'rEDM' package (v1.2.3) that is available on https://github.com/cran/rEDM/releases/tag/1.2.3
library('rEDM') 
library('Kendall') # Kendall's tau test for the convergence of CCM
library('tidyr')
library('dplyr')

###########################
# Set path for CCM analyses
setwd("D:\\data\\Meta_analysis\\Ecosystem network\\R_code_demo\\R_code_demo_20211025") 
seed=25647
set.seed(seed)

####Function for time series standardization (normalization, detrend and deseason)
# x: time series data (vector form)
# normalization: normalization to zero mean & unit variance 
# dseason: logic argument for Deseasonalization by monthly mean
# season_sd: Deseasonalization by both monthly mean and monthly standard deviation
# sea: The period of seasonality (e.g. 12 months for monthly data)
# dtrend: logic argument for detrend
# dTtype: Type of detrend by first difference (first) or linear regression (linear)
nomz=function(x, normalization=T, dseason=T, season_sd=T, sea=12, dtrend=T, dTtype="linear"){
  x=as.numeric(x)
  xt=x
  # Detrend
  if(dtrend==T & dTtype=="first"){xt=diff(xt)} else if (dtrend==T & dTtype=="linear"){
    lm.t=lm(xt~c(1:length(xt)))
    xt=xt-(lm.t$coefficients[1]+lm.t$coefficients[2]*c(1:length(xt)))}
  # Deseason
  if(dseason==T){
    xs=as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,mean,na.rm=T))
    xsd=as.numeric(apply(matrix(xt[1:(sea*length(xt)%/%sea)],ncol=sea,byrow=T),2,sd,na.rm=T))
    xt=xt-c(rep(xs,1+length(xt)%/%sea))[1:length(xt)]
    if(season_sd==T){xt=xt/(c(rep(xsd,1+length(xt)%/%sea))[1:length(xt)])}}
  # Normalization (zero mean & unity variance)
  if(normalization==T){xt=(xt-mean(xt,na.rm=T))/sd(xt,na.rm=T)}
  return(xt)
  
}

####### Function for generating lag time series 
laf=function(x,y,lagf){
  n <- NROW(x)
  x.t=x;y.t=y
  if(lagf<=0){x.t=x.t[(1-lagf):n];y.t=y.t[1:(n+lagf)]} # if lagf<0, y is leading
  if(lagf>0){x.t=x.t[1:(n-lagf)];y.t=y.t[(1+lagf):n]}  # if lagf>0, x is leading           
  return(cbind(x.t,y.t))
}

# data & function loading ends

# CCM causality network
######################################################################################
###### The reconstruction of causality networks ######################################
######################################################################################
# Loading time series data for Lake Kasumigaura station 9 (Ks9) 
# The example data is obtained from the Lake Kasumigaura Database (https://db.cger.nies.go.jp/gem/moni-e/inter/GEMS/database/kasumi/index.html). 
# The users should strictly follow the terms of use (https://db.cger.nies.go.jp/gem/moni-e/inter/GEMS/database/kasumi/contents/terms.html) 
# and contact the database administrator (biodiv.data@nies.go.jp)

month.dat=read.csv("Demo_Ks9_tsdata.csv",header=T,fill=T)
n=nrow(month.dat) # time series length
vb=c('Year','Month','Richness','Biomass','Temp','NO3','PO4')# rename the variables
colnames(month.dat)=vb

istd=NULL

# Detrend + deseason time series
sdat=data.frame(month.dat[,c(1:2)],apply(month.dat[,-c(1:2)],2,nomz))

##The index for testing causal links (a total of 12 links)
indmat=matrix(0,12,2);colnames(indmat)=c('Effect','Cause')
indmat[,1]=c('Biomass','NO3','PO4','Richness','NO3','PO4','Richness','Biomass','Richness','Biomass','Richness','Biomass')
indmat[,2]=c('Richness','Richness','Richness','Biomass','Biomass','Biomass','Temp','Temp','NO3','NO3','PO4','PO4')
indmat

# Determine the embedding dimensions (En) in CCM by try-and-error with best hindcast (tp=-1) skill
Emax=20
En=NULL
for(i in 1:nrow(indmat)){
  E.test=NULL
  for(E.t in 2:Emax){
    cmxy.t <- ccm(sdat, E = E.t,
                  lib_column = indmat[i,1], target_column = indmat[i,2],
                  lib_sizes = n, tp=-1,random_libs = F) 
    E.test=c(E.test,mean(cmxy.t$rho))
  }
  # Select the embedding dimension that makes the model with the highest hindcast predictive skill
  En=c(En,which.max(E.test)+1) 
}


################################################################################
### CCM analysis for all causality testlinks
lib_siz=sort(c(5,10,20,30,40,seq(50,n,50),n)) # a sequence of library size
#lib_siz=sort(c(5,10,20,30,40,seq(50,n,50),184,292,n)) # a sequence of library size
ccmda=uncertain.i=NULL
rho.i=NULL
for(i in 1:nrow(indmat)){
  ccmda.t=NULL
  
  for(j in 0:-3){
    da.t=laf(sdat[,indmat[i,1]],sdat[,indmat[i,2]],lagf=j) # Varying time lags
    colnames(da.t)=indmat[i,]
    # CCM analysis cross-mapping from one effect variable to its cause
    x_xmap_y <- ccm(da.t, E = En[i], # The embedding dimension E for each link were determined in previous step
                    lib_column = indmat[i,'Effect'], target_column = indmat[i,'Cause'],
                    lib_sizes = lib_siz, tp=0,RNGseed = seed,
                    num_samples = 100,replace=F)
    
    # Take average for the predictive skill under each library size
    aveg=cbind(unique(x_xmap_y$lib_size),aggregate(x_xmap_y[,c('rho')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'],
               aggregate(x_xmap_y[,c('mae')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'],
               aggregate(x_xmap_y[,c('rmse')], by=list(as.factor(x_xmap_y$lib_size)), mean)[,'x'])
    ccm_mean=data.frame(lag=rep(j,nrow(aveg)),x_xmap_y[1:nrow(aveg),]);
    ccm_mean[,c('lib_size','rho','mae','rmse')]=aveg
    ccm_mean[ccm_mean[,'rho']<0,'rho']=0
    
    ###########################
    # Convergence test in CCM
    ###########################
    # Fisher's delta rho Z test
    rho.Lmax=ccm_mean$rho[which.max(ccm_mean$lib_size)]
    rho.Lmin=ccm_mean$rho[1]
    ns=min(sum(!is.na(sdat[,indmat[i,1]])),sum(!is.na(sdat[,indmat[i,2]])))
    delta_rho=rho.Lmax-rho.Lmin
    z=abs(0.5*(log((1+rho.Lmax)/(1-rho.Lmax))-log((1+rho.Lmin)/(1-rho.Lmin)))*(2/(ns-3))^-0.5)
    z.p=(1-pnorm(z))
    # Kendall's tau test
    if(length(ccm_mean$rho)>3){
      kend=MannKendall(ccm_mean$rho)
      kend.tau=kend$tau[1]
      kend.p=kend$sl[[1]]
    }else{
      kend.tau=NA
      kend.p=NA
    }
    
    # Compile all the testing results
    ccmda.t=rbind(ccmda.t,
                  unlist(c(ccm_mean[nrow(ccm_mean),c(1:5,8)],rho_Lmax=rho.Lmax,rho_Lmin=rho.Lmin,
                           Z=z,p_Z=z.p,Kendall_tau=kend.tau,Kendall_p=kend.p)))
  }
  # Select the CCM results based on predictive skills
  lag_ind=which.max(ccmda.t[,'rho_Lmax'])
  best.lag=as.numeric(ccmda.t[lag_ind,'lag'])
  ccmda.olag=ccmda.t[lag_ind,]
  da.t=laf(sdat[,indmat[i,1]],sdat[,indmat[i,2]],lagf=best.lag) # Varying time lags
  colnames(da.t)=indmat[i,]
  x_xmap_y2 <- ccm(da.t, E = En[i], # The embedding dimension E for each link were determined in previous step
                   lib_column = indmat[i,'Effect'], target_column = indmat[i,'Cause'],
                   lib_sizes = ccmda.olag['lib_size'], tp=0,RNGseed = seed,
                   num_samples = 500,replace=T)
  
  rho.it=x_xmap_y2[,'rho'];rho.it[rho.it<0]=0
  rho.i=cbind(rho.i,rho.it)
  uncertain.i=rbind(uncertain.i,c(lag=best.lag,SD=sd(rho.it,na.rm=T)))
  ccmda=rbind(ccmda,ccmda.olag)
  cat("\r",i/nrow(indmat)*100,"%")
}
ccmda=data.frame(indmat,ccmda)
rownames(ccmda)=NULL
uncertain.i=data.frame(system=rep('Ks9',nrow(indmat)),indmat,uncertain.i)
Convergence=ccmda$p_Z<0.05 & ccmda$Kendall_p<=0.05 & ccmda$Kendall_tau>0
ccmda=data.frame(ccmda,Convergence)

# Standardized linkage strength by dividing the maximal linkage strength within the a system
(linkM.std=data.frame(system=rep('Ks9',nrow(indmat)),ccmda[,c('Cause','Effect','rho_Lmax','p_Z','Kendall_tau','Kendall_p','Convergence')],
                     Std_L_strength=ccmda$rho_Lmax/max(ccmda$rho_Lmax[1:12]),
                     SD=uncertain.i[,'SD'],SDstd=uncertain.i[,'SD']/max(ccmda$rho_Lmax[1:12])))
write.csv(linkM.std,'ccmda_Ks9_demo.csv',row.names=F)

############################################################################################
###  By repeating all these analyses to reconstruct causality networks for all systems 
###  Then we carried out the following cross-system comparison.
###  The dataset from the other systems are available as explained in Supplementary Table S2 
############################################################################################
