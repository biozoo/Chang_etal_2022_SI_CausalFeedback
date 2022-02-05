### CHANG_ETAL_SUPPLEMENTARY_INFORMATION_RCODE - Full R Code
### Cross-system comparison among quantitative causal networks
### Updated in Feb 3, 2022
rm(list = ls())
library('rEDM') # Empirical dynamical modeling for CCM analysis
library('tidyr')
library('olsrr')
library('vegan')
library('gridExtra')
library('dplyr')
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

# Set working directory
setwd("D:\\data\\Meta_analysis\\Ecosystem network\\R_code_demo\\R_code_demo_20220203")
seed=25647
set.seed(seed)

####Function for correlation analysis between x:environment & y:causal links and feedbacks
cor.vab=function(y,x){
  cor.matrix=matrix(0,nrow=ncol(x),ncol=ncol(y))
  cor.matrix.p=matrix(1,nrow=ncol(x),ncol=ncol(y))
  for(i in 1:ncol(x)){
    for(j in 1:ncol(y)){
      cor.out=cor.test(x[,i],y[,j])
      cor.matrix[i,j]=cor.out$estimate
      cor.matrix.p[i,j]=cor.out$p.value
      colnames(cor.matrix)=colnames(cor.matrix.p)=colnames(y)
      rownames(cor.matrix)=rownames(cor.matrix.p)=colnames(x)
    }
  }
  return(list(estimate=cor.matrix,p.value=cor.matrix.p))
}

######################################################################
# Loading the dataset
# A dataset includes causal strengths derived from CCM 
istd=read.csv('CCMdata_Demo.csv',header=T)
# A dataset includes environmental characteristics 
# (Some of characteristics were obtained from longterm average of environmental factors)
indx.o = read.csv('LongtermMean_demo.csv',header=T)
indx=indx.o[,-1]
indx.s=apply(indx,2,scale)
# Select variable includes no missing data (Latitude was dropped because it highly correlated with temperature)
selectV=c('Richness','Biomass','Temp','NO3','PO4','Depth')
sys=as.character(indx.o[,1])


######################################################################################
# Data rearrangement: the twelve common causal links quantified by CCM for all systems
lik=c("Richness->Biomass","Richness->NO3","Richness->PO4",
      "NO3->Richness","NO3->Biomass","PO4->Richness",
      "PO4->Biomass","Biomass->Richness","Biomass->NO3",
      "Biomass->PO4",'Temp->Richness','Temp->Biomass')
lik.vab=matrix(unlist(strsplit(lik,'->')),ncol=2,byrow=T)
lik.ef.ind=c(1,5,7,12)
lik.bd.ind=c(4,6,8,11)
lik.name=c("BD->EF","BD->NO3","BD->PO4",
            "NO3->BD","NO3->EF","PO4->BD",
            "PO4->EF","EF->BD","EF->NO3",
            "EF->PO4",'Temp->BD','Temp->EF')

causalLK=causalLK.c=causalLK.sd=matrix(NA,length(sys),length(lik))
for(i in 1:length(lik)){
  causalLK[,i]=filter(istd,Cause==lik.vab[i,1]&Effect==lik.vab[i,2])[,'Std_L_strength']
  causalLK.c[,i]=filter(istd,Cause==lik.vab[i,1]&Effect==lik.vab[i,2])[,'Convergence']
  causalLK.sd[,i]=filter(istd,Cause==lik.vab[i,1]&Effect==lik.vab[i,2])[,'SD_std']
}
colnames(causalLK)=colnames(causalLK.c)=colnames(causalLK.sd)=lik.name
rownames(causalLK)=rownames(causalLK.c)=rownames(causalLK.sd)=sys
causalLK.c[is.na(causalLK.c)]=FALSE


####################################################
## The strength of causal links affecting Biomass
causalLK.ef=data.frame(sys,causalLK[,lik.ef.ind])
colnames(causalLK.ef)=c('system',lik.name[lik.ef.ind])
# Number of system with signifiacnt results
(out.significant.ef=apply(causalLK.c[,lik.ef.ind],2,sum))
# Number of system with rank 1 strength
(out.rank1.ef=table(lik.name[lik.ef.ind][apply(causalLK.ef[,-1],1,which.max)]))

## The strength of causal links affecting diversity
causalLK.bd=data.frame(sys,causalLK[,lik.bd.ind])
colnames(causalLK.bd)=c('system',lik.name[lik.bd.ind])
# Number of system with signifiacnt results
(out.significant.bd=apply(causalLK.c[,lik.bd.ind],2,sum))
# Number of system with rank 1 strength
(out.rank1.bd=table(lik.name[lik.bd.ind][apply(causalLK.bd[,-1],1,which.max)]))


#########################################################
### Violin plots of individual causal links  (Fig. 2a, b)
# Drivers of Biomass
my_datal.ef <- causalLK.ef %>% 
  dplyr::select('system',lik.name[lik.ef.ind]) %>%
  gather(lik.name[lik.ef.ind], 
         key = "Cause", value = "Std_L_strength")%>%
  arrange(Std_L_strength) %>%
  mutate(Cause = factor(Cause, levels=c("NO3->EF","BD->EF","PO4->EF","Temp->EF"))) 

g1=ggplot(my_datal.ef, aes(x=Cause, y=Std_L_strength)) + 
  geom_violin(trim=T, fill="gray")+
  ylim(0,1)+
  labs(title="",x="", y = "Standardized linkage strength")+
  geom_boxplot(width=0.1)+
  theme_classic()

# Drivers of diversity
my_datal.bd <- causalLK.bd %>% 
  dplyr::select('system',lik.name[lik.bd.ind]) %>%
  gather(lik.name[lik.bd.ind], 
         key = "Cause", value = "Std_L_strength")%>%
  arrange(Std_L_strength) %>%
  mutate(Cause = factor(Cause, levels=c("NO3->BD","EF->BD","PO4->BD","Temp->BD"))) 

g2=ggplot(my_datal.bd, aes(x=Cause, y=Std_L_strength)) + 
  geom_violin(trim=T, fill="gray")+
  ylim(0,1)+
  labs(title="",x="", y = "Standardized linkage strength")+
  geom_boxplot(width=0.1)+
  theme_classic()

windows(width = 30, height = 15)
pl = list(g1,g2)
margin = theme(plot.margin = unit(c(1,1,1,1), "cm"))
grid.arrange(grobs = lapply(pl, "+", margin),nrow=1,newpage = F)


#########################################################
########Bar plot strength list all systems (Fig. S3)
# Causal Effects on Biomass
#A function to add arrows on the chart
error.bar <- function(x, y, upper, lower=upper, rang=NULL, length=0.025,...){
  upp=y+upper;loo=y-lower;
  if(!is.null(rang)){upp[upp>rang[2]]=rang[2];loo[loo<rang[1]]=rang[1]}
  arrows(x,upp, x, loo, angle=90, code=3, length=length, ...)
}


win.graph(35,60)
par(mfcol=c(4,1),mar=c(4,4,1,1))
caca=c("orange","Yellow",'red',"black")
barmat=t(as.matrix(causalLK.ef[,-1]))[c("BD->EF","Temp->EF","NO3->EF","PO4->EF"),];colnames(barmat)=sys
barmat.sd=t(as.matrix(causalLK.sd))[c("BD->EF","Temp->EF","NO3->EF","PO4->EF"),];colnames(barmat.sd)=sys
ph1=barplot(barmat[,1:10],beside=T,col=caca,ylim=c(0,1.1),xlab="",ylab="Loop weight", main='Drivers of Biomass')
error.bar(ph1,barmat[,1:10], barmat.sd[,1:10],rang=c(0,1))
abline(h=0)
legend('topright',c("BD->EF","Temp->EF","NO3->EF","PO4->EF"),lty=1,col=caca,cex=0.75,lwd=5,bty='n')
ph2=barplot(barmat[,11:19],beside=T,col=caca,ylim=c(0,1.1),xlab="System",ylab="Loop weight")
error.bar(ph2,barmat[,11:19], barmat.sd[,11:19],rang=c(0,1))
abline(h=0)

# Causal Effects on Species Richness
caca=c("blue","Yellow",'red',"black")
barmat=t(as.matrix(causalLK.bd[,-1]))[c("EF->BD","Temp->BD","NO3->BD","PO4->BD"),];colnames(barmat)=sys
barmat.sd=t(as.matrix(causalLK.sd[,-1]))[c("EF->BD","Temp->BD","NO3->BD","PO4->BD"),];colnames(barmat.sd)=sys

ph1=barplot(barmat[,1:10],beside=T,col=caca,ylim=c(0,1.1),xlab="",ylab="Loop weight", main='Drivers of diversity')
error.bar(ph1,barmat[,1:10], barmat.sd[,1:10],rang=c(0,1))
abline(h=0)
legend('topright',c("EF->BD","Temp->BD","NO3->BD","PO4->BD"),lty=1,col=caca,cex=0.75,lwd=5,bty='n')
ph2=barplot(barmat[,11:19],beside=T,col=caca,ylim=c(0,1.1),xlab="System",ylab="Loop weight")
error.bar(ph2,barmat[,11:19], barmat.sd[,11:19],rang=c(0,1))
abline(h=0)




#############################
#Pairwise feedback
PFD.name=c('Richness<->Biomass','Richness<->NO3','Richness<->PO4','NO3<->Biomass','PO4<->Biomass')
Pair.name=matrix(unlist(strsplit(PFD.name,'<->')),ncol=2,byrow=T)
pfd=cpfd=dpfd=matrix(NA,nrow=length(sys),ncol=length(PFD.name))
for(i in 1:length(sys)){
  for(j in 1:nrow(Pair.name)){
    pfd1=filter(istd,Cause==Pair.name[j,1],Effect==Pair.name[j,2],system==sys[i])[,'Std_L_strength']
    pfd2=filter(istd,Cause==Pair.name[j,2],Effect==Pair.name[j,1],system==sys[i])[,'Std_L_strength']
    pfd.c1=filter(istd,Cause==Pair.name[j,1],Effect==Pair.name[j,2],system==sys[i])[,'Convergence']
    pfd.c2=filter(istd,Cause==Pair.name[j,2],Effect==Pair.name[j,1],system==sys[i])[,'Convergence']
    
    pfd[i,j]=sqrt(pfd1*pfd2)
    dpfd[i,j]=pfd1-pfd2
    cpfd[i,j]=pfd.c1*pfd.c2
  }
}

PFD.name2=c('BD<->EF','BD<->NO3','BD<->PO4','NO3<->EF','PO4<->EF')
cpfd[is.na(cpfd)]=0  
colnames(cpfd)=PFD.name2

pfd=data.frame(system=sys,pfd)
dpfd=data.frame(system=sys,dpfd)
colnames(pfd)=colnames(dpfd)=c('system',PFD.name2)
pfd.c=pfd;pfd.c[,-1]=pfd.c[,-1]*cpfd
pfd.c[pfd.c==0]=NA

# Number of system with signifiacnt results
(out.significant.pfd=apply(cpfd,2,sum))
# Number of system with rank 1 strength
(out.rank1.bd=table(PFD.name2[unlist(apply(pfd.c[,-1],1,which.max))]))


########################
### Triangular feedbacks
TFD.name=rbind(
 'Richness->Biomass->NO3->Richness',
 'Biomass->Richness->NO3->Biomass',
 'Richness->Biomass->PO4->Richness',
 'Biomass->Richness->PO4->Biomass'
)
TFD.cause=matrix(unlist(lapply(strsplit(TFD.name,'->'),function(x){return(x[1:3])})),ncol=3,byrow=T)
TFD.effect=matrix(unlist(lapply(strsplit(TFD.name,'->'),function(x){return(x[2:4])})),ncol=3,byrow=T)

TFD.name2=c('Type I-N', 'Type II-N' ,'Type I-P', 'Type II-P')

trifd=ctfd=matrix(NA,nrow=length(sys),ncol=length(TFD.name))
for(i in 1:length(sys)){
  for(j in 1:length(TFD.name)){
    tfd1=filter(istd,Cause==TFD.cause[j,1],Effect==TFD.effect[j,1],system==sys[i])[,'Std_L_strength']
    tfd2=filter(istd,Cause==TFD.cause[j,2],Effect==TFD.effect[j,2],system==sys[i])[,'Std_L_strength']
    tfd3=filter(istd,Cause==TFD.cause[j,3],Effect==TFD.effect[j,3],system==sys[i])[,'Std_L_strength']
    
    tfd.c1=filter(istd,Cause==TFD.cause[j,1],Effect==TFD.effect[j,1],system==sys[i])[,'Convergence']
    tfd.c2=filter(istd,Cause==TFD.cause[j,2],Effect==TFD.effect[j,2],system==sys[i])[,'Convergence']
    tfd.c3=filter(istd,Cause==TFD.cause[j,3],Effect==TFD.effect[j,3],system==sys[i])[,'Convergence']
    

    trifd[i,j]=(tfd1*tfd2*tfd3)^(1/3)
    ctfd[i,j]=(tfd.c1*tfd.c2*tfd.c3)
  }
}

ctfd[is.na(ctfd)]=0
colnames(ctfd)=TFD.name2
trifd=data.frame(sys,trifd)
colnames(trifd)=c('system',TFD.name2)
trifd.c=trifd;trifd.c[,-1]=trifd.c[,-1]*ctfd
trifd.c[trifd.c==0]=NA

# Number of system with signifiacnt results
(out.significant.tfd=apply(ctfd,2,sum))

# Number of system with at least one signifiacnt triangular feedback
sum(apply(ctfd==1,1,any))

# Number of system with rank 1 strength
(out.rank1.tfd=table(TFD.name2[unlist(apply(trifd.c[,-1],1,which.max))]))


#########################################################
### Violin plots for more complex feedback  (Fig. 2c, d)
data.pfd <- pfd %>% 
  dplyr::select('system',PFD.name2) %>%
  gather(PFD.name2, 
         key = "Cause", value = "Std_L_strength")%>%
  arrange(Std_L_strength) %>%
  mutate(Cause = factor(Cause, levels=c("BD<->EF","BD<->NO3","BD<->PO4","NO3<->EF","PO4<->EF"))) 

fd1=ggplot(data.pfd, aes(x=Cause, y=Std_L_strength)) + 
  geom_violin(trim=T, fill="gray")+
  ylim(0,1)+
  labs(title="",x="", y = "Loop weight")+
  geom_boxplot(width=0.1)+
  theme_classic()


data.tfd <- trifd %>% 
  dplyr::select('system',TFD.name2) %>%
  gather(TFD.name2, 
         key = "Cause", value = "Std_L_strength")%>%
  arrange(Std_L_strength) %>%
  mutate(Cause = factor(Cause, levels=c("Type I-N","Type II-N","Type I-P","Type II-P"))) 

fd2=ggplot(data.tfd, aes(x=Cause, y=Std_L_strength)) + 
  geom_violin(trim=T, fill="gray")+
  ylim(0,1)+
  labs(title="",x="", y = "Loop weight")+
  geom_boxplot(width=0.1)+
  theme_classic()

################################################################
## Directional bias in pairwise feedbacks  Fig. 3
windows(width = 30, height = 15)
pl = list(fd1,fd2)
margin = theme(plot.margin = unit(c(1,1,1,1), "cm"))
grid.arrange(grobs = lapply(pl, "+", margin),nrow=1,newpage = F)

##################################################
# Strength differences in simple pairwise feedback 
my_datal <- dpfd %>% 
  dplyr::select('system',PFD.name2) %>%
  gather(PFD.name2, 
         key = "Cause", value = "Std_L_strength")%>%
    arrange(Std_L_strength) %>%
    mutate(Cause = factor(Cause, levels=c('BD<->EF','BD<->NO3','BD<->PO4','NO3<->EF','PO4<->EF'))) 

pfd.d=ggplot(my_datal, aes(x=Cause, y=Std_L_strength)) + 
  geom_violin(trim=T, fill="gray")+
  ylim(-0.7,0.7)+
  geom_hline(yintercept = 0, linetype=2)+
  labs(title="",x="", y = "Difference in linkage strength ")+
  geom_boxplot(width=0.1)+
  theme_classic()

windows()
pfd.d
# the number of systems showing positive difference in linkage strength ('->' is stronger than '<-')
apply(dpfd[,-1]>0,2,sum)




#################################################################################
### Multiple regression based on forward variable selection 
# Environmental factors affecting the strength of diversity effects on Biomass
X=indx.s[,selectV]
Y=causalLK[,'BD->EF']
selec.v.ef=NULL
model <- lm(Y ~ ., data = data.frame(X))
select.v=ols_step_forward_p(model, penter =0.1)
(selec.v.ef=union(selec.v.ef,select.v$predictors))
model2 <- lm(Y ~., data = data.frame(X[,select.v$predictors]))
summary(model2)

# Environmental factors affecting the strength of Biomass feedback on diversity 
Y=causalLK[,'EF->BD']
selec.v.bd=NULL
model <- lm(Y ~ ., data = data.frame(X))
select.v=ols_step_forward_p(model, penter =0.1)
selec.v.bd=union(selec.v.bd,select.v$predictors)
model2 <- lm(Y ~., data = data.frame(X[,select.v$predictors]))
summary(model2)


#####################################
# Direct gradient analysis
#####################################
#####################################
# RDA analysis of Biomass drivers
X.o=indx[,selectV]
X=indx.s[,selectV]

Y.o=causalLK[,lik.name[lik.ef.ind]]
Y=decostand(Y.o,'normalize')
sord=metaMDS(Y)$points

sig.v=names(which(apply(cor.vab(Y.o,X.o)$p.value<=0.05,1,sum)>=(0.5*ncol(Y.o))))
model.1 <- lm(sord[,1] ~ ., data = data.frame(X))
select.v1=ols_step_forward_p(model.1)
model.2 <- lm(sord[,2] ~ ., data = data.frame(X))
select.v2=ols_step_forward_p(model.2)
(select.v3=sort(unique(c(select.v1$predictors,select.v2$predictors,sig.v))))

rda.out=rda(Y,X[,select.v3])
#Data variations explained by RDA (%)
(Variex.r.ef=round(summary(rda.out)$cont[[1]][2,1:2]*100,1))
permutest(rda.out,permutations = how(nperm=10000))
da.rda.ef=rda.out$CCA

###################################
# RDA analysis of diversity drivers
Y.o=causalLK[,lik.name[lik.bd.ind]]
Y=decostand(Y.o,'normalize')

sord=metaMDS(Y)$points
sig.v=names(which(apply(cor.vab(Y.o,X.o)$p.value<=0.05,1,sum)>=(0.5*ncol(Y.o))))
model.1 <- lm(sord[,1] ~ ., data = data.frame(X))
select.v1=ols_step_forward_p(model.1)
model.2 <- lm(sord[,2] ~ ., data = data.frame(X))
select.v2=ols_step_forward_p(model.2)
(select.v3=sort(unique(c(select.v1$predictors,select.v2$predictors,sig.v))))


rda.out=rda(Y,X[,select.v3])
#Data variations explained by RDA (%)
(Variex.r.bd=round(summary(rda.out)$cont[[1]][2,1:2]*100,1))
permutest(rda.out,permutations = how(nperm=10000))
da.rda.bd=rda.out$CCA

###################################
# RDA analysis of pairwise feedback
Y=decostand(pfd[,-1],'normalize')
Y.o=pfd[,-1]
sord=metaMDS(Y)$points

sig.v=names(which(apply(cor.vab(Y.o,X.o)$p.value<=0.05,1,sum)>=(0.5*ncol(Y.o))))
model.1 <- lm(sord[,1] ~ ., data = data.frame(X))
select.v1=ols_step_forward_p(model.1)
model.2 <- lm(sord[,2] ~ ., data = data.frame(X))
select.v2=ols_step_forward_p(model.2)
(select.v3=sort(unique(c(select.v1$predictors,select.v2$predictors,sig.v))))

rda.out=rda(Y,X[,select.v3])
#Data variations explained by RDA (%)
(Variex.r.pfd=round(summary(rda.out)$cont[[1]][2,1:2]*100,1))
permutest(rda.out,permutations = how(nperm=10000))
da.rda.pfd=rda.out$CCA

######################################
#  RDA analysis of triangular feedback
Y=decostand(trifd[,-1],'normalize')
Y.o=trifd[,-1]
sig.v=names(which(apply(cor.vab(Y.o,X.o)$p.value<=0.05,1,sum)>=(0.5*ncol(Y.o))))
sord=metaMDS(Y)$points

model.1 <- lm(sord[,1] ~ ., data = data.frame(X))
select.v1=ols_step_forward_p(model.1)
model.2 <- lm(sord[,2] ~ ., data = data.frame(X))
select.v2=ols_step_forward_p(model.2)
(select.v3=sort(unique(c(select.v1$predictors,select.v2$predictors,sig.v))))

rda.out=rda(Y,X[,select.v3])
#Data variations explained by RDA (%)
(Variex.r.tfd=round(summary(rda.out)$cont[[1]][2,1:2]*100,1))
permutest(rda.out,permutations = how(nperm=10000))
da.rda.tfd=rda.out$CCA

#################################################################################
# RDA biplots with color scales based on longterm average of production (Fig. 4)
da.rda=da.rda.ef
mag.f=1.05
mag.f2=0.8
g.ef=ggplot(data.frame(da.rda$wa,indx), aes(x = RDA1, y = RDA2, colour = log(Biomass,10))) +
  geom_point(size = 3, shape=19) +
  scale_color_gradient(low="green", high="darkgreen")+
  labs(tag='a', x = paste("RDA 1 (", Variex.r.ef[1], '%)'), y = paste("RDA 2 (", Variex.r.ef[2], '%)'))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) +
  geom_vline(xintercept = 0, colour = 'grey') + 
  geom_hline(yintercept = 0, colour = 'grey') +
  geom_segment(data=data.frame(da.rda$biplot*mag.f2), aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="blue", arrow=arrow(length=unit(0.15,"inches"))) +
  geom_segment(data=data.frame(da.rda$v*mag.f2), aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="red", arrow=arrow(length=unit(0,"inches"))) +
  annotate(geom="text", x=da.rda$wa[,1]*mag.f, y=da.rda$wa[,2]*mag.f, label=sys,color="black")+
  annotate(geom="text", x=da.rda$v[,1]*mag.f2, y=da.rda$v[,2]*mag.f2, label=rownames(da.rda$v),color="red")+
  annotate(geom="text", x=da.rda$biplot[,1]*mag.f2, y=da.rda$biplot[,2]*mag.f2, label=rownames(da.rda$biplot),color="blue")

da.rda=da.rda.bd
g.bd=ggplot(data.frame(da.rda$wa,indx), aes(x = RDA1, y = RDA2, colour = log(Biomass,10))) +
  geom_point(size = 3, shape=19) +
  scale_color_gradient(low="green", high="darkgreen")+
  labs(tag='b', x = paste("RDA 1 (", Variex.r.bd[1], '%)'), y = paste("RDA 2 (", Variex.r.bd[2], '%)'))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) +
  geom_vline(xintercept = 0, colour = 'grey') + 
  geom_hline(yintercept = 0, colour = 'grey') +
  geom_segment(data=data.frame(da.rda$biplot*mag.f2), aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="blue", arrow=arrow(length=unit(0.15,"inches"))) +
  geom_segment(data=data.frame(da.rda$v*mag.f2), aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="red", arrow=arrow(length=unit(0,"inches"))) +
  annotate(geom="text", x=da.rda$wa[,1]*mag.f, y=da.rda$wa[,2]*mag.f, label=sys,color="black")+
  annotate(geom="text", x=da.rda$v[,1]*mag.f2, y=da.rda$v[,2]*mag.f2, label=rownames(da.rda$v),color="red")+
  annotate(geom="text", x=da.rda$biplot[,1]*mag.f2, y=da.rda$biplot[,2]*mag.f2, label=rownames(da.rda$biplot),color="blue")

da.rda=da.rda.pfd;
# For consistency with other RDA plots, the direction of gradient changed
da.rda$wa=-da.rda$wa;da.rda$v=-da.rda$v;da.rda$biplot=-da.rda$biplot;
g.pfd=ggplot(data.frame(da.rda$wa,indx), aes(x = RDA1, y = RDA2, colour = log(Biomass,10))) +
  geom_point(size = 3, shape=19) +
  scale_color_gradient(low="green", high="darkgreen")+
  labs(tag='c', x = paste("RDA 1 (", Variex.r.pfd[1], '%)'), y = paste("RDA 2 (", Variex.r.pfd[2], '%)'))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) +
  geom_vline(xintercept = 0, colour = 'grey') + 
  geom_hline(yintercept = 0, colour = 'grey') +
  geom_segment(data=data.frame(da.rda$biplot*mag.f2), aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="blue", arrow=arrow(length=unit(0.15,"inches"))) +
  geom_segment(data=data.frame(da.rda$v*mag.f2), aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="red", arrow=arrow(length=unit(0,"inches"))) +
  annotate(geom="text", x=da.rda$wa[,1]*mag.f, y=da.rda$wa[,2]*mag.f, label=sys,color="black")+
  annotate(geom="text", x=da.rda$v[,1]*mag.f2, y=da.rda$v[,2]*mag.f2, label=rownames(da.rda$v),color="red")+
  annotate(geom="text", x=da.rda$biplot[,1]*mag.f2, y=da.rda$biplot[,2]*mag.f2, label=rownames(da.rda$biplot),color="blue")

da.rda=da.rda.tfd
g.tfd=ggplot(data.frame(da.rda$wa,indx), aes(x = RDA1, y = RDA2, colour = log(Biomass,10))) +
  geom_point(size = 3, shape=19) +
  scale_color_gradient(low="green", high="darkgreen")+
  labs(tag='d', x = paste("RDA 1 (", Variex.r.tfd[1], '%)'), y = paste("RDA 2 (", Variex.r.tfd[2], '%)'))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(), plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1)) +
  geom_vline(xintercept = 0, colour = 'grey') + 
  geom_hline(yintercept = 0, colour = 'grey') +
  geom_segment(data=data.frame(da.rda$biplot*mag.f2), aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="blue", arrow=arrow(length=unit(0.15,"inches"))) +
  geom_segment(data=data.frame(da.rda$v*mag.f2), aes(x=0, xend=RDA1, y=0, yend=RDA2), 
               color="red", arrow=arrow(length=unit(0,"inches"))) +
  annotate(geom="text", x=da.rda$wa[,1]*mag.f, y=da.rda$wa[,2]*mag.f, label=sys,color="black")+
  annotate(geom="text", x=da.rda$v[,1]*mag.f2, y=da.rda$v[,2]*mag.f2, label=rownames(da.rda$v),color="red")+
  annotate(geom="text", x=da.rda$biplot[,1]*mag.f2, y=da.rda$biplot[,2]*mag.f2, label=rownames(da.rda$biplot),color="blue")

windows(width = 30, height = 22.5)
grid.arrange(g.ef, g.bd, g.pfd, g.tfd, nrow = 2)
