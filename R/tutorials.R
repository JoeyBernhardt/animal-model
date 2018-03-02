
library(tidyverse)
library(MCMCglmm)


					 

# read in data and get it in order  ------------------------------------------------------------
Data<-as.data.frame(read.table(file= "data-raw/gryphon.txt",header=TRUE))
names(Data)[1]<-"animal"


Data$animal<-as.factor(Data$animal)
Data$MOTHER<-as.factor(Data$MOTHER)
Data$BYEAR<-as.factor(Data$BYEAR)
Data$SEX<-as.factor(Data$SEX)
Data$BWT<-as.numeric(Data$BWT)
Data$TARSUS<-as.numeric(Data$TARSUS)
head(Data)


Ped <- as.data.frame(read.table(file= "data-raw/gryphonped.txt",header=TRUE))
for(x in 1:3) Ped[,x]<-as.factor(Ped[,x])
head(Ped)


## specify priors for the animal effect and the residual variance
p.var <- var(Data$BWT,na.rm=TRUE)
prior1.1<-list(G=list(G1=list(V=matrix(p.var/2),n=1)),
			   R=list(V=matrix(p.var/2),n=1))

model1.1_ped <- MCMCglmm(BWT ~ 1, pedigree = Ped, random = ~ animal, data = Data, prior = prior1.1)
model1.1_no_ped <- MCMCglmm(BWT ~ 1, random = ~ animal, data = Data, prior = prior1.1)

plot(model1.1_ped$Sol)
plot(model1.1_ped$VCV)

plot(model1.1_no_ped$Sol)
plot(model1.1_no_ped$VCV)

plot(model1.1b$Sol)
plot(model1.1b$VCV)

 model1.1b<-MCMCglmm(BWT~1,random=~animal,
					 	pedigree=Ped,data=Data,
					 nitt=65000,thin=50,burnin=15000,
					 prior=prior1.1,verbose=FALSE)
 
 autocorr(model1.1b$VCV)
 
 posterior.mode(model1.1b$VCV)
 
 HPDinterval(model1.1b$VCV)
 
prior1.1.2<-list(G=list(G1=list(V=matrix(p.var*0.95),n=1)), R=list(V=matrix(p.var*0.05),n=1))
model1.1.2<-MCMCglmm(BWT~1,random=~animal,pedigree=Ped,data=Data,prior=prior1.1.2,
					 nitt=65000,thin=50,burnin=15000,verbose=FALSE)


posterior.mode(model1.1b$VCV)

posterior.mode(model1.1.2$VCV)


posterior.heritability1.1 <- model1.1b$VCV[,"animal"]/(model1.1b$VCV[,"animal"]+model1.1b$VCV[,"units"])
HPDinterval(posterior.heritability1.1,0.95)

posterior.mode(posterior.heritability1.1)
plot(posterior.heritability1.1)


# bivariate animal model --------------------------------------------------


phen.var<-matrix(c(var(Data$BWT,na.rm=TRUE),0,0, var(Data$TARSUS,na.rm=TRUE)),2,2)
prior2.1<-list(G=list(G1=list(V=phen.var/2,n=2)), R=list(V=phen.var/2,n=2))

model2.1 <- MCMCglmm(cbind(BWT,TARSUS)~trait-1,
	random=~us(trait):animal,
	rcov=~us(trait):units,
	family=c("gaussian","gaussian"),
	pedigree=Ped,data=Data,
	prior=prior2.1,verbose=FALSE)

plot(model2.1$VCV[,"traitTARSUS:traitTARSUS.animal"])

model2.1$VCV

# model2.1<-dget(file="~/Desktop/JAE_MCMCglmm/model2point1LongRun.Rdat") > # the autocorrelation of the genetic variance of TARSUS at Lag 5
autocorr(model2.1$VCV)[,,"animal.trait.TARSUS"][3,4]
 
posterior.mode(model2.1$VCV)


# working with Pierre de Villemereuil’s tutorial --------------------------


## read in data
pedigreemulti<-read.table('data-raw/pedigreemulti.txt',header=T)
datamulti<-read.table('data-raw/datamulti.txt',header=T)
head(datamulti)

## define the prior
prior<- list(R=list(V=diag(2)/2,nu=2),
		 G=list(G1=list(V=diag(2)/2,nu=2)))


modelmulti <- MCMCglmm(cbind(phen1,phen2)~trait-1, ## the notation trait-1 allows us to estimate the mean of each trait, instead of the contrast between the two traits
					 random=~us(trait):animal, ##  define the structure of the variance-covariance matrix for the random effects 
					 rcov=~us(trait):units,  ##  define the structure of the variance-covariance matrix for the residual variances
					 family=c("gaussian","gaussian"), ## family requires a vector with the data distribution of each trait
					 prior=prior,
					 pedigree=pedigreemulti, ## specifies the pedigree; one row for each individual, each column specifies the parents for each individual
					 data=datamulti, ## includes two traits for each individual, i.e. phen1 and phen 2
					 nitt=10000,
					 burnin=1000,
					 thin=10)

summary(modelmulti)

head(modelmulti$VCV)

herit1<-modelmulti$VCV[,'traitphen1:traitphen1.animal']/ (modelmulti$VCV[,'traitphen1:traitphen1.animal']+ modelmulti$VCV[,'traitphen1:traitphen1.units'])
herit2<-modelmulti$VCV[,'traitphen2:traitphen2.animal']/(modelmulti$VCV[,'traitphen2:traitphen2.animal']+modelmulti$VCV[,'traitphen2:traitphen2.units'])
mean(herit1)
mean(herit2)


corr.gen<-modelmulti$VCV[,'traitphen1:traitphen2.animal']/
	sqrt(modelmulti$VCV[,'traitphen1:traitphen1.animal']*modelmulti$VCV[,'traitphen2:traitphen2.animal'])
	mean(corr.gen)

## now run the same thing, but just don't give it a pedigree (i.e. remove that argument altogether)

modelmulti_no_ped <- MCMCglmm(cbind(phen1,phen2)~trait-1, ## the notation trait-1 allows us to estimate the mean of each trait, instead of the contrast between the two traits
					   random=~us(trait):animal, ##  define the structure of the variance-covariance matrix for the random effects 
					   rcov=~us(trait):units,  ##  define the structure of the variance-covariance matrix for the residual variances
					   family=c("gaussian","gaussian"), ## family requires a vector with the data distribution of each trait
					   prior=prior,
					   # pedigree=pedigreemulti, ## specifies the pedigree; one row for each individual, each column specifies the parents for each individual
					   data=datamulti, ## includes two traits for each individual, i.e. phen1 and phen 2
					   nitt=10000,
					   burnin=1000,
					   thin=10)



summary(modelmulti_no_ped)


herit1_no_ped<-modelmulti_no_ped$VCV[,'traitphen1:traitphen1.animal']/ (modelmulti_no_ped$VCV[,'traitphen1:traitphen1.animal']+ modelmulti_no_ped$VCV[,'traitphen1:traitphen1.units'])
herit2_no_ped<-modelmulti_no_ped$VCV[,'traitphen2:traitphen2.animal']/(modelmulti_no_ped$VCV[,'traitphen2:traitphen2.animal']+modelmulti_no_ped$VCV[,'traitphen2:traitphen2.units'])
mean(herit1_no_ped)
mean(herit2_no_ped)


## now run the same thing, but give it a fake, clonal pedigree, with all the parents having NA's

pedrigree_fake <- pedigreemulti
pedrigree_fake$sire <- seq(from = 1001, to = 2000, by = 1 )
pedrigree_fake$dam <- seq(from = 1501, to = 2500, by = 1 )
pedrigree_fake$dam <- NA

### hmm, MCMCglmm doesn't seem to like any of my fake pedigrees
modelmulti_fake <- MCMCglmm(cbind(phen1,phen2)~trait-1, ## the notation trait-1 allows us to estimate the mean of each trait, instead of the contrast between the two traits
							  random=~us(trait):animal, ##  define the structure of the variance-covariance matrix for the random effects 
							  rcov=~us(trait):units,  ##  define the structure of the variance-covariance matrix for the residual variances
							  family=c("gaussian","gaussian"), ## family requires a vector with the data distribution of each trait
							  prior=prior,
							  pedigree=pedrigree_fake, ## specifies the pedigree; one row for each individual, each column specifies the parents for each individual
							  data=datamulti, ## includes two traits for each individual, i.e. phen1 and phen 2
							  nitt=10000,
							  burnin=1000,
							  thin=10)
?MCMCglmm()


### Lynn's code
mod <-MCMCglmm(cbind(AgeMat, HMI, PI, SizeMat, SizeNeo, SpineSizeNeo, Fec) ~ trait - 1, 
			   data = full.data, 
			   prior = prior,
			   random=~us(trait):Clone,
			   rcov=~us(trait):units,
			   family=c("gaussian","gaussian","gaussian","gaussian","gaussian","gaussian","gaussian"),
			   nitt=500000,burnin=50000,thin=100,verbose=T)


source("R/Function_DecompG.R")
mod.list <- list(modelmulti$VCV,modelmulti_no_ped$VCV)
decomp.all <- decompG(mod.list = mod.list, n.traits = 2, trait.names = c("phen1", "phen2"), sample.size = NULL)
Glist <- decomp.all$Gmatrices

Glist
