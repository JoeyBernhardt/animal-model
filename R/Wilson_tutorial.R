# Wilson et al 2010 tutorial ----------------------------------------------

# read in data from the Wilson tutorial and get it in order  ------------------------------------------------------------
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


