
### Goal here is to learn about the animal model, and figure out how you can run it without a pedigree.
### Use sample data from Pierre de Villemereuil’s tutorial with and without a pedigree. Compare results.


# Load packages -----------------------------------------------------------

library(MCMCglmm)



# working with Pierre de Villemereuil’s tutorial (http://devillemereuil.legtux.org/wp-content/uploads/2012/12/tuto_en.pdf) --------------------------

## read in data
pedigreemulti <- read.table('data-raw/pedigreemulti.txt',header=T)
datamulti <- read.table('data-raw/datamulti.txt',header=T)
head(datamulti)

## read in processed data

models_run <- load("data-processed/modelmulti.RData") 

## specify priors for the animal effect and the residual variance
prior <- list(R=list(V=diag(2)/2,nu=2),
		 G=list(G1=list(V=diag(2)/2,nu=2)))


modelmulti <- MCMCglmm(cbind(phen1,phen2)~trait-1, ## the notation trait-1 allows us to estimate the mean of each trait, instead of the contrast between the two traits
					 random=~us(trait):animal, ##  define the structure of the variance-covariance matrix for the random effects 
					 rcov=~us(trait):units,  ##  define the structure of the variance-covariance matrix for the residual variances
					 family=c("gaussian","gaussian"), ## family requires a vector with the data distribution of each trait
					 prior=prior,
					 pedigree=pedigreemulti, ## specifies the pedigree; one row for each individual, each column specifies the parents for each individual
					 data=datamulti, ## includes two traits for each individual, i.e. phen1 and phen 2
					 nitt=100000,
					 burnin=1000,
					 thin=10)

summary(modelmulti)

head(modelmulti$VCV)

## calculate heritabilities
herit1<-modelmulti$VCV[,'traitphen1:traitphen1.animal']/ (modelmulti$VCV[,'traitphen1:traitphen1.animal']+ modelmulti$VCV[,'traitphen1:traitphen1.units'])
herit2<-modelmulti$VCV[,'traitphen2:traitphen2.animal']/(modelmulti$VCV[,'traitphen2:traitphen2.animal']+modelmulti$VCV[,'traitphen2:traitphen2.units'])
mean(herit1)
mean(herit2)

## calculate genetic correlations among the two traits
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
					   nitt=100000,
					   burnin=1000,
					   thin=10)


summary(modelmulti)
summary(modelmulti_no_ped)

## calculate heritabilities
herit1_no_ped<-modelmulti_no_ped$VCV[,'traitphen1:traitphen1.animal']/ (modelmulti_no_ped$VCV[,'traitphen1:traitphen1.animal']+ modelmulti_no_ped$VCV[,'traitphen1:traitphen1.units'])
herit2_no_ped<-modelmulti_no_ped$VCV[,'traitphen2:traitphen2.animal']/(modelmulti_no_ped$VCV[,'traitphen2:traitphen2.animal']+modelmulti_no_ped$VCV[,'traitphen2:traitphen2.units'])
mean(herit1_no_ped)
mean(herit2_no_ped)

## calculate genetic correlations among the two traits
corr.gen_no_ped<-modelmulti_no_ped$VCV[,'traitphen1:traitphen2.animal']/
	sqrt(modelmulti_no_ped$VCV[,'traitphen1:traitphen1.animal']*modelmulti_no_ped$VCV[,'traitphen2:traitphen2.animal'])
mean(corr.gen_no_ped)


save(modelmulti, modelmulti_no_ped, file = "data-processed/modelmulti.rdata")


## now run the same thing, but give it a fake, clonal pedigree, with all the parents having NA's, or all the parents being different

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




