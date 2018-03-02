## Functions to calculate the decomposition of the G matrices
library(MCMCglmm)
## Insert list of the VCV of your models
decompG <- function(mod.list, n.traits = 7, trait.names = c("AgeMat", "HMI", "PI", "SizeMat", "SizeNeo", "SpineSizeNeo", "Fec"), sample.size = NULL){

	if(is.null(sample.size)){sample.size <- nrow(mod.list[[1]])}

	## Construct estimates of the G-matrices
	n <- length(mod.list)
	G.list <- array(, dim = c(n.traits, n.traits, n))
	for(i in 1:n){
		G.list[,,i] <- matrix(posterior.mode(mcmc(mod.list[[i]][,1:(n.traits)^2])),n.traits,n.traits)
	}
	dimnames(G.list) = list(trait.names, trait.names, NULL)

	eig.analysis <- apply(G.list,3,eigen)

	## Construct MCMC samples of the G-matrices
	G.MCMC <- array(, dim = c(n.traits, n.traits, n, sample.size))
	for(i in 1:n){
		for(k in 1:sample.size){
			G.MCMC[,,i,k] <- matrix(as.numeric(mod.list[[i]][k,1:(n.traits)^2]),n.traits,n.traits)
		}
	}

	## Project eigenvectors of estimated G through MCMC samples of G to obtain a distribution
	## around your eigenvalues
	eigVal.tmp <- array(, dim = c(sample.size, n.traits, n))
	for(i in 1:n){
		for(k in 1:sample.size){
			for(j in 1:n.traits){
				eigVal.tmp[k,j,i] <- t(eig.analysis[[i]]$vectors[,j])%*%G.MCMC[,,i,k]%*%eig.analysis[[i]]$vectors[,j]
			}
		}
	}

	eigVal.mode <- apply(eigVal.tmp, 2:3, function(x) posterior.mode(mcmc(x)))
	eigVal.HPD <- apply(eigVal.tmp, 2:3, function(x) HPDinterval(mcmc(x)))
	dimnames(eigVal.HPD) = list(c('lower', 'upper'), NULL, NULL)

	## Return the G-matrices, the eigendecomposition of each of the G-matrices,
	## the mode of the eigenvalues (of projecting the eigenvectors of estimated G through MCMC samples,
	## and the HPD interval of the eigenvalues
	return(list(G.MCMCsamples = G.MCMC, Gmatrices = G.list, eigendecomp = eig.analysis, 
		distr.eigval = list(mode = eigVal.mode, HPDinterval = eigVal.HPD)))
}






