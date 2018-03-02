## Functions to calculate the 4th order tensor
library(MCMCglmm)
library(matrixcalc)

## Function to calculate the tensor as in Robinson & Beckerman (2013)
covten <- function(matrices, n.traits = 7, n.G = 3){ #for n.G matrices
	n <- n.G
	covcov <- matrix(0, nrow = n.traits*(n.traits+1)/2, ncol = n.traits*(n.traits+1)/2)
	
	## upper left quadrant of S
	for(i in 1:n.traits){
		for(j in 1:n.traits){
		    covcov[i,j]=cov(matrices[i,i,1:n],matrices[j,j,1:n])
		}
	}

	## lower left and upper right quadrants of S
	for(k in 1:n.traits){
		count <- 1
	  	for(i in 1:(n.traits-1)){
    			for(j in (i+1):n.traits){
		      	covcov[k,count+n.traits]=sqrt(2)*cov(matrices[k,k,1:n],matrices[i,j,1:n])
      			covcov[count+n.traits,k]=sqrt(2)*cov(matrices[k,k,1:n],matrices[i,j,1:n])
				count=count+1
			}
		}
	}

	## lower right quadrant of S
	countx <- 1
	county <- 1
	for(k in 1:(n.traits-1)){
		for(l in (k+1):n.traits){
			for(i in 1:(n.traits-1)){
				for(j in (i+1):n.traits){
					covcov[countx+n.traits,county+n.traits]=2*cov(matrices[k,l,1:n],matrices[i,j,1:n])
					covcov[county+n.traits,countx+n.traits]=2*cov(matrices[k,l,1:n],matrices[i,j,1:n])
					county=county+1}
			}
		countx=countx+1
		county=1
		}
	}
return(covcov)
}

## Construct eigentensors for each fourth order tensor S
## Eigentensor are arranged in one big matrix
## So first eigentensor is 
eigtensor <- function(S, n.traits = 7, m = 3){ #m is for the amount of populations you are using
	eigvecs <- eigen(S)$vectors
	eigvals <- eigen(S)$values
	eigten <- matrix(nrow = n.traits*n.traits*(n.traits+1)/2, ncol = n.traits)

	for(eval in 1:(n.traits*(n.traits+1)/2)){
		for(i in 1:n.traits){ #The diagonal
			eigten[i+n.traits*(eval-1),i] <- eigvecs[i,eval]}
	count <- 1
		for(i in 1:(n.traits-1)){
			for(j in (i+1):n.traits){
				eigten[i+n.traits*(eval-1),j] <- 1/sqrt(2)*eigvecs[count+n.traits,eval]
				eigten[j+n.traits*(eval-1),i] <- 1/sqrt(2)*eigvecs[count+n.traits,eval]
				count <- count+1
			}
		}
	}
return(list(eigvalues = eigvals, eigvectors = eigvecs, eigtensor = eigten))
}

## Function to calculate the 4th order tensor and everything else
## Using the S estimated from the G-matrices themselves
## Note! Another option would be to estimate S from the MCMC samples, 
## however, eigendecomposition do not result here in the expected amount of eigenvalues = 0!!
tensor <- function(mod.list, n.traits = 7, trait.names = c("AgeMat", "HMI", "PI", "SizeMat", "SizeNeo", "SpineSizeNeo", "Fec"), sample.size = NULL){

	if(is.null(sample.size)){sample.size <- nrow(mod.list[[1]])}

	## Construct estimates of the G-matrices
	n <- length(mod.list)
	G.list <- array(, dim = c(n.traits, n.traits, n))
	for(i in 1:n){
		G.list[,,i] <- matrix(posterior.mode(mcmc(mod.list[[i]][,1:(n.traits)^2])),n.traits,n.traits)
	}
	dimnames(G.list) = list(trait.names, trait.names, NULL)

	## Calculate covariance tensor on the based G matrices
	## And perform its eigendecomposition
	covten.obs <- covten(G.list, n.traits = n.traits, n.G = n)
	eigdecompten <- eigen(covten.obs)

	## Construct MCMC samples of the S matrix
	G.MCMC <- array(, dim = c(n.traits, n.traits, n, sample.size))
	for(i in 1:n){
		for(k in 1:sample.size){
			G.MCMC[,,i,k] <- matrix(as.numeric(mod.list[[i]][k,1:(n.traits)^2]),n.traits,n.traits)
		}
	}
	S.MCMC <- array(, dim = c(n.traits*(n.traits+1)/2, n.traits*(n.traits+1)/2, sample.size))
	for(k in 1:sample.size){
		S.MCMC[,,k] <- covten(G.MCMC[,,,k], n.traits = n.traits, n.G = n)
	}

	## Estimate matrix S from the MCMC samples
	S.mean <- apply(S.MCMC, 1:2, function(x) posterior.mode(mcmc(x)))
	eigdecomp.Smean <- eigen(S.mean)

	## Project eigenvectors of estimated S through MCMC samples of S to obtain a distribution
	## around its eigenvalues
	eigVal.tmp <- array(, dim = c(sample.size, n.traits*(n.traits+1)/2))
	for(k in 1:sample.size){
		for(j in 1:(n.traits*(n.traits+1)/2)){
			eigVal.tmp[k,j] <- t(eigdecompten$vectors[,j])%*%S.MCMC[,,k]%*%eigdecompten$vectors[,j]
		}
	}

	eigVal.mode <- apply(eigVal.tmp, 2, function(x) posterior.mode(mcmc(x)))
	eigVal.HPD <- apply(eigVal.tmp, 2, function(x) HPDinterval(mcmc(x)))
	dimnames(eigVal.HPD) = list(c('lower', 'upper'), NULL)

	## Construct the eigentensors of the observed S and of the MCMC samples
	eigten.obs <- eigtensor(covten.obs)$eigtensor
	eigten.val <- array(, dim = c(n.traits*(n.traits+1)/2, n.traits))
	eigten.vec <- array(, dim = c(n.traits, n.traits,n.traits*(n.traits+1)/2))
	for(i in 1:(n.traits*(n.traits+1)/2)){
		eigten.val[i,] <- eigen(eigten.obs[(n.traits*(i-1)+1):(i*n.traits),])$values
		eigten.vec[,,i] <- eigen(eigten.obs[(n.traits*(i-1)+1):(i*n.traits),])$vectors
	}

	eigtensor.MCMC <- array(, dim = c(n.traits*n.traits*(n.traits+1)/2, n.traits, sample.size))
	for(k in 1:sample.size){
		eigtensor.MCMC[,,k] <- eigtensor(S.MCMC[,,k])$eigtensor
	}

	## Find the C^i,j coordinates of the observed ones and of the samples
	## See how each matrix relates to the space described by a specific eigentensor
	Ccoord.obs <- array(, dim = c(n.traits*(n.traits+1)/2, n))
	for(i in 1:(n.traits*(n.traits+1)/2)){
		for(j in 1:n){
			Ccoord.obs[i,j] <- frobenius.prod(eigten.obs[(n.traits*(i-1)+1):(i*n.traits),], G.list[,,j])
		}
	}

	Ccoord.MCMC <- array(, dim = c(n.traits*(n.traits+1)/2, n, sample.size))
	for(k in 1:sample.size){
		for(i in 1:(n.traits*(n.traits+1)/2)){
			for(j in 1:n){
				Ccoord.MCMC[i,j,k] <- frobenius.prod(eigten.obs[(n.traits*(i-1)+1):(i*n.traits),], G.MCMC[,,j,k])
			}
		}
	}

	Ccoord.mode <- apply(Ccoord.MCMC, 1:2, function(x) posterior.mode(mcmc(x))) 
	Ccoord.HPD <- apply(Ccoord.MCMC, 1:2, function(x) HPDinterval(mcmc(x))) 
	dimnames(Ccoord.HPD) <- list(c('lower', 'upper'), NULL, NULL)

	## Return the observed S-matrix, the eigendecomposition of each of the S-matrix, the samples of S-matrix, 
	## the mode of the eigenvalues (of projecting the eigenvectors of estimated S through MCMC samples,
	## and the HPD interval of the eigenvalues, the observed eigentensors (calculated from the observed S),
	## the eigendecomposition of these observed eigentensors, the MCMC samples of the eigentensors of all MCMC samples of S,
	## the C coordinates: observed values, mode and HPD interval, and the MCMC samples of the C coordinates
	return(list(Smatrix = covten.obs, eigendecompS = eigdecompten, S.samples = S.MCMC,
		distr.eigvalS = list(mode = eigVal.mode, HPDinterval = eigVal.HPD),
		eigten.obs = eigten.obs, eigtendecomp = list(values = eigten.val, vectors = eigten.vec),
		eigten.samples = eigtensor.MCMC, Ccoord = list(obs = Ccoord.obs, mode = Ccoord.mode, HPD = Ccoord.HPD),
		Ccoord.MCMC = Ccoord.MCMC))
}

plotCcoord <- function(Ccoord.val, Ccoord.HPD, sign.ten = 2, ymin = NULL, ymax = NULL){
	n <- ncol(Ccoord.val)

	if(is.null(ymin)){
		ymin <- min(c(Ccoord.HPD[,1:sign.ten,]))
		ymax <- max(c(Ccoord.HPD[,1:sign.ten,]))
	} else { ymin <- ymin
		ymax <- ymax
	}

	symbol <- rep(21:25, 3)

	plot(0:(n+1), rep(-100,n+2), ylim = c(ymin,ymax))
	for(j in 1:sign.ten){
		for(i in 1:n){
			points(i+j/10, Ccoord.val[j,i], pch = symbol[j], cex = 1.5)
			#points(i+j/10, Ccoord.HPD[1,j,i], pch = '-', cex = 1.5)
			#points(i+j/10, Ccoord.HPD[2,j,i], pch = '-', cex = 1.5)
			segments(i+j/10, Ccoord.HPD[1,j,i], i+j/10, Ccoord.HPD[2,j,i])
		}
		for(i in 1:(n-1)){
			segments(i+j/10, Ccoord.val[j,i], i+j/10+1, Ccoord.val[j,i+1])
		}
	}
}

plotCcoord_withoutHPD <- function(Ccoord.val, sign.ten = 2, ymin = NULL, ymax = NULL){
	n <- ncol(Ccoord.val)

	if(is.null(ymin)){
		ymin <- min(c(Ccoord.val[1:sign.ten,]))
		ymax <- max(c(Ccoord.val[1:sign.ten,]))
	} else { ymin <- ymin
		ymax <- ymax
	}

	symbol <- rep(21:25, 3)

	plot(0:(n+1), rep(-100,n+2), ylim = c(ymin,ymax))
	for(j in 1:sign.ten){
		for(i in 1:n){
			points(i+j/10, Ccoord.val[j,i], pch = symbol[j], cex = 1.5)
		}
		for(i in 1:(n-1)){
			segments(i+j/10, Ccoord.val[j,i], i+j/10+1, Ccoord.val[j,i+1])
		}
	}
}

