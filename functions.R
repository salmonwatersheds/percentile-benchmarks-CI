###############################################################################

###############################################################################
# Function to draw n realized harvest rates, with mean targetHarvest:
realizedHarvest <- function(n, targetHarvest, sigmaHarvest){
	
	beta1 <- (targetHarvest^2 - targetHarvest^3 - sigmaHarvest^2*targetHarvest)/(sigmaHarvest^2)
	beta2 <- (targetHarvest * (1 - targetHarvest)^2 - sigmaHarvest^2*(1 - targetHarvest))/(sigmaHarvest^2)
	
	harvestRate <- rbeta(n = n, shape1 = beta1, shape2 = beta2)
	
	return(harvestRate)
}

###############################################################################

###############################################################################
# Function to simulate spawner data:
simSpawners <- function(
	n = 50, 
	a = 1.4, 
	b = 0.0001, 
	sigma = 0.1, 
	tau = 0.8, 
	targetHarvest = 0.42, 
	sigmaHarvest = 0.13, 
	seed = NULL) {
	
	spawners <- numeric(n + 1)			# vector to store spawner abundance
	spawners[1] <- 0.2 * (a / b)		# initiate spawner abundance at 20% Seq
	
	phi <- numeric(n + 1)						# vector to store recruitment deviates
	phi[1] <- rnorm(1, 0, sigma)		# initiate recruitment deviates assuming phi[0] = 0
	recruits <- numeric(n + 1)			# vector to store recruits
	
	if(is.null(seed) == FALSE) set.seed(seed)
	
	harvestRate <- realizedHarvest(n, targetHarvest, sigmaHarvest)
	
	for(i in 1:n){
		phi[i + 1] <- tau * phi[i] + rnorm(1, 0, sigma) 
		recruits[i + 1] <-spawners[i] * exp(a - b * spawners[i]) * exp(phi[i + 1])
		spawners[i + 1] <- (1 - harvestRate[i]) * recruits[i + 1]
	}
	
	return(cbind(spawners, recruits, phi, c(0, harvestRate)))
}


###############################################################################

###############################################################################

blockBoot <- function(series, blockLength = 10, nBoot = 10000){
	
	# Number of blocks to be created
	nBlocks <- ceiling(length(series)/blockLength) 
	# Note: If the length of series (n) does not divide evenly by blockLength,
	# select the first 1:n elements of the new series
	
	# Matrix to store bootstrapped series
	obs.star <- matrix(nrow = n, ncol = nBoot)
	
	# Matrix to store bootstrapped values of CI
	HS_benchBoot <- matrix(
		NA, 
		nrow = nBoot, 
		ncol = 2, 
		dimnames = list(c(1:nBoot), c("S25", "S50")))
	
	# Block bootstrap:
	for(i in 1:nBoot){
		
		j <- sample(1:(length(series) - blockLength + 1), 
								size = nBlocks, replace = TRUE)
		newIndex <- rep(j, each = blockLength) + rep(0:(blockLength - 1), nBlocks)
		
		obs.star[, i] <- series[newIndex[1:length(series)]]
		
		HS_benchBoot[i, ] <- quantile(obs.star[, i], c(0.25, 0.5))
	} # end bootstrap loop
	
	HS_benchCI <- apply(HS_benchBoot, 2, quantile, c(0.025, 0.975))
	
	return(list(HS_benchCI, obs.star))
	
}

###############################################################################

###############################################################################
modelBoot <- function(series, numLags = 1, nBoot = 10000){
	# numLags is the lag for the autocorrelation; default is just 1 year
	
	n <- length(series)
	# Fit model to estimate autocorrelation
	ar.fit <- ar.mle(
		log(series), # spawner time series
		demean = TRUE, # Estimate mean spawners
		intercept = FALSE, # Intercept = 0 for autocorrelation
		order.max = numLags, # lag for autocorrelation
		aic = FALSE # estimate autocorrelation for all numLags 
	) # standard OLS
	
	# Matrices to store bootstrapped residuals and spawners (obs)
	res.star <- matrix(nrow = n, ncol = nBoot)
	res.star[, ] <- sample(na.omit(ar.fit$resid), n * nBoot, replace = TRUE)
	
	obs.star.log <- matrix(nrow = n + numLags, ncol = nBoot) 
	
	# Matrix to store bootstrapped values of CI
	HS_benchBoot <- matrix(
		NA, 
		nrow = nBoot, 
		ncol = 2, 
		dimnames = list(c(1:nBoot), c("S25", "S50")))
	
	# Model bootstrap:
	for(i in 1:nBoot){
		
		# Initialize the simualted time series using the true data for the first
		# 1:numLags points, starting from a random place in the timeseries
		j.init <- sample(1:n, 1) # starting point for intialization
		while((j.init + numLags - 1) > n) j.init <- sample(1:n, 1)
		obs.star.log[1:numLags, i] <- log(y[j.init:(j.init + numLags - 1)])
		
		for (j in 1:n){ # For each timepoint in the simulated series
			obs.star.log[(numLags + j), i] <- ar.fit$x.mean + ar.fit$ar %*% (obs.star.log[j:(j + numLags - 1), i] - ar.fit$x.mean) + res.star[j, i]
		} #end j
		
		HS_benchBoot[i, ] <- quantile(exp(obs.star.log[(numLags + 1):(numLags + n), i]), c(0.25, 0.5))
	} # end bootstrap loop
	
	HS_benchCI <- apply(HS_benchBoot, 2, quantile, c(0.025, 0.975))
	
	obs.star <- exp(tail(obs.star.log, n))
	
	return(list(HS_benchCI, obs.star))
	
}