#______________________________________________________________________________
#' Function to draw n realized harvest rates, with mean targetHarvest 
#' 
#' This function generates harvest rates bound between zero and one,
#' incorporating beta error around a target harvest rate, h'.
#'
#' @param targetHarvest A numeric value or vector giving the target harvest 
#' rate(s)
#' @param sigmaHarvest A single numeric value for the standard deviation in 
#' error around the target harvest rate
#' @param n The number of realized harvest rates to return. If 
#' targetHarvest is a vector, then n must equal length(targetHarvest). May be 
#' the number of years or number of populations.
#' 
#' @return Returns a numeric vector of realized harvest rates for each n
#'
#' @examples
#'
#' @export
realizedHarvest <- function(n, targetHarvest, sigmaHarvest){
	
	beta1 <- (targetHarvest^2 - targetHarvest^3 - sigmaHarvest^2*targetHarvest)/(sigmaHarvest^2)
	beta2 <- (targetHarvest * (1 - targetHarvest)^2 - sigmaHarvest^2*(1 - targetHarvest))/(sigmaHarvest^2)
	
	harvestRate <- rbeta(n = n, shape1 = beta1, shape2 = beta2)
	
	return(harvestRate)
}

#______________________________________________________________________________
#' Function to simulate spawner data
#' 
#' This function simualtes a spawner-recruit dataset using the Ricker model, 
#' incorporating stochastic realized harvest rates (see function 
#' \code{realizedHarvest}) and temporal autocorrelation in the recruitment deviates.
#' The time series is seeded with a a spawner abundance of 20% Seq, or
#' 0.2 * (a / b). Requires function \code{realizedHarvest}.
#'
#' @param n Number of years to simulate data
#' @param a A numeric value for alpha parameter in Ricker productivity at low spawner
#' abundance.
#' @param b A numeric vector of beta values, i.e. density dependence para-
#' meter.
#' @param sigma A numeric value for the standard deviation in recruitment deviates
#' (log scale) each year.
#' @param tau A numeric value for the temporal autocorrelation in recruitment 
#' deviates
#' @param targetHarvest A numeric value or vector giving the target harvest 
#' rate(s)
#' @param sigmaHarvest A single numeric value for the standard deviation in 
#' error around the target harvest rate
#' @param seed Option to set the seed for simulations to make output reproducible
#' (defaults to NULL)
#' 
#' @return Returns a matrix with columns for spawners, recruits, the annual 
#' recruitment deviates (phi), and the annual realized harvest rate.
#'
#' @examples
#'
#' @export

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
		recruits[i + 1] <- spawners[i] * exp(a - b * spawners[i]) * exp(phi[i + 1])
		spawners[i + 1] <- (1 - harvestRate[i]) * recruits[i + 1]
	}
	
	out <- cbind(spawners = spawners, recruits = recruits, phi = phi, harvestRate = c(0, harvestRate))
	return(out)
}


#______________________________________________________________________________
#' Function to bootstrap confidence intervals on the historical spawners
#' benchmarks
#' 
#' This function takes a time series of spawner abundances and bootstraps
#' 95% confidence intervals on the 25th and 50th percentiles of historical 
#' spawner abundance. Bootstrapping can either be naiive (blockLength = 1), 
#' which means that temporal autocorerlation in the series is ignored when
#' calculating confidence invervals, or block boostrapping can be used when 
#' \code{blockLength > 1}. For information on boostrapping, including the block
#' boostrap, see 
#' http://anson.ucdavis.edu/~peterh/sta251/bootstrap-lectures-to-may-16.pdf
#' 
#' 
#' @param series Time series of spawner abundance data. Missing values should be
#' entered as NA
#' @param blockLength The number of years for each block when implementing block
#' bootstrapping
#' @param nBoot The number of permutations of the timeseries to be used for 
#' calculating the bootstrap confidence intervals on benchmarks
#' @param benchmarks The quantiles of the historical spawners series to be used
#' as the upper and lower benchmarks. Defaults to the 25th and 50th percentiles (
#' code{benchmarks = c(0.25, 0.5)}).
#' 
#' @return Returns a list; the first element is a matrix with the 95% CI 
#' (rows) on the lower and upper benchmarks (columns). The second element
#' is a matrix of dimension number of timesteps (\code{length(series)}) by 
#' \code{nBoot}containing the permuted timeseries (columns).
#'
#' @examples
#'
#' @export

blockBoot <- function(
	series, 
	blockLength = 1, 
	nBoot = 10000, 
	benchmarks = c(0.25, 0.5)
	){
	
	n <- length(series)
	
	# Number of blocks to be created
	nBlocks <- ceiling(n/blockLength) 
	# Note: If the length of series (n) does not divide evenly by blockLength,
	# select the first 1:n elements of the new series
	
	# Matrix to store bootstrapped series
	obs.star <- matrix(nrow = n, ncol = nBoot)
	
	# Matrix to store bootstrapped values of CI
	HS_benchBoot <- matrix(
		NA, 
		nrow = nBoot, 
		ncol = 2, 
		dimnames = list(c(1:nBoot), c("lower", "upper")))
	
	# Block bootstrap:
	for(i in 1:nBoot){
		
		j <- sample(1:(n - blockLength + 1), 
								size = nBlocks, replace = TRUE)
		newIndex <- rep(j, each = blockLength) + rep(0:(blockLength - 1), nBlocks)
		
		obs.star[, i] <- series[newIndex[1:n]]
		
		HS_benchBoot[i, ] <- quantile(obs.star[, i], benchmarks)
	} # end bootstrap loop
	
	HS_benchCI <- apply(HS_benchBoot, 2, quantile, c(0.025, 0.975))
	
	return(list(CI = HS_benchCI, permutedSeries = obs.star))
	
}

#______________________________________________________________________________
#' Function to boostrap confidence intervals based on modelled timeseries of
#' residuals 
#' 
#' **Description needed.**
#' 
#' 
#' @param series Time series of spawner abundance data. Missing values should be
#' entered as NA
#' @param numLags The number of years for each block when implementing block
#' bootstrapping
#' @param nBoot The number of permutations of the timeseries to be used for 
#' calculating the bootstrap confidence intervals on benchmarks
#' @param benchmarks The quantiles of the historical spawners series to be used
#' as the upper and lower benchmarks. Defaults to the 25th and 50th percentiles (
#' code{benchmarks = c(0.25, 0.5)}).
#' 
#' @return Returns a list; the first element is a matrix with the 95% CI 
#' (rows) on the lower and upper benchmarks (columns). The second element
#' is a matrix of dimension number of timesteps (\code{length(series)}) by 
#' \code{nBoot}containing the simulated timeseries (columns).
#'
#' @examples
#'
#' @export


modelBoot <- function(
	series, 
	numLags = 1, # numLags is the lag for the autocorrelation; default is just 1 year
	nBoot = 10000, 
	benchmarks = c(0.25, 0.5)
	){
	
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
		dimnames = list(c(1:nBoot), c("lower", "upper")))
	
	# Model bootstrap:
	for(i in 1:nBoot){
		
		# Initialize the simualted time series using the true data for the first
		# 1:numLags points, starting from a random place in the timeseries
		j.init <- sample(1 : (n - numLags + 1), 1) # starting point for intialization
		obs.star.log[1:numLags, i] <- log(series[j.init:(j.init + numLags - 1)])
		
		for (j in 1:n){ # For each timepoint in the simulated series
			obs.star.log[(numLags + j), i] <- ar.fit$x.mean + ar.fit$ar %*% (obs.star.log[j:(j + numLags - 1), i] - ar.fit$x.mean) + res.star[j, i]
		} #end j
		
		HS_benchBoot[i, ] <- quantile(exp(obs.star.log[(numLags + 1):(numLags + n), i]), benchmarks)
	} # end bootstrap loop
	
	HS_benchCI <- apply(HS_benchBoot, 2, quantile, c(0.025, 0.975))
	
	obs.star <- exp(tail(obs.star.log, n))
	
	return(list(CI = HS_benchCI, simulatedSeries = obs.star))
	
}