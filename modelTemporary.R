dat <- data.frame(spawners = spawners[2:(n + 1)], year = 1:n)

# simple mean model
mod1 <- lm(spawners ~ 1, data = dat)
summary(mod1)

par(mfrow=c(2,2))
plot(mod1)

par(mfrow=c(1,1))
plot(residuals(mod1), type = "o")
abline(h = 0, lty=2)

acf(residuals(mod1))


# Is there autocorrelation in the real data?
dat.real <- read.csv("Data/PSE_chumSpawnerAbundanceAndTotalRunSize.csv")

i <- 1

dat.i <- subset(dat.real, dat.real$location == unique(dat.real$location)[i] & dat.real$parameter == "LGL counts")

plot(dat.i$year, dat.i$datavalue)

plot(dat.i$year, log(dat.i$datavalue), type="o", ylim = c(0, 14))
quantile(dat.i$datavalue, c(0.25, 0.75))

# simple mean model
mod1 <- lm(datavalue ~ 1, data = dat.i)
summary(mod1)

par(mfrow=c(1,1))
plot(residuals(mod1), type = "o")
abline(h = 0, lty=2)

acf(residuals(mod1))

autocor <- matrix(NA, nrow = 15, ncol = 6)
for(i in 1:15){
	dat.i <- subset(dat.real, dat.real$location == unique(dat.real$location)[i] & dat.real$parameter == "LGL counts")
	mod1 <- lm(datavalue ~ 1, data = dat.i)
	a <- acf(residuals(mod1))
	autocor[i, ] <- a[[1]][1:6]
}


plot(1:6, autocor[1,], "o", ylim=range(autocor, na.rm=TRUE))
for(i in 2:15){
	points(1:6, autocor[i, ], "o")
}
 abline(h = 0.2, lty=2, col=4)

 t.test(autocor[,2]) 
 
 
 
 ##############################################################################
 # Simulated data
 ##############################################################################
 
# see https://eranraviv.com/bootstrapping-time-series-r-code/
y <- spawners[2:(n+1)]
 
numLags <- 1
ar1 <- ar.mle (y, demean = TRUE, intercept = FALSE, order.max = numLags, aic = FALSE) # standard OLS

cbind(y - ar1$x.mean, ar1$resid)

# Sample residuals
res.star <- matrix(nrow = n, ncol = nBoot)
res.star[, ] <- sample(na.omit(scaled.res), n * nBoot, replace = TRUE)

obs.star <- matrix(nrow = n + numLags, ncol = nBoot) # will hold the bootstrapped series

for (i in 1:nBoot){
	# Initialize the simualted time series using the true data for the first
	# 1:numLags points, starting from a random place in the timeseries
	j.init <- sample(1:n, 1) # starting point for intialization
	while((j.init + numLags - 1) > n) j.init <- sample(1:n, 1)
	obs.star[1:numLags, i] <- y[j.init:(j.init + numLags - 1)] 
	
	for (j in 1:n){ # For each timepoint in the simulated series
		obs.star[(numLags + j), i] <- ar1$x.mean + ar1$ar %*% (obs.star[j:(j + numLags - 1), i] - ar1$x.mean) + res.star[j, i]
		} #end j

} # end boot

	# Sanity check:
plot(1:n, y, type="o", col=grey(0.8))
lines(1:n, obs.star[(numLags + 1:n), i], type = "o") 
abline(h = ar1$x.mean)

#------------------------------------------------------------
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


HS_benchCI_model <- modelBoot(series = spawners[2:(n+1)], numLags = 1, nBoot = 10000)


# plot simulated timeseries
par(mfrow = c(1,1), mar=c(4,4,2,1), oma=rep(0,4), mgp=c(3,1,0))
plot(1:n, spawners[2:(n+1)]*10^-3, "n", xlab="Time (years)", ylab = "Spawners (thousands)", las=1, col=2)
for(i in sample(1:nBoot, 10)){
	lines(1:n, tail(HS_benchCI_model[[2]][, i], n)*10^-3, col="#00000030")
	abline(h = 0)
}
lines(1:n, spawners[2:(n+1)]*10^-3, col=2)
