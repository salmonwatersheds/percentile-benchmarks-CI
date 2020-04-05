source("functions.R")
###############################################################################
# Calculate "true" confidence intervals
###############################################################################

#------------------------------------------------------------------------------
# Parameters

a <- 1.4							# productivity
b <- 0.0001 					# density dependence
sigma <- 0.1  				# variance in recruitment deviates
tau <- c(0.1, 0.3, 0.6, 0.9) # temporal autocorrelation in recruitment deviates
n <- 50								# number of years of data to simulate 
benchmarks <- c(0.25, 0.5) # HS benchmarks
nBoot <- 10000				# Number of iterations in bootstrap and trueCI calculation

#------------------------------------------------------------------------------
# Set up matrices to store time series and true CI

spawnersTrue <- matrix(NA, nrow = n, ncol = length(tau))
HS_bench <- matrix(NA, nrow = 2, ncol = length(tau), dimnames = list(c("lower", "upper"), NULL))
trueCI <- array(rep(NA, 2*2*length(tau)), dim = c(2, 2, length(tau)), dimnames = list(c("2.5%", "97.5%"), c("lower", "upper"), NULL))
	
for(k in 1:length(tau)){
	
	#------------------------------------------------------------------------------
	# Simulate spawner recruit data
	simDat <- simSpawners(seed = 3958, tau = tau[k])
	spawnersTrue[, k] <- simDat[2:(n + 1), 'spawners']
	
	HS_bench[, k] <- quantile(spawnersTrue[, k], benchmarks)
	
	# True CI: Simulate time series using the Ricker parameters
	trueBench <- matrix(nrow = nBoot, ncol = 2)
	for(i in 1:nBoot){
		simDat.i <- simSpawners(tau = tau[k])
		trueBench[i, ] <- quantile(simDat.i[2:(n + 1), 'spawners'], benchmarks)
	}
	trueCI[, , k] <- apply(trueBench, 2, quantile, c(0.025, 0.975))
}

###############################################################################
# Calculate confidence intervals using different approaches
###############################################################################

CI.all <- list(); length(CI.all) <- length(tau)

for(k in 1:length(tau)){
	CI.all[[k]] <- list(
		naiveCI = blockBoot(series = spawnersTrue[, k], blockLength = 1, nBoot = nBoot),
		block5CI = blockBoot(series = spawnersTrue[, k], blockLength = 5, nBoot = nBoot),
		block10CI = blockBoot(series = spawnersTrue[, k], blockLength = 10, nBoot = nBoot),
		block15CI = blockBoot(series = spawnersTrue[, k], blockLength = 15, nBoot = nBoot),
		modelCI = modelBoot(series = spawnersTrue[, k], numLags = 1, nBoot = nBoot)
	)
}
	

###############################################################################
# Plot examples of the bootstrapped series
###############################################################################

quartz(width = 7.5, height = 4.5, pointsize = 10)
m <- 8
M <- matrix(rep(1:15, each = m), nrow = 3, byrow = TRUE)
layout(cbind(M[, 1:m], rep(max(M)+1, 3), M[, (m+1):(m*4)], rep(max(M) + 2, 3), M[, (m*4 +1):(m*5)]))

par(mar = c(1,1,0,0), oma=c(2,2,4,2))
for(i in 1:3){
	for(j in 1:5){
		plot(1:n, spawners*10^-3, "n", xaxt="n", las=1, ylab="", xlab="", yaxt="n")
		lines(1:n, spawners*10^-3, col=grey(0.8), lwd=1.2, lty= 3)#, "o", pch=21, bg="white")
		abline(h = c(S25, S50)* 10^-3, col= statusCols[c('r', 'g')], lwd=1.2, lty = 3)
		# abline(h = mean(spawners*10^-3), col=grey(0.8))
		lines(1:n, CI.all[[j]][[2]][, i]*10^-3)#, "o", pch=19)
		abline(h = quantile(CI.all[[j]][[2]][, i], c(0.25, 0.5))* 10^-3, col= statusCols[c('r', 'g')])
		
		if(j == 1 & i == 1) mtext(side = 3, "Naive", line=2.5)
		if(j == 3 & i == 1) mtext(side = 3, "Block", line=2.5)
		if(j == 5 & i == 1) mtext(side = 3, "Model", line=2.5)
		
		if(j == 1 & i == 1) mtext(side = 3, "(b = 1)", line=1)
		if(j == 2 & i == 1) mtext(side = 3, "b = 5", line=1)
		if(j == 3 & i == 1) mtext(side = 3, "b = 10", line=1)
		if(j == 4 & i == 1) mtext(side = 3, "b = 15", line=1)
	}
	mtext(side = 4, line = 1, paste("Example ", i))
}	
mtext(side = 1, outer = TRUE, "Year")
mtext(side = 2, outer=TRUE, "Spawners")

###############################################################################
# Compare confidence intervals
###############################################################################

# Compile output
CI.compiled <- list(
	S25 = matrix(nrow = 6, ncol = 2),
	S50 = matrix(nrow = 6, ncol = 2))
CI.compiled[['S25']][1, ] <- trueCI[, 1]
CI.compiled[['S50']][1, ] <- trueCI[, 2]
for(j in 1:5){
	CI.compiled[['S25']][j + 1, ] <- CI.all[[j]][[1]][, 1]
	CI.compiled[['S50']][j + 1, ] <- CI.all[[j]][[1]][, 2]
}

quartz(width = 4, height = 3, pointsize = 9)
par(mfrow = c(1,1), mar=c(4,4,2,1), oma = c(0,0,0,0))

plot(c(0, 1, 2:4 + 0.5, 6), rep(c(S25, S50)[i], 6), "n", xlim=c(-0.5, 6.5), xaxt="n", xlab="", ylab="Benchmark", las=1, ylim=range(CI.compiled))
for(i in 1:2){ # for S25 and S50
	plotCI(c(0, 1, 2:4 + 0.5, 6)+c(-0.05, 0.05)[i], rep(c(S25, S50)[i], 6), pch = NA, li = CI.compiled[[i]][, 1], ui = CI.compiled[[i]][, 2], gap = 0, xlim=c(-0.5, 6.5), add=TRUE, sfrac = 0.005, lwd=1.2, col = statusCols[c('r', 'g')][i])
	abline(h = c(S25, S50)[i], col = statusCols[c('r', 'g')][i])
}

axis(side = 1, at = c(0, 1, 2.5, 3.5, 4.5, 6), labels = c("True", "Naive", "b=5", "b=10", "b=15", "Model"), cex.axis=1)
u <- par('usr')
l1 <- u[3] - (u[4] - u[3])/6
l2 <- u[3] - (u[4] - u[3])/4.5
segments(x0 = 2, x1 = 5, y0 = l1, y1=l1, xpd=NA)
text(3.5, l2, "Block", xpd=NA)


###############################################################################
# On real data
###############################################################################
