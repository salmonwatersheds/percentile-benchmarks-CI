###############################################################################
# This code compares three approaches to computing confidence intervals on
# HS benchmarks:
# 1) Naive bootstrap (i.e., block bootstrap with block length of b = 1)
# 2a) Block bootstrap with b = 5
# 2a) Block bootstrap with b = 10
# 2a) Block bootstrap with b = 15
# 3) Model-based approach (bootstrapping residuals)
# Trying three different lengths of blocks for block bootstrapping, we compared
# five different methods in total. See README for further details.
###############################################################################

library(gplots)
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

for(k in 1:length(tau)){ # for each value of autocorrelation
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

k <- 3 # Which level of autocorrelation to use? Default tau = 0.6

# Layout for plot
pdf(file = "Figures/BootstrapSeries.pdf", width = 8, height = 4.5, pointsize = 10)
m <- 8
M <- matrix(rep(1:15, each = m), nrow = 3, byrow = TRUE)
layout(t(cbind(M[, 1:m], rep(max(M)+1, 3), M[, (m+1):(m*4)], rep(max(M) + 2, 3), M[, (m*4 +1):(m*5)])))
par(mar = c(1,1,0,0), oma=c(2,2,4,6))

# Plot
for(i in 1:3){ # for three different bootstrap iterations (columns)
	for(j in 1:5){ # for each CI approach
		plot(1:n, spawnersTrue[, k]*10^-3, "l", xaxt="n", las=1, ylab="", xlab="", yaxt="n", col=grey(0.8), lwd=1.2, lty= 1)
		
		# True benchmarks
		abline(h = HS_bench[, k]*10^-3, col= statusCols[c('r', 'g')], lwd=1, lty = 2)
		
		# Spawners for that boostrap iteration
		lines(1:n, CI.all[[k]][[j]][[2]][, i]*10^-3, lwd = 0.8)
		
		# Benchmarks for that bootstrap  iteration
		abline(h = quantile(CI.all[[k]][[j]][[2]][, i], benchmarks)* 10^-3, col= statusCols[c('r', 'g')])
		
		if(j == 1 & i == 3) mtext(side = 4, "Naive\n(b=1)", line=1, las=1)
		# if(j == 3 & i == 3) mtext(side = 4, "Block", line=1)
		if(j == 5 & i == 3) mtext(side = 4, "Model", line=1, las=1)
		
		# if(j == 1 & i == 3) mtext(side = 4, "(b = 1)", line=1, las=1)
		if(j == 2 & i == 3) mtext(side = 4, "Block\nb = 5", line=1, las=1)
		if(j == 3 & i == 3) mtext(side = 4, "Block\nb = 10", line=1, las=1)
		if(j == 4 & i == 3) mtext(side = 4, "Block\nb = 15", line=1, las=1)
		
		if(j==1) mtext(side = 3, line = 1, paste("Example ", i))
	}
	
}	
mtext(side = 1, outer = TRUE, "Year")
mtext(side = 2, outer=TRUE, "Spawners")
dev.off()

###############################################################################
# Compare confidence intervals
###############################################################################

# Compile output
CI.compiled <- list(
	lower = array(NA, dim = c(6, 2, length(tau))),
	upper = array(NA, dim = c(6, 2, length(tau))))

for(k in 1:length(tau)){
	# First row is the true CI from re-simulated data
	CI.compiled[['lower']][1, , k] <- trueCI[, 1, k]
	CI.compiled[['upper']][1, , k] <- trueCI[, 2, k]
	
	# Last 5 rows are the different approches
	for(j in 1:5){ # for each approach
		CI.compiled[['lower']][j + 1, , k] <- CI.all[[k]][[j]][[1]][, 1]
		CI.compiled[['upper']][j + 1, , k] <- CI.all[[k]][[j]][[1]][, 2]
	} # end approach j
} # end tau

# Figure

pdf(file = "Figures/CI_comparison.pdf", width = 7, height = 3, pointsize = 10)
par(mfrow = c(1, length(tau)), mar=c(1,1,1,1), oma = c(6,4,2,0))

for(k in 1:length(tau)){
	plot(c(1, 2:4 + 0.5, 6), rep(0, 5), "n", xlim=c(0.5, 6.5), xaxt="n", xlab="", ylab="Benchmark", las=1, ylim=range(CI.compiled), yaxt = "n")
	if(k == 1){
		axis(side = 2, las =1) 
		mtext(side = 2, "Spawners", line = 3.5)
		} else axis(side = 2, labels = FALSE)
	
	for(i in 1:2){ # for S25 and S50
		polygon(x = c(0, 7, 7, 0), rep(CI.compiled[[i]][1, , k], each = 2), border=NA, col = paste(statusCols[c('r', 'g')][i], "30", sep=""))
		
		abline(h = HS_bench[i, k], col = statusCols[c('r', 'g')][i])
		
		plotCI(c(1, 2:4 + 0.5, 6)+c(-0.05, 0.05)[i], rep(HS_bench[i, k], 5), pch = NA, li = CI.compiled[[i]][2:6, 1, k], ui = CI.compiled[[i]][2:6, 2, k], gap = 0, xlim=c(-0.5, 6.5), add=TRUE, sfrac = 0.005, lwd=1.2, col = statusCols[c('r', 'g')][i])
		
	}
	
	axis(side = 1, at = c(1, 2.5, 3.5, 4.5, 6), labels = c("Naive", "b=5", "b=10", "b=15", "Model"), cex.axis=1, las=2)
	u <- par('usr')
	l1 <- u[3] - (u[4] - u[3])/4.5
	l2 <- u[3] - (u[4] - u[3])/3.5
	segments(x0 = 2, x1 = 5, y0 = l1, y1=l1, xpd=NA)
	text(3.5, l2, "Block", xpd=NA)
	
	mtext(side = 3, line = 1, substitute(paste(tau == U), list(U = tau[k])))
}
mtext(side = 1, outer=TRUE, "Approach", line=5)
dev.off()

###############################################################################
# On real data
###############################################################################

# Compare CIs from real data
# Downloaded from https://data.salmonwatersheds.ca/data-library/export.aspx
# April 5, 2020
# with parameters:
#		regions(s) = Central Coast
#		species = Chum
#		dataset = Spawner Abundance and Total Run Size for Salmon Populations
#		parameter = LGL counts
#		locations = Douglas-Gardner
#		dates = ALL
#		output type = csv file

data <- read.csv("Data/Apr_5_2020_10_28_25_QueryResults.csv")

benchmarks <- c(0.25, 0.75) # Use 0.75 upper benchmark to be consistent with Salmon Explorer

# Plot time series and benchmarks
pdf(file = "Figures/DouglasGardnerAbund.pdf", width = 5, height = 2.7, pointsize = 10)
par(mfrow = c(1,1), mar = c(4,4,1,1))
plot(data$year, data$datavalue*10^-6, "o", las = 1, cex = 0.6, xlab = "", ylab = "Estimated spawners (millions)")
abline(h = quantile(data$datavalue, benchmarks)*10^-6, col = statusCols[c('r', 'g')])
abline(h = prod(tail(data$datavalue, 4))^(1/4)*10^-6, lty = 3)
dev.off()

# Calculate CI based on different metrics
CI.dat <- list(
	naiveCI = blockBoot(series = data$datavalue, blockLength = 1, nBoot = nBoot, benchmarks = benchmarks),
	block5CI = blockBoot(series = data$datavalue, blockLength = 5, nBoot = nBoot, benchmarks = benchmarks),
	block10CI = blockBoot(series = data$datavalue, blockLength = 10, nBoot = nBoot, benchmarks = benchmarks),
	block15CI = blockBoot(series = data$datavalue, blockLength = 15, nBoot = nBoot, benchmarks = benchmarks),
	modelCI = modelBoot(series = data$datavalue, numLags = 1, nBoot = nBoot, benchmarks = benchmarks)
)

CI.dat2 <- array(NA, dim = c(5, 2, 2))
for(j in 1:5){
	CI.dat2[j, , ] <- CI.dat[[j]][[1]]
}

B <- quantile(data$datavalue, benchmarks)*10^-6
C <- prod(tail(data$datavalue, 4))^(1/4)*10^-6
ylims <- range(CI.dat2)*10^-6

pdf(file = "Figures/DouglasGardnerComparison.pdf", width = 3, height = 2.7, pointsize = 10)
par(mfrow = c(1, 1), mar = c(6,4,1,1))

plot(1, 1, "n", xlim = c(0.5, 6.5), ylim = ylims,  las = 1, xaxt="n", ylab = "Spawners (millions)", xlab = "")
polygon(x = c(0, 7, 7, 0), y = rep(c(0, B[1]), each = 2), col = statusCols['r'], border = NA)
polygon(x = c(0, 7, 7, 0), y = rep(c(B[1], B[2]), each = 2), col = statusCols['a'], border = NA)
polygon(x = c(0, 7, 7, 0), y = rep(c(3, B[2]), each = 2), col = statusCols['g'], border = NA)
abline( h = C, lty = 3)

for(i in 1:2){ # for lower and upper
		plotCI(c(1, 2:4 + 0.5, 6)+c(-0.05, 0.05)[i], rep(B[i], 5), pch = NA, li = CI.dat2[, 1, i]*10^-6, ui = CI.dat2[, 2, i]*10^-6, add=TRUE, gap = 0) # lower
}

	axis(side = 1, at = c(1, 2.5, 3.5, 4.5, 6), labels = c("Naive", "b=5", "b=10", "b=15", "Model"), cex.axis=1, las=2)
	u <- par('usr')
	l1 <- u[3] - (u[4] - u[3])/2.5
	l2 <- u[3] - (u[4] - u[3])/2
	segments(x0 = 2, x1 = 5, y0 = l1, y1=l1, xpd=NA)
	text(3.5, l2, "Block", xpd=NA)
	
dev.off()