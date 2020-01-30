library(boot)
library(DistributionUtils)



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


HS_benchCI_block <- blockBoot(series = spawners[2:(n+1)], blockLength = 50, nBoot = 10000)

HS_benchCI_block
HS_benchCI_naive

# quartz(width = 6.3, height = 7, pointsize = 10)
plot(1:n, spawners[2:(n+1)], "o", ylim=c(0, 10000))
abline(h = quantile(spawners[2:(n+1)], c(0.25, 0.5)), col = statusCols[c('r', 'g')], lwd=2)
# Naive bootstrap
polygon(x = c(-10, -10, 60, 60), y = c(HS_benchCI_naive[, 1], HS_benchCI_naive[2:1,1]), col=paste(statusCols['r'], "30", sep=""), border=NA)
polygon(x = c(-10, -10, 60, 60), y = c(HS_benchCI_naive[, 2], HS_benchCI_naive[2:1,2]), col=paste(statusCols['g'], "30", sep=""), border=NA)

# Block bootstrap
abline(h = HS_benchCI_block[,1], lty = 2, col=statusCols['r'])
abline(h = HS_benchCI_block[,2], lty = 2, col=statusCols['g'])

###############################################################################
# What blockLength to use?
###############################################################################

blockLengths <- c(1, 2, 5, 8, 10, 12, 15, 18, 20, 25, 30, 40)
HS_benchCI_blockTrial <- list(); length(HS_benchCI_blockTrial) <- length(blockLengths)

for(i in 1:length(blockLengths)){
	HS_benchCI_blockTrial[[i]] <- blockBoot(series = spawners[2:(n+1)], blockLength = blockLengths[i], nBoot = 10000)
}

quartz(width = 6.3, height = 7, pointsize = 10)
par(mfrow = c(4,3), mar=c(1,1,0,0), oma=c(4,4,2,1))

for(i in 1:length(blockLengths)){
	plot(1:n, spawners[2:(n+1)], "o", ylim=c(0, 10000), xaxt="n", yaxt="n", col=grey(0.8))
	axis(side = 1, labels = FALSE)
	axis(side = 2, labels = FALSE)
	abline(h = quantile(spawners[2:(n+1)], c(0.25, 0.5)), col = statusCols[c('r', 'g')], lwd=2)
	# Naive bootstrap
	polygon(x = c(-10, -10, 60, 60), y = c(HS_benchCI_naive[, 1], HS_benchCI_naive[2:1,1]), col=paste(statusCols['r'], "30", sep=""), border=NA)
	polygon(x = c(-10, -10, 60, 60), y = c(HS_benchCI_naive[, 2], HS_benchCI_naive[2:1,2]), col=paste(statusCols['g'], "30", sep=""), border=NA)
	
	# Block bootstrap
	abline(h = HS_benchCI_blockTrial[[i]][,1], lty = 2, col=statusCols['r'])
	abline(h = HS_benchCI_blockTrial[[i]][,2], lty = 2, col=statusCols['g'])
	
	mtext(side = 3, adj=0, line = -1.5, paste(" ", letters[i], ") b = ", blockLengths[i], sep=""))
}


# Different appraoch
blockLengths <- c(1:(n/2))
HS_benchCI_blockTrial2 <- matrix(NA, nrow = length(blockLengths), ncol = 4)

for(i in 1:length(blockLengths)){
	HS_benchCI_blockTrial.i <- blockBoot(series = spawners[2:(n+1)], blockLength = blockLengths[i], nBoot = 10000)
	
	HS_benchCI_blockTrial2[i, 1] <- HS_benchCI_blockTrial.i[1,1]
	HS_benchCI_blockTrial2[i, 2] <- HS_benchCI_blockTrial.i[1,2]
	HS_benchCI_blockTrial2[i, 3] <- HS_benchCI_blockTrial.i[2,1]
	HS_benchCI_blockTrial2[i, 4] <- HS_benchCI_blockTrial.i[2,2]
}

par(mfrow = c(1,2), mar=c(4,4,1,0), oma=c(0,1,1,1), mgp=c(3.5, 1, 0))

# S25
plot(blockLengths, HS_benchCI_blockTrial2[, 1], "o", col = statusCols['r'], cex=0.5, las=1, ylim=range(HS_benchCI_blockTrial2[, 1:2]), ylab="", xlab="")
points(blockLengths, HS_benchCI_blockTrial2[, 2], "o", col = statusCols['r'], cex=0.5)
abline(h = quantile(spawners[2:(n+1)], 0.25), col = paste(statusCols['r'], 50, sep=""), lwd=3)
mtext(side = 3, adj = 0, line = 0.5, expression(paste("a) ", italic(S)[25])))
mtext(side = 2, line = 3.5, "Confidence interval on benchmarks")
legend("topright", lty = c(1,1), pch = c(NA, 1), pt.cex = c(NA, 0.5), col = c(paste(statusCols['r'], 50, sep=""), statusCols['r']), lwd=c(3, 1), c("Benchmark", "95% CI"), bty="n")
abline(h = HS_benchCI_model[[1]][, 1], col = statusCols['r'], lty=3)

# S50
plot(blockLengths, HS_benchCI_blockTrial2[, 3], "o", col = statusCols['g'], cex=0.5, las=1, ylim=range(HS_benchCI_blockTrial2[, 3:4]), ylab="", xlab="")
points(blockLengths, HS_benchCI_blockTrial2[, 4], "o", col = statusCols['g'], cex=0.5)
abline(h = quantile(spawners[2:(n+1)], 0.5), col = paste(statusCols['g'], 50, sep=""), lwd=3)
mtext(side = 3, adj = 0, line = 0.5, expression(paste("b) ", italic(S)[50])))
abline(h = HS_benchCI_model[[1]][, 2], col = statusCols['g'], lty=3)

legend("topright", lty = c(1,1), pch = c(NA, 1), pt.cex = c(NA, 0.5), col = c(paste(statusCols['g'], 50, sep=""), statusCols['g']), lwd=c(3, 1), c("Benchmark", "95% CI"), bty="n")

mtext(side = 1, outer = TRUE, "Block length (b)", line = -1)
