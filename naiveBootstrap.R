# Code to bootstrap percentile benchmarks for status assessments
# 3 May 2017
require(bootstrap)
S<-rnorm(20)  #Time-series of observed spawners
percentiles<-function(x){quantile(x, probs=c(0.25,0.75))}	 #Function to calculate 25th percentile
benchmarks<-percentiles(S)	# S25th and S75th benchmarks
results <- bootstrap(S,nboot=1000,theta=percentiles)$thetastar	#Non-parametric bootstrap values
CIs.S25th<-quantile(results[1,], probs=c(0.025, 0.975))	# Upper and lower CIs for S25th
CIs.S75th<-quantile(results[2,], probs=c(0.025, 0.975))  # Upper and lower CIsfor S75th
benchmarks
CIs.S25th
CIs.S75th