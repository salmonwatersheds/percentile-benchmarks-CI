Constructing confidence intervals on Historical Spawners benchmarks
================
Stephanie Peacock
2020-01-15

This project aims to compare different methods to constructing 95% confidence intervals on Historical Spawners (HS) benchmarks for inclusion in the [Pacific Salmon Explorer](www.salmonexplorer.ca).

Background
----------

To date, the assessments of biological status in the Salmon Explorer have provided confidence intervals for the Stock-Recruitment (SR) benchmarks, but not for the HS benchmarks. The HS benchmarks are calculated as the 25th and 50th percentiles of historical spawner abundance for the CU. (Note that older assessments under the Pacific Salmon Explorer considered the upper HS benchmark to be the 75th percentile of historical spawner abundance.)

![Fig. 1: Screen shot of the biological status assessment of chum salmon CU Douglas-Gardner from the Pacific Salmon Explorer from April 5, 2020. Note that the spawner-recruitment benchmarks have confidence intervals, but the historic(al) spawners benchmarks do not.](Figures/SalmonExplorer_screenShotApr52020.tiff)

In their study of data-limited Chum-salmon Conservation Units (CUs) in Southern BC, [Holt et al. (2018)](http://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2018/2018_011-eng.html) constructed confidence intervals on HS benchmarks using naive bootstrapping "by resampling the time-series with replacement to generate a distribution of lower and upper benchmarks". However, they recognized that this approach may may over-estimate confidence intervals if time-series are autocorrelated, which is likely the case for time-series of spawner abundances. They recommended that methods that account for temporal autocorrelation should be considered in the future.

Methods
-------

We consider two alternative approaches to generating confidence intervals on the HS benchmarks that account for autocorrelation in the time-series of spawner abundances. The first is known as "block bootstrapping" and the second is a model-based approach. We compare these two approaches along with the naive bootstapping applied by [Holt et al. (2018)](http://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2018/2018_011-eng.html).

### Simulating data

To quantify the difference between approaches with increasing temporal autocorrelation in spawner abundance, we tested the different approaches for computing confidence intervals on simulated data. We simulated a time-series of spawner abudnances using the Ricker equation, including temporal autocorrelation in residuals. The parameters for the simulated data were:

``` r
a <- 1.4                            # productivity
b <- 0.0001                     # density dependence
sigma <- 0.1                # variance in recruitment deviates
tau <- 0.6                      # temporal autocorrelation in recruitment deviates
n <- 50                             # number of years of data to simulate 
```

We also simulated harvest to have a realistic time-series. We assumed a constant target harvest rate with beta-distributed error around the realized harvest rate each year:

``` r
targetHarvest <- 0.42                           # target harvest rate
sigmaHarvest <- 0.13                            # standard deviation in realized harvest rate

# Draw n realized harvest rates, with mean targetHarvest:
beta1 <- (targetHarvest^2 - targetHarvest^3 - sigmaHarvest^2*targetHarvest)/(sigmaHarvest^2)
beta2 <- (targetHarvest * (1 - targetHarvest)^2 - sigmaHarvest^2*(1 - targetHarvest))/(sigmaHarvest^2)
harvestRate <- rbeta(n = n, shape1 = beta1, shape2 = beta2)
```

We initialized the simualted time-series using a spawner abundance that was 20% of the equilibrium spawner abundance, and then simulated the dynamics of recruitment and harvest over the next 50 years:

``` r
# Simulate spawner data:
spawners <- numeric(n + 1)          # vector to store spawner abundance
spawners[1] <- 0.2 * (a / b)        # initiate spawner abundance at 20% Seq
    
phi <- numeric(n + 1)                       # vector to store recruitment deviates
phi[1] <- rnorm(1, 0, sigma)        # initiate recruitment deviates assuming phi[0] = 0
recruits <- numeric(n + 1)          # vector to store recruits

for(i in 1:n){
    phi[i + 1] <- tau * phi[i] + rnorm(1, 0, sigma) 
    recruits[i + 1] <-spawners[i] * exp(a - b * spawners[i]) * exp(phi[i + 1])
    spawners[i + 1] <- (1 - harvestRate[i]) * recruits[i + 1]
}
```

### Calculating benchmarks

The HS benchmarks for the simulated time-series are:

``` r
benchmarks <- c(0.25, 0.5)

# Remove initialization year when calculating benchmarks:
series <- spawners[2 : (n + 1)] 

HS_bench <- quantile(series, benchmarks)
AvgS <- prod(tail(series, 4))^(1/4)

HS_bench
```

    ##      25%      50% 
    ## 6096.471 7529.093

![](README_files/figure-markdown_github/unnamed-chunk-5-1.png)

### Calculating confidence intervals

We applied three different types of bootstrapping approaches: (1) naive boostrap, not accounting for temporal autocorrelation, (2) block bootstrap, which samples blocks of the original timeseries to include the autocorrelation structure in the permuted timeseries, and (3) model-based, which fits a model to estimate the magnitude of autocorrelation and then resamples (with replacement) from the fitted residuals to simulate a new dataset with temporal autocorrelation.

#### Naive bootstrapping

Naive bootstrapping does not account for the autocorrelation in the timeseries. The spawner abundances are simply sampled with replacement, and the HS benchmarks calculated from the re-sampled data. This is repeated many (thousands) times and then the 2.5% and 97.5% quantiles of the resulting bootsrapped benchmarks give the 95% confidence interval. This "naive bootstrapping" approach was applied by [Holt et al. (2018)](http://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2018/2018_011-eng.html), who recognized that it may overestimate confidence intervals if time-series are autocorrelated.

``` r
nBoot <- 10000 # number of bootstrap iterations

HS_benchBoot <- matrix(NA, nrow = nBoot, ncol = 2, dimnames = list(c(1:nBoot), c("lower", "upper")))

for(i in 1:nBoot){
    series.i <- sample(series, replace = TRUE) # re-sample data
    HS_benchBoot[i, ] <- quantile(series.i, benchmarks) # save benchmark
}

# Calculate 2.5% and 97.5% quantiles of bootstrapped benchmarks
HS_benchCI_naive <- apply(HS_benchBoot, 2, quantile, c(0.025, 0.975))
HS_benchCI_naive
```

    ##          lower    upper
    ## 2.5%  5176.684 6944.610
    ## 97.5% 7028.533 8293.627

### Block boostrapping

Block bootstrapping aims to avoid bias in the estimates from naiive bootstrapping by sampling the original timeseries in "blocks", thus preserving (at least in part) the temporal autocorrelation in the original timeseries. There are different ways to choose the blocks being sampled; we adopted one of the simplest approaches that allows blocks to overlap (i.e., the "moving block" approch; [see here](http://anson.ucdavis.edu/~peterh/sta251/bootstrap-lectures-to-may-16.pdf)). We kept block length constant, but randomly chose the starting point for each block between one and `n - blockLength + 1`, where `n` is the length of the timeseries and `blockLength` is the length of the block. Determining the optimal block length is one of the cruxes of block bootstrapping, and we investigated the outcome under different block lengths for comparison (see below). For illustration, we use a block length of 10 here. If the length of series `n` did not divide evenly by `blockLength`, we used the first `1:n` elements of the new series when calculating benchmarks.

``` r
blockLength <- 10
nBlocks <- ceiling(n/blockLength) # Number of blocks to be created

# Matrix to store bootstrapped values of CI
HS_benchBoot2 <- matrix(
        NA, 
        nrow = nBoot, 
        ncol = 2, 
        dimnames = list(c(1:nBoot), c("lower", "upper")))
    
# Block bootstrap:
for(i in 1:nBoot){ # for each iteration
    
    # Select starting points for blocks
    j <- sample(1:(n - blockLength + 1), 
                            size = nBlocks, replace = TRUE)
    
    # Indices of new time series to be constructed
    newIndex <- rep(j, each = blockLength) + rep(0:(blockLength - 1), nBlocks)
    
    # Block-sampled spawner timeseries
    series.i <- series[newIndex[1:n]]
    
    HS_benchBoot2[i, ] <- quantile(series.i, benchmarks)
} # end bootstrap loop
    
HS_benchCI_block <- apply(HS_benchBoot2, 2, quantile, c(0.025, 0.975))
HS_benchCI_block
```

    ##          lower    upper
    ## 2.5%  5746.770 7005.276
    ## 97.5% 7217.146 8396.342

We created the function `blockBoot` (see Appendix) return the bootstrapped confidence intervals on the historical spawners benchmarks for a timeseries of spawner abundances. Note that setting `blockLength = 1` will yield the naive bootstrap result from the previous subsection.

### Model-based approach

The model-based approach is another form of bootstrapping, but rather than resampling the from time series, we first fit a simple model to the time series of log spawner abundances to estimate the autocorrelation. We calculate the temporal autocorrelation for a lag of one year (`numLags = 1`), but further consideration may want to be given to including other lags.

``` r
numLags <- 1

# Fit model to estimate autocorrelation
ar.fit <- ar.mle(
    log(series), # spawner time series
    demean = TRUE, # Estimate mean spawners
    intercept = FALSE, # Intercept = 0 for autocorrelation
    order.max = numLags, # lag for autocorrelation
    aic = FALSE # estimate autocorrelation for all numLags 
) # standard OLS
```

Next, we extract the residuals from the fitted series and generate `nBoot` new timeseries of residuals, `res.star`.

``` r
# Matrices to store bootstrapped residuals and spawners (obs)
    res.star <- matrix(nrow = n, ncol = nBoot)
    res.star[, ] <- sample(na.omit(ar.fit$resid), n * nBoot, replace = TRUE)
```

Finally, we simulate new timeseries centered on the mean spawner abundance (`ar.fit$x.mean`), using the boostrapped residuals and estimatd autocorrelation.

``` r
# Matrix to store (log) simulated spawner abundances
obs.star.log <- matrix(nrow = n + numLags, ncol = nBoot) 
    
# Matrix to store bootstrapped values of CI
HS_benchBoot3 <- matrix(
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
    
    HS_benchBoot3[i, ] <- quantile(exp(obs.star.log[(numLags + 1):(numLags + n), i]), benchmarks)
} # end bootstrap loop

HS_benchCI_model <- apply(HS_benchBoot3, 2, quantile, c(0.025, 0.975))

obs.star <- exp(tail(obs.star.log, n))
```

We created the function `modelBoot` to return the model-bootstrapped confidence intervals on the historical spawners benchmarks for a timeseries of spawner abundances (see Appendix).

### Comparing approaches

#### Simulated data

We compared the two approaches to calculating confidence intervals (block boostrap and model-based) to the "true" confidence intervals in the upper and lower historical-spawners benchmarks. We calculated the block boostrap using block lengths of 1 (i.e., naiive boostrap), 5, 10, and 15 years. This yielded five different sets of confidence intervals being compared (block lengths 1, 5, 10, 15, and model-based).

We estimated the "true" confidence intervals by simulating 10,000 independent timeseries of spawner abundance using the same underlying Ricker parameters and calculating the HS benchmarks of each. The "true" confidence intervals were the 2.5th and 97.5th quantiles of these resulting 10,000 HS benchmarks.

The temporal autocorrelation in recruitment deviates may affect the relative performance of the different appraoches, and so we compared confidence intervals under autocorrelation coefficients (*τ*) of 0.1 (low), 0.3, 0.6, and 0.9 (high).

R code for simulations and comparison plots can be found in `comparison.R`.

#### True data

We also calculated the five different confidence intervals on HS benchmarks from estimated spawner abundance from a real central-coast chum salmon Conservation Unit: Douglas-Gardner. According to the [Pacific Salmon Explorer](https://salmonexplorer.ca/#!/central-coast/chum/douglas-gardner&to=population&pop=BENCHMARK_STATUS&pop-detail=1), this CU has a current spawner abundance of 152,789 and an amber status under the HS benchmarks.

The true data analysis could be easily applied to other CUs; see the code in `comparison.R`.

### Results

#### Simulated data

![Fig. 2: Three examples of bootstrapped time series (black lines) under the different approaches to calculating CIs. The original (true) time series is in grey in each panel. Horizontal solid red and green lines are the lower and upper HS benchmarks of 25th and 50th percentile of historical spawner abundance, and the dashed horizontal lines are the corresponding benchmarks from the bootstrapped time series. Here, just three examples are shown, but to calculate confidence intervals on benchmarks, this would be repeated thousands of times and the 2.5% and 97.5% quantiles of the resulting bootstrapped benchmarks (dashed lines) would give the confidence interval.](Figures/BootstrapSeries.pdf)

#### True data

### Conclusion

The model-based confidence intervals most closely matched the "true" confidence intervals from the simulated data. It is therefore recommended that the model-based approach, which accounts for autocorrelation in the time series, be used when calcualting confidence intervals on Historical Spawners benchmarks.

The model-based approach we applied is very simple, estimating just the mean abundance and lag-1 autocorrelation in the time series. Further work may incorporate life-history details and/or consider different time-lags for the autocorrelation used in simulating time series.

Appendix: Functions for calculating confidence intervals
--------------------------------------------------------

The functions below are in the functions.R file, which also contains two functions used for simulating the spawner-recruit data.

``` r
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
```

``` r
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
```
