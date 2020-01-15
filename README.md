Constructing confidence intervals on Historical Spawners benchmarks
================
Stephanie Peacock
2020-01-15

This project aims to compare different methods to constructing 95% confidence intervals on Historical Spawners (HS) benchmarks for inclusion in the [Pacific Salmon Explorer](www.salmonexplorer.ca).

Background
----------

To date, the assessments of biological status in the Salmon Explorer have provided confidence intervals for the Stock-Recruitment (SR) benchmarks, but not for the HS benchmarks. The HS benchmarks are calculated as the 25th and 50th percentiles of historical spawner abundance for the CU. (Note that older assessments under the Pacific Salmon Explorer considered the upper HS benchmark to be the 75th percentile of historical spawner abundance.)

In their study of data-limited Chum-salmon Conservation Units (CUs) in Southern BC, [Holt et al. (2018)](http://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2018/2018_011-eng.html) constructed confidence intervals on HS benchmarks using naive bootstrapping "by resampling the time-series with replacement to generate a distribution of lower and upper benchmarks". However, they recognized that this approach may may over-estimate confidence intervals if time-series are autocorrelated, which is likely the case for time-series of spawner abundances. They recommended that methods that account for temporal autocorrelation should be considered in the future.

Approach
--------

We consider two alternative approaches to generating confidence intervals on the HS benchmarks that account for autocorrelation in the time-series of spawner abundances. The first is a model-based approach and the second is known as "block bootstrapping".
