# mmln: Mixed Multinomial Logistic Normal modeling
## Last Update: March 2023 (needs work)
Applied to forecasting baseball batting outcomes for players transitioning between NPB and MLB
For paper based on this work, and to cite usage, please see:
Gerber, E., Craig, B. (2021). A mixed effects multinomial logistic-normal model for forecasting baseball performance. *Journal of Quantitative Analysis in Sports, 17*(3), 221-239. https://doi.org/10.1515/jqas-2020-0007.

## mmln_helpers.R
Contains functions necessary for running the algorithms defined in mmln.R

## mmln.R
Contains two main sets of functions
 - mixgibbs
     - the mmln function (need to rename and clean up/make more efficient): fits by default a mixed effects multinomial logistic normal model (can fit a fixed effects, ignoring specified random effects matrix, with option `mixed=FALSE`)
     - arguments defined in comments (should move them here...)
 - multmixgibbs_b0prior
     - a Bayesian mixed effects multinomial logistic regression model, with proposal covariance matrix based on Polya-Gamma auxiliary variables

Both main functions (above) have arguments used for checking simulation studies (for defining "true" values of the various parameters) which are not used by default. They also both have log-likelihood functions that take the MCMC chains and allow plotting of the likelihood over the iterations to assess convergence. Related functions for using the chains to build the predictions are included in the mmln_predict.R file.

## mmln_predict.R
Contains functions necessary for estimating random effects and making predictions for out-of-sample (test) data based on model runs.

## yoshida_run_march23.R
The most recent run/work with the model; ran at the beginning of the 2023 season to generate predictions for Mastaka Yoshida, and the several players who transitioned from MLB to NPB for the first time.

## Data Files
### NZPAcrosslist_Final.csv
Updated at the end of the 2022 season with the 2022 season stats for all players still playing in MLB and NPB to have played in both.

### NZPAcrosslist_Final_test.csv
Updated at the beginning of the 2023 season with the career stats for all players making the transition between the two leagues. Includes at least a couple "fake guys" (players who are needed only because of the inefficient way I coded cleaning the data; i.e. if there are no players making the transition to the AL, will need to make a line in the test set for a player who is). These players are removed before performing any actual analysis, i.e. before estimating player random effects.
