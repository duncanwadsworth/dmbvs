
###################################################################
# Example Analysis of simulated data using a Dirichlet-Multinomial
# Bayesian variable selection model
###################################################################

# -- README -------------------------------------------------------------------
# This script gives an example of how to run the code on a simulated
# dataset. It assumes several things:
# 1. The Gnu Scientific Library is installed
# 2. The paths in the Makefile are correct
# 3. The data is available
# -----------------------------------------------------------------------------

# -- simulate data ------------------------------------------------------------
source(file.path("code", "helper_functions.R"))
simdata = simulate_dirichlet_multinomial_regression(n_obs = 100, n_vars = 50,
                                                    n_taxa = 50, n_relevant_vars = 5,
                                                    n_relevant_taxa = 5)
Ysim = simdata$YY
Xsim = simdata$XX

# -- run analysis -------------------------------------------------------------
# preliminaries
source(file.path("code", "wrapper.R"))
system("cd code; make")
executable_location = "code/dmbvs.x"
save_prefix = "simulation"

# prepare and check data
YY = as.matrix(Ysim)
XX = scale(as.matrix(Xsim[,-1]), center = T, scale = T)
colnames(YY) = paste0("taxa", 1:ncol(YY))
colnames(XX) = paste0("covariate", 1:ncol(XX))
dim(YY)
dim(XX)

# MCMC and hyperparameters
#GG = 301L; thin = 2L; burn = 101L; # for testing
GG = 11001L; thin = 10L; burn = 1001L;
cat("number of kept iterations:", (GG - burn)/thin, "\n")
bb_alpha = 0.02; bb_beta = 2 - bb_alpha
cat("Beta-Binomial mean:", bb_alpha/(bb_alpha + bb_beta), "\n")
proposal_alpha = 0.5; proposal_beta = 0.5
slab_variance = 10; intercept_variance = 10

# run the algorithm
results = dmbvs(XX = XX, YY = YY, intercept_variance = intercept_variance,
                slab_variance = slab_variance, bb_alpha = bb_alpha,
                bb_beta = bb_beta, GG = GG, thin = thin, burn = burn,
                init_beta = "warmstart", init_alpha = "warmstart",
                proposal_alpha = proposal_alpha, proposal_beta = proposal_beta,
                exec = executable_location,
                output_location = ".")
params = data.frame(GG, burn, thin, intercept_variance,
                    slab_variance, bb_alpha, bb_beta,
                    proposal_alpha, proposal_beta)
save(results, params, XX, YY,
     file = paste0("results-", save_prefix, "-", Sys.Date(), ".RData"))

# quick check results
source(file.path("code", "helper_functions.R"))
mppi = colMeans((results$beta != 0) + 0)
(blfdrate = bfdr(mppi, threshold = 0.1)$threshold)
MPPI = data.frame(expand.grid(covariates = colnames(results$hyperparameters$inputdata$XX),
                              taxa = colnames(results$hyperparameters$inputdata$YY)),
                  mppi = mppi,
                  beta = colMeans(results$beta),
                  truebeta = c(t(simdata$betas[,-1])))
subset(MPPI, mppi > blfdrate | truebeta != 0)
plot(mppi, type = "h", ylab = "MPPI",
     xlab = "beta index", main = "Manhattan plot")

# active variable traceplot
plot.ts(rowSums((results$beta != 0) + 0), main = "Active variables traceplot",
        ylab = "number of betas in the model", xlab = "iteration")

# some of the selected beta traceplots
selected = which(mppi > 0.5)
fortraces = selected[sample(length(selected), 10)]
plot.ts(results$beta[,fortraces], main = "Some selected beta traceplots",
        xlab = "iteration", ylab = "")
