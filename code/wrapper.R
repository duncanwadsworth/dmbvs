#' an R wrapper to C code for spike-and-slab Dirichlet--Multinomial
#' Bayesian variable selection
#'
#' @param XX covariate matrix (without intercept)
#' @param YY count matrix
#' @param intercept_variance a scalar for the prior variance on the intercept
#' of the log-linear predictors
#' @param slab_variance a scalar for the prior variance on the slab of the
#' spike-and-slab
#' @param bb_alpha a scalar for the alpha hyperparameter of the Beta-Bernoulli
#' spike inclusion prior
#' @param bb_beta a scalar for the beta hyperparameter of the Beta-Bernoulli
#' spike inclusion prior
#' @param GG the total number of MCMC iterations
#' @param thin the MCMC thinning interval
#' @param burn the number of MCMC iterations out of GG that will be discarded
#' @param proposal_alpha initial value, either a scalar or a vector of length
#' ncol(YY), if a scalar, that value is used for all proposals on alpha
#' @param proposal_beta initial value, either a scalar or a matrix with
#' ncol(XX) columns and ncol(YY) rows, if a scalar, that value is used for all
#' proposals on beta
#' @param init_alpha either a scalar or a vector of size ncol(YY)
#' @param init_beta either a scalar or a matrix with ncol(YY) rows and
#' ncol(XX) columns, inclusion initialization uses non-zero elements of
#' init_beta
#' @param exec the path to the C executable
#' @param output_location if NULL, output goes to an the directory output/
#' created in the working directory
#' @param r_seed an integer seed to pass to GSL's random number generator
#' @param selection_type either "ss" (Stochastic Search) or "Gibbs", determines
#' the MCMC mechanism for variable selection
#'
#' @return alpha: a matrix with iterations in the rows and the alphas in the
#' columns
#' @return alpha_accept: the Metropolis-Hastings acceptance ratio for the alphas
#' @return beta: a matrix with iterations in the rows and the per-iteration beta
#' matrix -- flattened by rows -- in the columns
#' @return beta_accept: the Metropolis-Hastings acceptance ratio for the betas
#' @return hyperparameters: a list containing the hyperparameters, the MCMC
#' parameters, and the data from the original function call
#'
#' @export
dmbvs = function(XX, YY, intercept_variance, slab_variance, bb_alpha, bb_beta,
                 GG, thin, burn, proposal_alpha = 0.5, proposal_beta = 0.5,
                 init_alpha = 0, init_beta = 0, exec = file.path(".","dmbvs.x"),
                 output_location = NULL, r_seed = NULL, selection_type = "Gibbs"){

  # data dimensions
  n_cats = ncol(YY)
  n_obs = nrow(XX)
  n_vars = ncol(XX)

  # argument checking
  #if(!is.integer(YY)){warning("YY is not integer typed: is it a count matrix?")}
  if(bb_alpha > bb_beta | bb_alpha == bb_beta){warning("Beta-Bernoulli prior is symmetric or left skewed: right skewed enforces sparsity")}
  if(all(XX[,1] == 1)){stop("XX appears to have an intercept column")}
  if(n_obs != nrow(YY)){stop("dimensions of XX and YY do not match")}
  if(burn > GG){stop("burnin greater than the number of iterations")}
  if(any(intercept_variance < 0, slab_variance < 0, bb_alpha < 0, bb_beta < 0)){stop("check hyperparameter values: at least one is negative")}
  if((length(proposal_beta) != 1) & (nrow(as.matrix(proposal_beta)) != n_cats) & (ncol(as.matrix(proposal_beta)) != n_vars)){stop("bad dimension for proposal_beta")}
  if(!(selection_type %in% c("ss", "Gibbs"))){stop("unrecognized variable selection type: choose ss or Gibbs")}

  # set output location
  # for extremely long runs it's necessary to leave output in a subdirectory of the
  # working directory since temp directories seems to get cleaned out occasionally
  if(is.null(output_location)){
    # for parallel simulations need unique output directories
    dir.create(paste0(getwd(), "/output-pid-", Sys.getpid()))
    out_dir = paste0(getwd(), "/output-pid-", Sys.getpid())
  }else{
    dir.create(paste0(output_location, "/output-pid-", Sys.getpid()))
    out_dir = paste0(output_location, "/output-pid-", Sys.getpid())
  }

  # text files for data
  utils::write.table(as.matrix(XX), file.path(out_dir, "covariates.txt"),
                     row.names = F, col.names = F)
  utils::write.table(as.matrix(YY), file.path(out_dir, "count_matrix.txt"),
                     row.names = F, col.names = F)

  # text files for proposals and initialization
  # intercept initialization
  if(is.null(init_alpha)){
    utils::write.table(rep(0, times = n_cats), file.path(out_dir, "init_alpha.txt"),
                       row.names = F, col.names = F)
  }else if(length(init_alpha) == 1 & is.numeric(init_alpha)){
    utils::write.table(rep(init_alpha, times = n_cats), file.path(out_dir, "init_alpha.txt"),
                       row.names = F, col.names = F)
  }else if(length(init_alpha) == n_cats & is.numeric(init_alpha)){
    utils::write.table(init_alpha, file.path(out_dir, "init_alpha.txt"),
                       row.names = F, col.names = F)
  }else if(init_alpha == "warmstart"){
    utils::write.table(scale(log(colSums(YY))), file.path(out_dir, "init_alpha.txt"),
                       row.names = F, col.names = F)
  }else{
    stop("init_alpha not recognized")
  }
  # regression parameter initialization
  if(is.null(init_beta)){
    # empty initialization
    utils::write.table(rep(0, times = n_cats * n_vars),
                       file.path(out_dir, "init_beta.txt"), row.names = F,
                       col.names = F)
  }else if(length(init_beta) == 1 & is.numeric(init_beta)){
    # scalar initialization
    utils::write.table(rep(init_beta, times = n_cats * n_vars),
                       file.path(out_dir, "init_beta.txt"), row.names = F,
                       col.names = F)
  }else if(!is.null(dim(init_beta)) & all(dim(init_beta) == c(n_cats, n_vars))){
    # matrix initialization
    utils::write.table(c(t(as.matrix(init_beta))),
                       file.path(out_dir, "init_beta.txt"), row.names = F,
                       col.names = F)
  }else if(is.character(init_beta) | init_beta == "warmstart"){
    # false discovery rate on correlation test initialization
    cormat = matrix(0, n_cats, n_vars)
    pmat = matrix(0, n_cats, n_vars)
    yy = YY/rowSums(YY) # compositionalize
    for(rr in 1:n_cats){
      for(cc in 1:n_vars){
        pmat[rr, cc] = stats::cor.test(XX[, cc], yy[, rr], method = "spearman",
                                       exact = F)$p.value
        cormat[rr, cc] = stats::cor(XX[, cc], yy[, rr], method = "spearman")
      }
    }
    # defaults to 0.2 false discovery rate
    pm = matrix((stats::p.adjust(c(pmat), method = "fdr") <= 0.2) + 0, n_cats,
                n_vars)
    betmat = cormat * pm
    utils::write.table(c(t(betmat)), file.path(out_dir, "init_beta.txt"),
                       row.names = F, col.names = F)
  }else{
    stop("init_beta not recognized, if using a matrix, ensure proper dimensions")
  }

  # intercept proposal
  if(length(proposal_alpha) == 1){
    utils::write.table(rep(proposal_alpha, times = n_cats),
                       file.path(out_dir, "proposal_alpha.txt"), row.names = F,
                       col.names = F)
  }else if(length(proposal_alpha) == n_cats){
    utils::write.table(proposal_alpha, file.path(out_dir, "proposal_alpha.txt"),
                       row.names = F, col.names = F)
  }else{
    stop("problem with proposal_alpha")
  }
  # regression proposal
  if(length(proposal_beta) == 1){
    utils::write.table(rep(proposal_beta, times = n_cats * n_vars),
                       file.path(out_dir, "proposal_beta.txt"), row.names = F,
                       col.names = F)
  }else if((nrow(proposal_beta) == n_cats) & (ncol(proposal_beta) == n_vars)){
    utils::write.table(c(t(proposal_beta)), file.path(out_dir, "proposal_beta.txt"),
                       row.names = F, col.names = F)
  }else{
    stop("problem with proposal_beta")
  }

  # pass the external seed to GSL which will have its own corresponding seed (see the .out file)
  if(is.null(r_seed)){
    a_random_seed = sample(1e6, 1)
  }else if(is.integer(r_seed)){
    a_random_seed = r_seed
  }else{
    stop("problem with external seed to pass to GSL")
  }
  # run compiled code
  if(selection_type == "ss"){
    command = paste(exec, GG, thin, burn, intercept_variance, slab_variance,
                    bb_alpha, bb_beta, n_cats, n_obs, n_vars, a_random_seed,
                    out_dir, 0, ">", paste0("dmbvs-pid-", Sys.getpid(), ".out"))
  } else if(selection_type == "Gibbs"){
    command = paste(exec, GG, thin, burn, intercept_variance, slab_variance,
                    bb_alpha, bb_beta, n_cats, n_obs, n_vars, a_random_seed,
                    out_dir, 1, ">", paste0("dmbvs-pid-", Sys.getpid(), ".out"))
  }
  system(command)
  # read output
  aa = utils::read.table(file.path(out_dir, "alpha.out"))
  aaa = utils::read.table(file.path(out_dir, "alpha_acceptance.out"))
  # bb is read as a n_vars * n_cats long list with each element have n_iters
  # so when it's unlist()ed you get the first variable for n_iters, then the
  # second, etc
  bb = utils::read.table(file.path(out_dir, "beta.out"))
  bba = utils::read.table(file.path(out_dir, "beta_acceptance.out"))

  # variables in the columns, iterations in the rows
  return(list(alpha = t(matrix(unlist(aa), nrow = n_cats, byrow = T)),
              alpha_accept = aaa[,1],
              # after transpose beta has iterations in the rows
              beta = t(matrix(unlist(bb), nrow = (n_cats * n_vars), byrow = T)),
              beta_accept = bba[,1],
              hyperparameters = list(mcmc = data.frame(GG = GG, thin = thin, burn = burn,
                                                       proposal_alpha = proposal_alpha,
                                                       proposal_beta = proposal_beta,
                                                       random_seed = a_random_seed),
                                     priors = data.frame(intercept_variance = intercept_variance,
                                                         slab_variance = slab_variance,
                                                         bb_alpha = bb_alpha, bb_beta = bb_beta),
                                     inputdata = list(XX = XX, YY = YY))))
}

