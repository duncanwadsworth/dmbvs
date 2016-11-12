#' calculate the Bayesian False Discovery Rate
#'
#' @param mppi_vector A vector of marginal posterior probabilities of inclusion.
#' @param threshold The expected false discovery rate threshold
#'
#' @return selected: A boolean vector of selected (= T) and rejected (= F)
#' variables
#' @return threshold: The BFDR threshold
#'
#' @references
#' Newton, M. A., Noueiry, A., Sarkar, D., & Ahlquist, P. (2004). Detecting
#' differential gene expression with a semiparametric hierarchical mixture
#' method. Biostatistics, 5(2), 155-76. doi:10.1093/biostatistics/5.2.155
#'
#' @export
bfdr = function(mppi_vector, threshold = 0.1){
  # arg checking
  if(any(mppi_vector > 1 | mppi_vector < 0)){
    stop("Bad input: mppi_vector should contain probabilities")
  }
  if(threshold > 1 | threshold < 0){
    stop("Bad input: threshold should be a probability")
  }
  # pretty simple
  onek = ((1 - mppi_vector) < threshold) + 0
  # possible to select none
  if(sum(onek) == 0){
    selected = rep(FALSE, length(onek))
    thecut = 0
  }else{
    thecut = sum((1 - mppi_vector) * onek)/sum(onek)
    selected = (1 - mppi_vector) < thecut
  }
  return(list(selected = selected, threshold = 1 - thecut))
}

#' a bipartite graph showing the association between X and Y
#'
#' @note requires the igraph package
#'
#' @param MPPI a data.frame with n_vars x n_taxa rows and these four columns:
#'  1) covariate: the names of the columns of X
#'  2) taxa: the names of the columns of Y
#'  3) mppi: the marginal posterior probability of inclusion for each taxa
#'          by covariate parameter
#'  4) beta: a point estimate of for each taxa by covariate parameter
#' @param mppi_threshold the threshold for inclusion in the plot
#' @param inc_legend boolean to include the legend
#' @param lwdx a scalar multiplier for growing or shrinking the widths of edges
#' @param graph_layout either "circular" for a round layout or "bipartite" for 
#' a side-by-side layout
#' @param lab_dist a scalar argument for the distance between node centers and 
#' the node labels
#' @param ... passthrough arguments
#'
#' @return a plot
#'
#' @export
association_plot = function(MPPI, mppi_threshold = 0.5, inc_legend = F,
                            lwdx = 5, graph_layout = "circular", lab_dist = 2, 
                            ...){
  if(any(!(colnames(MPPI) %in% c("covariates", "taxa", "mppi", "beta")))){
    stop("please ensure the column names of MPPI are correct\n*** they must be: covariates, taxa, mppi, beta")
  }
  if(!require(igraph)){
    stop("please ensure the igraph package is installed")
  }
  # first munge the data into a format easily read by graph_from_data_frame()
  mm = subset(MPPI, mppi > mppi_threshold)
  if(nrow(mm) == 0){ # when there are no relevant associations
    graphics::plot(1, type = "n", axes = F, xlab = "", ylab = "",
                   xlim = c(-1, 1), ylim = c(-1, 1), ...)
    graphics::text(0, 0, "No Associations\nAbove MPPI\nThreshold")
  }else{ # when there are relevant associations
    mm$esign = round((sign(mm$beta) + 2)/3) + 1
    mm$emag = abs(mm$beta)
    mppi = mm[order(mm$covariates, mm$taxa),]
    # readin data.frame
    hh = igraph::graph_from_data_frame(mppi, directed = F)
    # modify graph characteristics
    igraph::E(hh)$weight = igraph::E(hh)$emag*lwdx
    #igraph::E(hh)$lty = igraph::E(hh)$esign
    igraph::E(hh)$lty = 1
    igraph::E(hh)$color = ifelse(igraph::E(hh)$esign == 2, 1, 2)
    igraph::V(hh)$type = c(rep(T, times = length(unique(mppi$covariates))),
                           rep(F, times = length(unique(mppi$taxa))))
    # layout as a bipartite graph
    if(graph_layout == "bipartite"){
      la = layout.bipartite(hh, hgap = 20)
      graphics::plot(hh, layout = la[,c(2,1)],
                     edge.width = igraph::E(hh)$weight,
                     edge.lty = igraph::E(hh)$lty,
                     edge.color = c("red", "blue")[igraph::E(hh)$color],
                     vertex.color = c("green", "yellow")[igraph::V(hh)$type + 1],
                     vertex.shape = c("circle", "square")[igraph::V(hh)$type + 1],
                     vertex.frame.color = "black", vertex.label.cex = 1.2,
                     vertex.label.color = "black", vertex.label.family = "sans",
                     vertex.label.font = 2,
                     vertex.label.degree = c(2*pi, pi)[igraph::V(hh)$type + 1],
                     vertex.label.dist = lab_dist, ...)
    }
    # layout as a circular graph
    if(graph_layout == "circular"){
      la = igraph::layout.circle(hh)
      # from http://stackoverflow.com/questions/23209802/placing-vertex-label-outside-a-circular-layout-in-igraph
      radian.rescale = function(x, start = 0, direction = 1) {
        c.rotate = function(x) (x + start) %% (2 * pi) * direction
        c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
      }
      lab.locs = radian.rescale(x = 1:nrow(la), start = 0, direction = -1)
      browser()
      graphics::plot(hh, layout = la,
                     edge.width = igraph::E(hh)$weight,
                     edge.lty = igraph::E(hh)$lty,
                     edge.color = c("red", "blue")[igraph::E(hh)$color],
                     vertex.color = c("green", "yellow")[igraph::V(hh)$type + 1],
                     vertex.shape = c("circle", "square")[igraph::V(hh)$type + 1],
                     vertex.frame.color = "black", vertex.label.cex = 1.2,
                     vertex.label.color = "black", vertex.label.family = "sans",
                     vertex.label.font = 2, vertex.label.dist = lab_dist,
                     vertex.label.degree = lab.locs)
    }
    # a black and white color scheme without vertex shapes
    # graphics::plot(hh, layout = la,
    #                edge.width = igraph::E(hh)$weight,
    #                edge.color = "black",
    #                vertex.color = c("grey", "white")[V(hh)$type + 1],
    #                vertex.label.cex = 1.3,
    #                vertex.label.color = "black", vertex.label.family = "sans",
    #                vertex.label.font = c(2, 3)[igraph::V(hh)$type + 1],
    #                vertex.size = 10, vertex.label.dist = 0.5,
    #                vertex.label.degree = lab.locs)
    if(inc_legend){
      #legend("topright", legend = c("positive", "negative"), lty = c(1, 2), col = "grey50", lwd = 3)
      graphics::legend("topright", legend = c("positive", "negative"), lwd = c(3, 3), col = c("blue","red"))
    }
  }
}

#' abundance table visualization
#'
#' @note requires the ggplot2 package
#'
#' @param count_matrix a matrix of integers
#' @param title an optional title
#'
#' @return a plot
#'
#' @export
abundance_plot = function(count_matrix, title = ""){
  if(!require(ggplot2)){
    stop("please ensure the ggplot2 package is installed")
  }
  # sort matrix by column means - make pretty
  ss = count_matrix[, order(colMeans(count_matrix), decreasing = T)]
  # coordinates
  xyz = cbind(expand.grid(1:dim(ss)[1], 1:dim(ss)[2]), as.vector(ss), as.vector(ss) > 0)
  names(xyz) = c("Sample.ID","Phylotype","Counts","Presence")
  print(ggplot2::ggplot(xyz, ggplot2::aes(y = Sample.ID, x = Phylotype, fill = log(Counts))) +
          ggplot2::geom_raster() + ggplot2::theme_bw() + ggplot2::ggtitle(title))
}


#' simulate data from a Dirichlet-Multinomial regression model
#'
#' @note Requires the dirmult and MASS packages
#'
#' @param n_obs: the number of samples
#' @param n_vars: number of covariates excluding the intercept
#' @param n_taxa: number of species
#' @param n_relevant_vars: number of relevant nutrients
#' @param n_relevant_taxa: number of relevant species
#' @param beta_min: minimum absolute value of the regression parameters
#' @param beta_max: maximum absolute value of the regression parameters
#' @param signoise: scalar multiplier on the regression parameters
#' @param n_reads_min: lower bound on uniform distribution for number of reads in each sample
#' @param n_reads_max: upper bound on uniform distribution for number of reads in each sample
#' @param theta0: the dispersion parameter
#' @param rho: the correlation between covariates
#'
#' @return XX: (design matrix) with intercept: n_obs * (n_vars + 1)
#' @return YY: (count matrix) rows: n_obs samples, columns: n_taxa species
#' @return alphas: simulated intercept vector
#' @return betas: simulated coefficient matrix n_taxa * (n_vars + 1)
#' @return n_reads_min, n_read_max: row sum parameters
#' @return theta0, phi, rho, signoise: simulation inputs
#'
#' @export
simulate_dirichlet_multinomial_regression = function(n_obs = 100,
                                                     n_vars = 100,
                                                     n_taxa = 40,
                                                     n_relevant_vars = 4,
                                                     n_relevant_taxa = 4,
                                                     beta_min = 0.5,
                                                     beta_max = 1.0,
                                                     signoise = 1.0,
                                                     n_reads_min = 1000,
                                                     n_reads_max = 2000,
                                                     theta0 = 0.01,
                                                     rho = 0.4){

  # check for required packages
  if(!require(dirmult)){
    stop("dirmult package required")
  }
  if(!require(MASS)){
    stop("MASS package required")
  }
  # covariance matrix for predictors
  Sigma = matrix(1, n_vars, n_vars)
  Sigma = rho^abs(row(Sigma) - col(Sigma))
  # include the intercept
  XX = cbind(rep(1, n_obs),
             scale(MASS::mvrnorm(n = n_obs, mu = rep(0, n_vars), Sigma = Sigma)))
  # empties
  YY = matrix(0, n_obs, n_taxa)
  betas = matrix(0, n_taxa, n_vars)
  phi = matrix(0, n_obs, n_taxa)
  # parameters with signs alternating
  st = 0
  low_side = beta_min
  high_side = beta_max
  if(n_relevant_taxa != 1){
    # warning if the lengths don't match
    coef = suppressWarnings(seq(low_side, high_side, len = n_relevant_taxa) * c(1, -1))
  }else{
    coef = (low_side + high_side) / 2
  }
  coef_g = rep(1.0, len = n_relevant_vars)
  for(ii in 1:n_relevant_vars){
    # overlap species
    betas[(st:(st + n_relevant_taxa - 1)) %% n_taxa + 1, 3 * ii - 2] = coef_g[ii] * sample(coef)[((ii - 1):(ii + n_relevant_taxa - 2)) %% n_relevant_taxa + 1]
    st = st + 1
  }
  # -2.3 and 2.3 so that the intercept varies over three orders of magnitude
  intercept = runif(n_taxa, -2.3, 2.3)
  Beta = cbind(intercept, signoise * betas)
  # row totals
  ct0 = sample(n_reads_min:n_reads_max, n_obs, rep = T)
  for(ii in 1:n_obs){
    thisrow = as.vector(exp(Beta %*% XX[ii, ]))
    phi[ii, ] = thisrow/sum(thisrow)
    YY[ii, ] = dirmult::simPop(J = 1, n = ct0[ii], pi = phi[ii, ], theta = theta0)$data[1, ]
  }

  return(list(XX = XX, YY = YY, alphas = intercept, betas = Beta,
              n_reads_min = n_reads_min, n_reads_max = n_reads_max,
              theta0 = theta0, phi = phi, rho = rho, signoise = signoise))

}

