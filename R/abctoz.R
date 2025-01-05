#' Aproximated Bayesian Computation to derive Z (ABC-to-z)
#'
#' @param p A numeric vector of p-values
#' @param operator A vector that contains either <, > or NA (exact p-value)
#' @param nsim Number of simulations
#' @param prior_mu Prior values for mu(Z) following a truncated (0) normal distribution, which needs to notated c(mean, sd)
#' @param prior_sd Prior values for sd(Z) following a gamma distribution, which needs to be notated c(shape, rate)
#' @param prior_threshold Prior values for the threshold for 'significance'
#' @param prior_cens Prior values for the fractions of censored samples following a beta distribution, which needs to be notated as c(alpha, beta)
#' @param prior_uncens Prior values for the fractions of uncensored samples following a beta distribution, which needs to be notated as c(alpha, beta)
#' @param prior_odds_cens The likelihood of the prior odds
#' @param distribution Either z-distribution or shifted t-distribution with 1 df
#' @param density_steps An argument that determines the number of steps of the density function to derive the fit of the simulated to the data
#' @param seed Set the seed
#'
#' @importFrom foreach foreach
#' @importFrom foreach %dopar%
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#'
#' @export
abctoz <- function(p, operator=NULL,
                   nsim=100000,
                   prior_mu=c(1, 1),
                   prior_sd=c(1, 2),
                   prior_threshold=c(1.96, 0.1),
                   prior_cens=c(1, 1.5, 0.02, 1),
                   prior_uncens=c(1, 1.5, 0, 0.02),
                   prior_odds_cens=0.97,
                   distribution="z",
                   density_steps=100,
                   min_dens_steps=20,
                   seed=1){

  #number of samples and prior array
  nz               <- length(p)

  #generate a data frame for the priors
  priors           <- as.data.frame(array(NA, dim=c(nsim, 7)))
  colnames(priors) <- c("mu", "sd", "threshold", "cens", "n1", "n2", "H0/H1")

  #based on the prior odds generate priors for the censoring
  selection <- rbinom(nsim, 1, prior_odds_cens)
  prior_cens_values <- mapply(function(x) {
    if(x == 1){
      return(truncdist::rtrunc(1, shape1=prior_cens[1], shape2=prior_cens[2], a=prior_cens[3], b=prior_cens[4], spec="beta"))
    }else{
      return(truncdist::rtrunc(1, shape1=prior_uncens[1], shape2=prior_uncens[2], a=prior_uncens[3], b=prior_uncens[4], spec="beta"))}
  }, selection)

  #generate all priors and place in df
  priors[,1] <- truncnorm::rtruncnorm(nsim, a = 0, b=Inf, mean=prior_mu[1], sd=prior_mu[2])
  priors[,2] <- rgamma(nsim, prior_sd[1], prior_sd[2])
  priors[,3] <- truncnorm::rtruncnorm(nsim, a = 0, b=Inf, mean=prior_threshold[1], sd=prior_threshold[2])
  priors[,4] <- prior_cens_values
  priors[,5] <- round(nz*priors[,4])
  priors[,6] <- nz-priors[,5]
  priors[,7] <- selection

  #if no operator is given
  if(is.null(operator)){
    z        <- qnorm(1 - p/2)

    #generate a density curve for the data
    xlen   <- base::seq(0, max(z), length.out=density_steps)
    dval   <- approx(density(z)$x, density(z)$y, xout = xlen)$y
  }else{
    z_list <- vector("list", nsim)}

  #p to z imputation function
  catp <- function(x, p) {
    if (is.na(x)) {
      return(qnorm(1 - p/2))
    } else if (x == "<") { p <- runif(1, 0, p)
    } else if (x == ">"){ p <- runif(1, p, 1)}
    return(qnorm(1 - p/2))}

  #store simulations
  sim_list     <- vector("list", nsim)

  #create some free space
  gc()

  #set cores to use (minus one for safety)
  n_cores <- parallel::detectCores() - 1
  cl      <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  #set seed
  set.seed(seed)
  results <- foreach::foreach(i = 1:nsim, .packages = c("truncnorm", "truncdist")) %dopar% {

    #impute pvalues if half reported
    if(!is.null(operator)){
      z           <- mapply(catp, operator,  p)
      xlen        <- base::seq(0, max(z), length.out=density_steps)
      dval        <- approx(density(z)$x, density(z)$y, xout = xlen)$y
      imputed_z   <- z}else{imputed_z <- NULL}

    #normal and t dist
    if(distribution=="z"){
      sim_vals         <- c(truncnorm::rtruncnorm(priors$n1[i], a=priors$threshold[i], b=Inf, mean=priors$mu[i], sd=priors$sd[i]),
                            truncnorm::rtruncnorm(priors$n2[i], a=0, b=Inf, mean=priors$mu[i], sd=priors$sd[i]))}
    else if(distribution=="t"){
      sim_vals         <- c(truncdist::rtrunc(priors$n1[i], a=priors$threshold[i], b=Inf, df=3, spec="t")*sqrt(priors$sd[i]^2 * (3-2)/3)+priors$mu[i],
                            truncdist::rtrunc(priors$n2[i], a=0, b=Inf, df=3, spec="t")*sqrt(priors$sd[i]^2 * (3-2)/3)+priors$mu[i])}

    #fit density distribution
    sim_density <- approx(density(sim_vals)$x, density(sim_vals)$y, xout = xlen)$y
    n_dens      <- length(xlen)-sum(is.na(sim_density))

    #calulcate distance
    dsim        <- approx(density(sim_vals)$x, density(sim_vals)$y, xout = xlen)$y
    if(n_dens<20){dist <- Inf}else{dist <- sum(abs(dval-dsim), na.rm = T)/n_dens}

    list(simulations=sim_vals, distance=dist, imputed_z=imputed_z)}

  sim_list           <- lapply(results, `[[`, "simulations")
  simulations        <- lapply(results, `[[`, "imputed_z")
  priors$distance    <- sapply(results, `[[`, "distance")

  if(!is.null(operator)){zvals <- simulations}else{zvals <- z}

  return(invisible(list(iterations=priors, data=zvals, raw_data=qnorm(1 - p/2), simulations=sim_list)))}
