#' Analytical approximate Bayesian Meta Analysis
#'
#' @param estimate A vector containing the effect-size
#' @param stderr A vector containing the standard error
#' @param prior_mu Prior for the mean
#' @param prior_se Prior for the se
#' @param prior weights for prior odds (default=1/number of priors)
#' @param interval Credibility intervals for the summary (default=0.9)
#' @param RE An argument indicating if RE or FE should be used (default RE=TRUE)
#'
#' @description
#' Analytical Bayesian meta-analysis with random-effect (RE) or fixed-effect (FE). For the estimation for
#' tau^2 the DSL DerSimonian and Laird (1986). This approach is easy yet slightly underestimates the 'true'
#' variance (O Bai et al. 2016). That being said the between study variance is therefore treated as having an
#' improper uniform prior, which simplifies most calculations.
#'
#' @export
abmeta <- function(estimate, stderr, prior_mu=0, prior_se=1000, prior_weights=NULL, interval=0.9, RE=T) {

  #Give warning if estimate length is 1
  if(length(estimate) == 1){RE <- F; warning("Number of estimates is 1 then RE is automatically set to FALSE.")}

  #Stop if estimate an se are not of the same length
  if(length(estimate) != length(stderr)){stop("Vector of estimates is not of the same length as the vector of standard errors.")}

  #If length weights is not length of priors then set to average
  if(length(prior_mu) && !is.null(prior_weights) != length(prior_weights)){ if(is.null(prior_weights)){prior_weights <- rep(1/length(prior_mu), length(prior_mu))}}

  if(is.null(prior_weights)){prior_weights <- rep(1/length(prior_mu), length(prior_mu))}

  #Analytical posterior based on conjugation
  postvals <- function(data_mu, data_sigma, prior_mu, prior_sigma){

    post_mu <- (prior_mu/prior_sigma^2+sum(data_mu / data_sigma^2))/
      (1/prior_sigma^2+sum(1/data_sigma^2))

    post_sigma <- sqrt(1/((1/prior_sigma^2)+sum(1/data_sigma^2)))

    return(c(mu = post_mu, sigma = post_sigma))}

  #Loop over all priors
  posteriors <- sapply(1:length(prior_mu), function(k){
    postvals(estimate, stderr, prior_mu[k], prior_se[k])})

  #Extract posterior
  post_mu    <- posteriors["mu", ]
  post_se    <- posteriors["sigma", ]

  #Calculate the  BMA
  pooled     <- sum(prior_weights * post_mu)
  se         <- sqrt(sum(prior_weights*(post_se^2+post_mu^2))-pooled^2)

  #RE Model
  if (RE == T){
    #weights and pooled mu
    w      <- 1/stderr^2
    pooled <- sum(estimate*w)/sum(w)

    #DL method for tau2
    Q      <- sum(w*(estimate-pooled)^2)
    tau2   <- max(0, (Q-(length(estimate)-1))/(sum(w)-sum(w^2)/sum(w)))

    #Use tau2 for new posteriors
    posteriors <- sapply(1:length(prior_mu), function(k) {
      postvals(estimate, sqrt(stderr^2+tau2), prior_mu[k], prior_se[k])})

    #Extract posterior
    post_mu <- posteriors["mu", ]
    post_se <- posteriors["sigma", ]

    #Calculate the  BMA
    pooled <- sum(prior_weights*post_mu)
    se     <- sqrt(sum(prior_weights*(post_se^2+post_mu^2))-pooled^2)}

  #Simple ETI intervals because its is normal ETI is HDI
  ci <- pooled + se * c(-1, 1) * qnorm(interval + ((1 - interval) / 2))

  return(round(c(mu=as.numeric(pooled), se=as.numeric(se), ll=ci[1], ul=ci[2]), 4))}
