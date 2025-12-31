#' Empirical Bayesian Meta Analysis
#'
#' @param estimate A vector containing the effect-size
#' @param stderr A vector containing the standard error
#' @param prior_mu Prior for the mean
#' @param prior_mu_se Prior for the se
#' @param prior_weights Weights for prior odds (default=1/number of priors)
#' @param tau_2 The method used to estimate tau^2 either a heuristic method 'HE0' or
#' 'DSL'=DerSimonian and Laird method (default = DSL).
#' @param interval Credibility intervals for the summary (default=0.9)
#' @param RE An argument indicating if RE or FE should be used (default RE=TRUE)
#'
#' @description
#' Empirical Bayesian meta-analysis with random-effect (RE) or fixed-effect (FE). For the estimation of
#' tau^2 an Empirical Bayesian method was used utilizing DerSimonian and Laird (1986) or DSL estimator.
#' This approach is easy but intervals are slightly smaller (O Bai et al. 2016). That being said the between study
#' variance is therefore treated as having an improper uniform prior, which simplifies most calculations.
#'
#' The heuristic methods was (naive) attempt to formulate an approximation to a closed form solution. I started with a
#' simple method of moment estimator but did not progress further. The HE is still rather empirical and results in similar results to the DSL and REML approach.
#' \eqn{\tau^2 = \max\left(0, \frac{1}{n} \sum_{i=1}^{n} (\hat{\theta_i} - \theta_{\text{pooled}})^2 - \frac{\sum_{i=1}^{n} (w_i \cdot se_i)}{\sum_{i=1}^{n} w_i} \right)}
#' However, HE will create wider intervals over the estimate under smaller sample sizes. Yet, by default DSL is still
#' preferred.
#
#' @export
abmeta <- function(estimate, stderr, prior_mu=0, prior_mu_se=1000, prior_weights=NULL,
                   a=-Inf, b=Inf, tau_2="DSL", interval=0.9, RE=T, warnings=F) {

  #Give warning if estimate length is 1
  if(length(estimate) == 1){RE <- F;   if(warnings==T){warning("Number of estimates is 1 then RE is automatically set to FALSE.")}}

  #Stop if estimate an se are not of the same length
  if(length(estimate) != length(stderr)){stop("Vector of estimates is not of the same length as the vector of standard errors.")}

  #If length weights is not length of priors then set to average
  if(length(prior_mu) && !is.null(prior_weights) != length(prior_weights)){
    if(is.null(prior_weights)) {prior_weights <- rep(1/length(prior_mu), length(prior_mu))}}

  #Assign equal weights to the priors 1/2
  if(is.null(prior_weights)){prior_weights <- rep(1/length(prior_mu), length(prior_mu))}

  #For a and b
  if(length(a)>1 || length(b)>1){stop("Vector of truncation limits can only have a length of 1.")}

  #Analytical posterior based on conjugation
  postvals <- function(data_mu, data_sigma, prior_mu, prior_sigma, a=-Inf, b=Inf){

    #Mean posterior
    post_mu <- (prior_mu/prior_sigma^2 + sum(data_mu/data_sigma^2))/
      (1/prior_sigma^2 + sum(1/data_sigma^2))

    #Sd posterior
    post_sigma <- sqrt(1/(1/prior_sigma^2+sum(1/data_sigma^2)))

    #Truncation
    if(!is.infinite(a) || !is.infinite(b)) {
      alpha <- (a-post_mu)/post_sigma
      beta  <- (b-post_mu)/post_sigma
      Z     <- pnorm(beta)-pnorm(alpha)

      if(Z <= .Machine$double.eps) {
        stop("Truncated posterior has negligible denstity within a and b.")
      } else {
        phi_a <- dnorm(alpha)
        phi_b <- dnorm(beta)

        mu_trunc  <- post_mu+post_sigma*(phi_a-phi_b)/Z
        var_trunc <- post_sigma^2 * (1 +
                                       if(is.finite(alpha)) alpha * phi_a / Z else 0 -
                                       if(is.finite(beta))  beta * phi_b / Z else 0 -
                                       ((phi_a - phi_b)/Z)^2)
        post_mu    <- mu_trunc
        post_sigma <- sqrt(var_trunc)
      }
    }

    return(c(mu = post_mu, sigma = post_sigma))
  }

  #Loop over all priors
  posteriors <- sapply(1:length(prior_mu), function(k){
    postvals(estimate, stderr, prior_mu[k], prior_mu_se[k], a, b)})

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

    if (tau_2 == "HE") {
      # Heuristic method for tau^2
      tau2 <- max(0, (1/length(estimate))*sum((estimate-pooled)^2)-(sum(w*stderr)/sum(w)))
    } else if (tau_2 == "DSL") {
      #DSL method for tau^2
      Q    <- sum(w*(estimate-pooled)^2)
      tau2 <- max(0, (Q-(length(estimate)-1))/(sum(w)-sum(w^2)/sum(w)))
    } else {
      stop("Not a correct method, either DSL or HE.")
    }

    #Use tau^2 for new posteriors
    posteriors <- sapply(1:length(prior_mu), function(k) {
      postvals(estimate, sqrt(stderr^2+tau2), prior_mu[k], prior_mu_se[k])})

    #Extract posterior
    post_mu <- posteriors["mu", ]
    post_se <- posteriors["sigma", ]

    #Calculate the  BMA
    pooled <- sum(prior_weights*post_mu)

    #Calculate SE
    se  <- sqrt(sum(prior_weights*(post_se^2+post_mu^2))-pooled^2)

    sdpi<- sqrt(se^2+tau2)
  }

  #Simple ETI intervals because its is normal ETI is HDI
  ci    <- pooled+se*c(-1, 1)*qnorm(interval+((1-interval)/2))

  if(RE==T){
      results <- c(mu=as.numeric(pooled), se=as.numeric(se), ll=ci[1], ul=ci[2],  sd=as.numeric(sdpi), tau2=tau2, I2=tau2/(se^2+tau2), a=a, b=b, n=length(estimate))
  }else{
      results <- c(mu=as.numeric(pooled), se=as.numeric(se), ll=ci[1], ul=ci[2], sd=as.numeric(se), tau2=0, I2=0, a=a, b=b, n=length(estimate))
  }

  return(round(results, 4))}
