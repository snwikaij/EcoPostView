#' Analytical Bayesian Meta Analysis
#'
#' @param estimate A vector containing the effect-size
#' @param stderr A vector containing the standard error
#' @param mu_prior Prior for the mean
#' @param se_prior Prior for the se
#' @param interval Credibility intervals for the summary (default=0.9)
#' @param RE An argument indicating if RE or FE should be used (default RE=TRUE)
#'
#' @description
#' Analytical Bayesian meta-analysis with random-effect (RE) or fixed-effect (FE).
#'
#' @export
abmeta <- function(estimate, stderr, mu_prior=0, se_prior=1, interval=0.9, RE=T){

  postvals <- function(mu_data, sigma_data, mu_prior, sigma_prior){

    post_mu <- (mu_prior / sigma_prior^2 + mu_data /sigma_data^2) / (1 / sigma_prior^2 + 1 / sigma_data^2)
    post_sigma <- sqrt(1/(1 / sigma_prior^2 + 1 / sigma_data^2))

    c(mu=post_mu, sigma=post_sigma)}

  post_set <- postvals(estimate, stderr, mu_prior=mu_prior, sigma_prior=se_prior)

  theta <- post_set[1:length(estimate)]
  w     <- 1/post_set[c((length(estimate)+1):(length(estimate)*2))]^2

  pooled <- sum(theta*w)/sum(w)
  se     <- 1/sqrt(sum(w))

  if(RE==T){
    Q <- sum(w*(estimate-pooled)^2)
    tau2 <- max(0, (Q-(length(estimate)-1))/(sum(w)-sum(w^2)/sum(w)))
    w2 <- 1/(1/w+tau2)

    pooled <- sum(theta*w2)/sum(w2)
    se     <- 1/sqrt(sum(w2))}

  ci <- pooled+se*c(-1,1)*qnorm(interval+((1-interval)/2))

  round(c(mu=pooled, se=se, ll=ci[1], ul=ci[2]),4)}
