#' Simple Generalized Linear Models with JAGS
#'
#' @param formula Satandard notation y~x
#' @param data A dataset
#' @param family A set of "norm_ident", "norm_log", "gamma_ident", "gamma_log", "pois_ident", "pois_log", "negbinom_ident"
#' "negbinom_log"
#' @param prior_mu The prior for the mean of
#' @param prior_mu_se The prior for the standard error of the parameter
#' @param prior_dispersion The prior for the dispersion parameters (sigma or size) set via the gamma distribution
#'
#' @examples
#'df  <- as.data.frame(array(rnorm(200), c(100,2)))
#'b0=2
#'b1=0.25
#'b2=0.1
#'
#'df$V3 <- rpois(nrow(df), exp(b0+b1*df$V1+b2*df$V2))
#'
#'results <- glmmJAGS(V3~V1+V2, data=df, family = "pois_log")
#'
#' @export
glmJAGS <- function(formula=NULL, data=NULL, family="norm_ident",
                     prior_mu=0, prior_mu_se=100, prior_dispersion=c(0.001, 0.001)){

  argument.call <- match.call()
  argument      <- gsub(" ", "", formula)

if(is.null(formula)){stop("No formula provided")}
if(is.null(data)){stop("No data provided")}
if(argument[1] != "~"){stop("Not correct formula provided")}
if(length(argument) != 3){stop("Not correct formula provided")}

if(family=="norm_ident"){

  model <- function(){

  for(i in 1:ni){
    y[i]    ~ dnorm(mu[i], tau)

    mu[i]   <- inprod(bn[],x[i,])}

  for(j in 1:nx){bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

  sigma ~ dgamma(prior_dispersion[1], prior_dispersion[2])
  tau   <- 1/sigma^2}

  return_parameters <- c("b0", "bn", "sigma")

}else if(family=="norm_log"){

  model <- function(){

    for(i in 1:ni){
      y[i]    ~ dnorm(mu[i], tau)

      mu[i]   <- inprod(bn[],x[i,])}

    for(j in 1:nx) {bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

    sigma ~ dgamma(prior_dispersion[1], prior_dispersion[2])
    tau   <- 1/sigma^2}

  return_parameters <- c("b0", "bn", "sigma")

}else if(family=="gamma_identity"){

  model <- function(){

    for(i in 1:ni){
      y[i]    ~ dgamma(mu[i]^2/sigma^2, mu[i]/sigma^2)
      mu[i]   <- b0+inprod(bn[],x[i,])}

    b0 ~ dnorm(prior_mu, prior_mu_se)

    for(j in 1:nx){bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

    sigma ~ dgamma(prior_dispersion[1], prior_dispersion[2])}

  return_parameters <- c("b0", "bn", "sigma")

}else if(family=="gamma_log"){

  model <- function(){

    for(i in 1:ni){
      y[i]         ~ dgamma(mu[i]^2/sigma^2, mu[i]/sigma^2)
      log(mu[i])   <- b0+inprod(bn[],x[i,])}

    b0 ~ dnorm(prior_mu, prior_mu_se)

    for(j in 1:nx){bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

    sigma ~ dgamma(prior_dispersion[1], prior_dispersion[2])}

  return_parameters <- c("b0", "bn", "sigma")

}else if(family=="pois_ident"){

  model <- function(){

    for(i in 1:ni){
      y[i]         ~ dpois(mu[i])
      mu[i]       <- b0+inprod(bn[],x[i,])}

    b0 ~ dnorm(prior_mu, prior_mu_se)

    for(j in 1:nx){bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

  sigma <- 1}

  return_parameters <- c("b0", "bn", "sigma")

}else if(family=="pois_log"){

  model <- function(){

    for(i in 1:ni){
      y[i]         ~ dpois(mu[i])
      log(mu[i])  <- b0+inprod(bn[],x[i,])}

    b0 ~ dnorm(prior_mu, prior_mu_se)

    for(j in 1:nx){bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

    sigma <- 1}

  return_parameters <- c("b0", "bn", "sigma")

}else if(family=="negbinom_ident"){

  model <- function(){

    for(i in 1:ni){
      y[i]         ~ dnegbin(size/(size+mu[i]), size)

      mu[i]       <- b0+inprod(bn[],x[i,])}

    b0 ~ dnorm(prior_mu, prior_mu_se)

    for(j in 1:nx){bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

    size ~ dgamma(prior_dispersion[1], prior_dispersion[2])}

  return_parameters <- c("b0", "bn", "size")

}else if(family=="negbinom_log"){

  model <- function(){

    for(i in 1:ni){
      y[i]         ~ dnegbin(size/(size+mu[i]), size)

      log(mu[i])   <- b0+inprod(bn[],x[i,])}

    b0 ~ dnorm(prior_mu, prior_mu_se)

    for(j in 1:nx){bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

    size ~ dgamma(prior_dispersion[1], prior_dispersion[2])}

  return_parameters <- c("b0", "bn", "size")

}else{stop("family not recognized")}

############################
#standard formula y~x+...+x#
############################

parts <- unlist(strsplit(argument[3], "\\+"))
x     <- as.matrix(data[colnames(data) %in% parts])
y     <- data[,argument[2]]

############################
#model data in final format#
############################

model_list  <- list(y=y, x=x, ni=nrow(x), nx=ncol(x), rx=rx,
                      prior_mu=prior_mu, prior_mu_se=prior_mu_se, prior_dispersion=prior_dispersion)

output <- jags.parallel(data = model_list, model.file = model,
              parameters.to.save = return_parameters,
              n.iter = 5500,
              n.thin = 20,
              n.burnin = 500,
              n.chains = 10)

summary <- output$BUGSoutput$summary
rownames(summary)[1:c(1+model_list$nx)] <- c("Intercept", parts)

cat("Call:\n")
print(argument.call)
print(summary)

list(summary=summary, JAGS_output=output, data=data)}
