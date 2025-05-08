
dtest <- as.data.frame(array(rnorm(200), c(100,2)))
b0=5
b1=0.5
b2=0.9
sigma <- 0.3
y <- rpois(nrow(dtest), b0+b1*dtest[,1]+b2*dtest[,2])#, sigma)

dtest[3] <- as.integer(y)
dtest <- dtest[dtest$V3 != 0,]

pairs(dtest)

plot(dtest$V2, dtest$V3)

results <- glmmJAGS(V3~V2+V1, data=dtest, family = "negbinom_log")

hist(results$JAGS_output$BUGSoutput$sims.list$residuals)


glmmJAGS <- function(formula=NULL, random=NULL, data=NULL, family="norm_ident",
                     prior_mu=0, prior_mu_se=100){

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
    mu[i]   <- b0+inprod(bn[],x[i,])}

  b0 ~ dnorm(prior_mu, prior_mu_se)

  for(j in 1:nx){bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

  sigma ~ dgamma(0.01, 0.01)
  tau   <- 1/sigma^2}

  return_parameters <- c("b0", "bn", "sigma")

}else if(family=="norm_log"){

  model <- function(){

    for(i in 1:ni){
      y[i]    ~ dnorm(mu[i], tau)
      log(mu[i])   <- b0+inprod(bn[],x[i,])}

    b0 ~ dnorm(prior_mu, prior_mu_se)

    for(j in 1:nx){bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

    sigma ~ dgamma(0.01, 0.01)
    tau   <- 1/sigma^2}

  return_parameters <- c("b0", "bn", "sigma")

}else if(family=="gamma_identity"){

  model <- function(){

    for(i in 1:ni){
      y[i]    ~ dgamma(mu[i]^2/sigma^2, mu[i]/sigma^2)
      mu[i]   <- b0+inprod(bn[],x[i,])}

    b0 ~ dnorm(prior_mu, prior_mu_se)

    for(j in 1:nx){bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

    sigma ~ dgamma(0.01, 0.01)}

  return_parameters <- c("b0", "bn", "sigma")

}else if(family=="gamma_log"){

  model <- function(){

    for(i in 1:ni){
      y[i]         ~ dgamma(mu[i]^2/sigma^2, mu[i]/sigma^2)
      log(mu[i])   <- b0+inprod(bn[],x[i,])}

    b0 ~ dnorm(prior_mu, prior_mu_se)

    for(j in 1:nx){bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

    sigma ~ dgamma(0.01, 0.01)}

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

    size ~ dgamma(0.01, 0.01)}

  return_parameters <- c("b0", "bn", "size")

}else if(family=="negbinom_log"){

  model <- function(){

    for(i in 1:ni){
      y[i]         ~ dnegbin(size/(size+mu[i]), size)

      log(mu[i])   <- b0+inprod(bn[],x[i,])}

    b0 ~ dnorm(prior_mu, prior_mu_se)

    for(j in 1:nx){bn[j] ~ dnorm(prior_mu, 1/prior_mu_se^2)}

    size ~ dgamma(0.01, 0.01)}

  return_parameters <- c("b0", "bn", "size")

}else{stop("family not recognized")}

parts <- unlist(strsplit(argument[3], "\\+"))
x     <- as.matrix(data[colnames(data) %in% parts])
y     <- data[,argument[2]]

model_list  <- list(y=y, x=x, ni=nrow(x), nx=ncol(x), prior_mu=prior_mu, prior_mu_se=prior_mu_se)

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
