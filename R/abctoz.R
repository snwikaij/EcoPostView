#' Aproximated Bayesian Computation to derive Z (ABC-to-z)
#'
#' @param z A numeric vector of z values
#' @param nsim Number of simulations
#' @param prior_mu Prior values for mu(Z) following a truncated (0) normal distribution, which needs to notated c(mean, sd)
#' @param prior_sd Prior values for sd(Z) following a gamma distribution, which needs to be notated c(shape, rate)
#' @param prior_threshold Prior values for the threshold for 'significance'
#' @param prior_cens Prior values for the fractions of censored samples following a beta distribution, which needs to be notated as c(alpha, beta)
#' @param print_progress Print the progress of simulations
#' @param seed Set the seed
#'
abctoz <- function(z,
                   nsim=50000,
                   prior_mu=c(0, 1),
                   prior_sd=c(1, 1),
                   prior_threshold=c(1.96, 0.1),
                   prior_cens=c(1, 2),
                   distribution="t",
                   print_progress=T,
                   seed=123){

#number of samples and prior array
z                <- abs(z)
nz               <- length(z)

priors           <- as.data.frame(array(NA, dim=c(nsim, 7)))
colnames(priors) <- c("mu", "sd", "threshold", "cens", "n1", "n2", "dist")

priors[,1] <- truncnorm::rtruncnorm(nsim, a = 0, b=Inf, mean=prior_mu[1], sd=prior_mu[2])
priors[,2] <- rgamma(nsim, prior_sd[1], prior_sd[2])
priors[,3] <- truncnorm::rtruncnorm(nsim, a = 0, b=Inf, mean=prior_threshold[1], sd=prior_threshold[2])
priors[,4] <- rbeta(nsim, prior_cens[1], prior_cens[2])
priors[,5] <- round(nz*priors[,4])
priors[,6] <- nz-priors[,5]

mod_dens <- function(x){d <- density(x); d$x[which.max(d$y)]}

#parameters of the data
data_mu  <- mean(z)
data_mod <- mod_dens(z)
data_sd  <- sd(z)

#store simulations
sim_list     <- vector("list", nsim)

#Create a distance function
dist_fun <- function(sim){
sum(abs(c(mean(sim), mod_dens(sim), sd(sim))-c(data_mu, data_mod, data_sd)))/3}

#create some free space (is it at all needed?)
gc()

#set seed
set.seed(seed)

#ABC-rejection algorithm
for(i in 1:nrow(priors)){

if(print_progress == T) print(i)

if(distribution=="z"){
sim_list[[i]]         <- c(truncnorm::rtruncnorm(priors$n1[i], a=priors$threshold[i], b=Inf, mean=priors$mu[i], sd=priors$sd[i]),
                           truncnorm::rtruncnorm(priors$n2[i], a=0, b=Inf, mean=priors$mu[i], sd=priors$sd[i]))}
else if(distribution=="t"){
sim_list[[i]]         <- c(truncdist::rtrunc(priors$n1[i], a=priors$threshold[i], b=Inf, df=3, spec="t")*sqrt(priors$sd[i]^2 * (3-2)/3)+priors$mu[i],
                           truncdist::rtrunc(priors$n2[i], a=0, b=Inf, df=3, spec="t")*sqrt(priors$sd[i]^2 * (3-2)/3)+priors$mu[i])}

priors[i,7]      <- dist_fun(sim_list[[i]])}

return(list(iterations=priors, data=z, simulations=sim_list))}
