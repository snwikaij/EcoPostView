#' Title Meta-Analysis over Regression coefficients (MAR)
#'
#' @param estimate The model estimates
#' @param stderr The standard error of the model estimates
#' @param parameter The model parameter b0 or b1
#' @param predictor The name of the predictor
#' @param link_function The link of the model
#' @param grouping The names of the groups
#' @param gradient_dependency A value that expresses the gradient dependency, mean, variance, min, max, cv, etc.
#' @param prior_mu Priors for the mean
#' @param prior_mu_se Priors for the standard error of the mean of the mean
#' @param prior_sd_max Maximum value for the variance of the prior (it is a uniform prior starting at 0).
#' @param prior_weights The weights attached to each prior (prior odds)
#' @param get_prior_only Argument when true a list is returned for adjustment
#' @param lower_trunc lower truncation
#' @param upper_trunc upper truncation
#' @param n_chain The number of chains
#' @param n_thin Thinning step
#' @param n_seed Seed
#' @param n_iter Number of iterations
#' @param n_burnin Length of the burn-in
#' @param Rhat_warn When higher returns an warning
#' @param Eff_warn When lower returns a warning
#'
mar <- function(estimate, stderr, parameter, predictor,
                             link_function, grouping, gradient_dependency,
                             prior_mu=0,
                             prior_mu_se=1,
                             prior_sd_max=5,
                             prior_weights=1,
                             get_prior_only = FALSE,
                             lower_trunc=-10e+6,
                             upper_trunc=10e+6,
                             n_chain = 2,
                             n_thin = 1,
                             n_seed = 666,
                             n_iter=10000,
                             n_burnin=1000,
                             Rhat_warn = 1.01,
                             Eff_warn = 1000){

  #Combine input to generate levels
  level <- paste(parameter, predictor, link_function, grouping, sep = "_")

  #Unique number of levels
  Ll <- length(unique(level))

  #If prior mu is given appoint to data
  if(is.numeric(prior_mu) && nrow(prior_mu) == 1 | is.numeric(prior_mu) && is.null(nrow(prior_mu))){
    prior_mu        <- rep(prior_mu, Ll)
  }else if(
    nrow(prior_mu) != Ll ){stop(paste0("The length of the priors for mu (n=", nrow(prior_mu), ") is not the same as the length of unique levels (n=", Ll, ")."))}

  #If prior se is given appoint to data
  if(is.numeric(prior_mu_se) && nrow(prior_mu_se) == 1 | is.numeric(prior_mu_se) && is.null(nrow(prior_mu_se))){
    prior_mu_se     <- rep(prior_mu_se, Ll)
  }else if(
    nrow(prior_mu_se) != Ll){stop(paste0("The length of the priors for se (n=", length(prior_mu_se),") is not the same as the length of unique levels (n=", Ll, ")."))}

  #If multiple priors and weights are given asses the number
  if(is.data.frame(prior_mu) && is.data.frame(prior_mu_se) &&
     any(apply(prior_mu,2,is.numeric)) && any(apply(prior_mu,2,is.numeric))  &&
     length(prior_weights) == ncol(prior_mu) && length(prior_weights) == ncol(prior_mu_se)){
    npw <- length(prior_weights)}else if(!is.data.frame(prior_mu) && !is.data.frame(prior_mu_se) &&
                                         is.numeric(prior_mu) && length(prior_weights) == 1 && is.numeric(prior_mu_se) && length(prior_weights) == 1){
      npw <- length(prior_weights)}else{stop(paste0("The length of the prior model weights (n=",length(prior_weights),") is not the same as the length of the different model priors provided."))}


  #a source of gradient dependency (e.g., CV, Mean or Variance)
  if(is.null(gradient_dependency)){gradient_dependency <- rep(0, length(estimate))}

  #Place all given data in a list
  mod_data <- list(est=as.numeric(estimate),
                   se=as.numeric(stderr),
                   N=length(estimate),
                   level=factor(level, levels = unique(level)[order(unique(level))]),
                   L=Ll,
                   Lc=as.numeric(table(factor(level, levels = unique(level)[order(unique(level))]))),
                   grad_expr=ifelse(is.na(gradient_dependency), 1,gradient_dependency),
                   Pm=prior_mu,
                   Pe=prior_mu_se,
                   Ps=prior_sd_max,
                   pw=prior_weights,
                   npw=npw,
                   l_trunc=lower_trunc,
                   u_trunc=upper_trunc)

  #If only the prior is needed return
  if(get_prior_only){return(data.frame(Levels=levels(mod_data$level),
                                       Prior_mu=mod_data$Pm,
                                       Prior_se=mod_data$Pe))}else if(get_prior_only == F){

                                         if(mod_data$npw>1){
                                           meta_analysis <- function(){

                                             ##likelihood
                                             for (i in 1:N){
                                               est[i]  ~ dnorm(mu2[i], tau2[i])
                                               mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                               mu1[i]  <-mu[level[i]] + adjust1[level[i]]*grad_expr[i] + adjust2[level[i]]*log(1/se[i])
                                               tau2[i] <-1/(se[i]^2)

                                               tvar[i] <- 1/tau2[i]+bvar[level[i]]
                                               I[i]    <- bvar[level[i]]/tvar[level[i]]}

                                             sigma_individual  ~ dunif(0, Ps)

                                             ##priors
                                             for(j in 1:L){
                                               adjust1[j]      ~ dnorm(0, 1/0.15^2)
                                               adjust2[j]      ~ dnorm(0, 1/0.15^2)
                                               for(k in 1:npw){
                                                 mu_M[j,k]    ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}

                                               d[j]            ~ dcat(pw[1:npw])

                                               mu[j]           <- inprod(mu_M[j, 1:npw], d[j] == 1:npw)

                                               tau1[j]         <- 1/sigma_individual^2

                                               ##variance between studies
                                               bvar[j]   <- 1 / tau1[j]}

                                             for(k in 1:L){I2[k] <- 1-sum(I[k])/Lc[L]}}}else{
                                               meta_analysis <- function(){

                                                 ##likelihood
                                                 for (i in 1:N){
                                                   est[i]  ~ dnorm(mu2[i], tau2[i])
                                                   mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                                   mu1[i]  <-mu[level[i]] + adjust1[level[i]]*grad_expr[i] + adjust2[level[i]]*log(1/se[i])
                                                   tau2[i] <-1/(se[i]^2)

                                                   tvar[i] <- 1/tau2[i]+bvar[level[i]]
                                                   I[i]    <- bvar[level[i]]/tvar[level[i]]}

                                                 sigma_individual  ~ dunif(0, Ps)

                                                 ##priors
                                                 for(j in 1:L){
                                                   adjust1[j]      ~ dnorm(0, 1/0.15^2)
                                                   adjust2[j]      ~ dnorm(0, 1/0.15^2)

                                                   mu[j]           ~ dnorm(Pm[j], 1/Pe[j]^2)



                                                   tau1[j]         <- 1/sigma_individual^2

                                                   ##variance between studies
                                                   bvar[j]   <- 1 / tau1[j]}

                                                 for(k in 1:L){I2[k] <- 1-sum(I[k])/Lc[L]}}}

                                         #Run the model
                                         model <- jags.parallel(data = mod_data,
                                                                model.file = meta_analysis,
                                                                parameters.to.save =  c("mu", "d", "adjust1", "adjust2", "I2"),
                                                                n.chains = n_chain,
                                                                n.thin = n_thin,
                                                                jags.seed = n_seed,
                                                                n.iter = n_iter,
                                                                n.burnin = n_burnin)

                                         #Warning messages for bad mixing chains
                                         if(any(c(model$BUGSoutput$summary[-1,8]>Rhat_warn, model$BUGSoutput$summary[-1,9]<Eff_warn))){warning("The Rhat or/and effective sample size for some >",Rhat_warn," or/and <",Eff_warn,". Chains might not be mixing.")}

                                         extract_chain <- function(chains, data){

                                           #Extract chains from the model
                                           mcmcdf             <- as.data.frame(chains)

                                           #Appoint names to the chains of mu
                                           colnames(mcmcdf)   <- levels(data$level)

                                           #Transform to long format
                                           mcmclong           <- tidyr::gather(mcmcdf)

                                           #Split names and appoint information to chains
                                           mcmclong           <- cbind.data.frame(do.call(rbind, strsplit(mcmclong$key, "_")), mcmclong[,2])
                                           colnames(mcmclong) <- c("parameter", "predictor", "link", "group", "estimate")
                                           return(mcmclong)}

                                         #Extract chains for means
                                         mcmc_mu         <- extract_chain(model$BUGSoutput$sims.list$mu, mod_data)

                                         #Extract chains for adjustment1
                                         mcmc_adj1        <- extract_chain(model$BUGSoutput$sims.list$adjust1, mod_data)

                                         #Extract chains for adjustment2
                                         mcmc_adj2        <- extract_chain(model$BUGSoutput$sims.list$adjust2, mod_data)

                                         #Extract chains for I2
                                         mcmc_I2         <- extract_chain(model$BUGSoutput$sims.list$I2, mod_data)

                                         #Extract chains for posterior weights
                                         if(mod_data$npw>1){
                                         mcmc_podd       <- extract_chain(model$BUGSoutput$sims.list$d, mod_data)

                                         #Total support per model
                                         support         <- table(mcmc_podd$estimate)/sum(table(mcmc_podd$estimate))/mod_data$pw}


                                         return(list(Estimates=split(mcmc_mu, mcmc_mu$parameter),
                                                     Gain=ifelse(mod_data$npw>1, split(mcmc_podd, mcmc_mu$parameter), NA),
                                                     N_level=table(mod_data$level),
                                                     Chains_mu=mcmc_mu,
                                                     Chains_adj1=mcmc_adj1,
                                                     Chains_adj2=mcmc_adj2,
                                                     Chains_podd=ifelse(mod_data$npw>1, mcmc_podd, NA),
                                                     Chains_I2=mcmc_I2,
                                                     Poster_to_prior_odds=ifelse(mod_data$npw>1, support, NA),
                                                     Prior_weight=mod_data$pw,
                                                     JAGS_model=model,
                                                     Data=mod_data,
                                                     Priors=data.frame(Levels=levels(mod_data$level),
                                                                       Prior_mu=mod_data$Pm,
                                                                       Prior_se=mod_data$Pe)))}}

