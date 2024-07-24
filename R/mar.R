#' Title
#'
#' @param estimate Parameter b0 or b1 estimated from a LM or GLM
#' @param stderr Standard error belonging to b0 and b1
#' @param parameter A category (name) either b0 or b1
#' @param predictor A predictor name (e.g., salinity) as multiple predictors can be handled.
#' @param link_function The link function of the LM or GLM currently only 'identity', 'logit' and 'log' are supported
#' @param grouping A category (name) for the group as multiple groups can be handled.
#' @param PEESE An argument indicating if PEESE should be used (default  PEESE=TRUE)
#' @param RE An argument indicating if RE or FE should be used (default RE=TRUE)
#' @param prior_mu Prior for the mean which can be vector or matrix
#' @param prior_mu_se Prior for the se which can be vector or matrix
#' @param prior_sigma_max Prior for sigma is not used when RE=FALSE and is a uniform prior starting at 0 restricted and given value (default=5)
#' @param prior_weights Prior weights when Bayesian model averaging is applied
#' @param interval The interval level for the summary (default=0.9)
#' @param get_prior_only If it is unclear how many levels and how to formulate multiple priors for each level this argument will only return a data frame of priors
#' @param n_chain Number of chains
#' @param n_thin Thinning interval of the chains
#' @param n_seed Seed
#' @param n_iter Number of iterations
#' @param n_burnin Burn-in period
#' @param Rhat_warn Warning leve for Rhat
#' @param Eff_warn Warning level for Effective sample size
#'
#' @return
#' Returns a summary, all chains per level and a list of all model specifications
#'
#' @description
#' Full Bayesian meta-analytic method using Bayesian Model Averaging, with random-effect (RE), fixed-effect (FE) in combination
#' with Precision Effect Estimate with Standard Error (PEESE; Stanley  and Doucouliagos (2013)). The goals of this function is to
#' analyse the parameter estimates intercepts (b0) and regression coefficients (b1) of LM or GLM models with 'identity', 'log' or 'logit'
#' link.
#'
#' This was package was especially develop to analyse the  parameters when the independent variable (predictor variable) was log (natural log)
#' transformed and either the link functions are 'log' or 'logit' are used, or the dependent variable (target or response variable) is log transformed and
#' 'identity' link is used. Under such conditions the estimate parameter is named the (semi-)elasticity coefficient. The
#' elasticity coefficient indicates the percentage change in the target variable per 1 percent increase in the predictor variable.
#' The advantage is that the units of the regression coefficient (b1) are retained: g(E(y|x))=b0+b1*log(x)+error. If the link function (g) is a
#' log-link the log(E(y|x))=b0+b1*log(x) and so b1=log(E(y|x))/log(x) the units for x (e.g.,TP mg/L) and y (e.g., species richness) are retained.
#' And b1 still has its interpretation on the log-scale (log(species richness)/log(TP mg/L)). Therefore if b1=-0.2 this indicates
#' a 0.2% decline in species richness per increase in 1% TP. Furthermore, the decline in species richness can still be predicted
#' as a function of the increase in log(x) because 'species  richness' = exp(b1*log(x)). This also makes broad comparisons
#' between different 'groupings' possible due to the scaling of the parameters.
#'
#' This does not mean that other estimated parameters cannot be analysed. However, I specifically designed this for comparisons
#' (e.g., using density plots in the pdplot function) among such parameters or to generate  and finally to generate HOP-lines
#' using the (hop function). This can all be performed without losing the interpretation of the parameter which then tells us something
#' ecological senseful between predictor and target variable. This 'sense' is lost when looking at correlation coefficients, or standardized
#' effect sizes (Tukey 1969; Baguley 2009; Correl et al., 2020). Hence, what does Cohens' D tell us about the impact of oxygen depletion on EPT-taxa
#' When it is 1 that the mean difference between two groups divided by its pooled standard deviation is 1. An how do we judge this to
#' be ecologically 'relevant'. Relatively little considering there is a gradual relation between O2 and the number of EPT-taxa, which is only
#' retained within the units of the parameter e.g., -0.2 EPT-taxa/log(O2 mg/L). Additionally is our data so noisy that any deviation is
#' expected due to the noise alone. More importantly we cannot either asses our long-run estimations or used them in our prior or compare
#' our posterior. Thus transforming to standardized-effect sizes is an increadible waste of information if the underlying data not is
#' provided.
#'
#'
mar <- function(estimate, stderr, parameter, predictor,
                link_function, grouping,
                PEESE=TRUE, RE=TRUE,
                prior_mu=0,
                prior_mu_se=0.5,
                prior_sigma_max=5,
                prior_weights=1,
                interval=0.9,
                get_prior_only = FALSE,
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

  #a random effect data frame
  if(!is.null(random) && is.data.frame(random)){
    if(nrow(random)!=length(estimate)){stop("Random  effect  needs to be a data frame with
                                           factors of same length as the estimates")}}

  #Place all given data in a list
  mod_data <- list(est=as.numeric(estimate),
                   se=as.numeric(stderr),
                   N=length(estimate),
                   level=factor(level, levels = unique(level)[order(unique(level))]),
                   L=Ll,
                   Lc=as.numeric(table(factor(level, levels = unique(level)[order(unique(level))]))),
                   Pm=prior_mu,
                   Pe=prior_mu_se,
                   Ps=prior_sigma_max,
                   pw=prior_weights,
                   npw=npw,
                   PEESE=ifelse(PEESE==TRUE, 1, 0))

  #If only the prior is needed return
  if(get_prior_only){return(data.frame(Levels=levels(mod_data$level),
                                       Prior_mu=mod_data$Pm,
                                       Prior_se=mod_data$Pe))}else if(get_prior_only == F){

                                       if(RE==TRUE){
                                         if(mod_data$npw>1){
                                           ##1.##RE model with npw>1
                                           meta_analysis <- function(){

                                             ##likelihood
                                             for (i in 1:N){
                                               est[i]  ~ dnorm(mu2[i], tau2[i])
                                               mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                               mu1[i]  <- mu[level[i]] + ifelse(PEESE==1, beta_PEESE[level[i]]*log(1/se[i]), 0)
                                               tau2[i] <- 1/(se[i]^2)}

                                             sigma  ~ dunif(0, Ps)

                                             ##priors
                                             for(j in 1:L){
                                               beta_PEESE[j]~ dnorm(0, 1/0.5^2)
                                               tau1[j]      <- 1/sigma^2

                                             for(k in 1:npw){
                                               mu_M[j,k]    ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}
                                               d[j]         ~ dcat(pw[1:npw])
                                               mu[j]        <- inprod(mu_M[j, 1:npw], d[j] == 1:npw)}

                                             #I2 per level
                                             for(l in 1:L){
                                               I2[l] <- (1/tau1[level[l]])/(1/tau2[l]+1/tau1[level[l]])}}}
                                         else{
                                           ##2.##RE model with npw=0
                                           meta_analysis <- function(){

                                           ##likelihood
                                           for (i in 1:N){
                                             est[i]  ~ dnorm(mu2[i], tau2[i])
                                             mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                             mu1[i]  <- mu[level[i]] + ifelse(PEESE==1, beta_PEESE[level[i]]*log(1/se[i]), 0)
                                             tau2[i] <- 1/(se[i]^2)}

                                           sigma  ~ dunif(0, Ps)

                                           ##priors
                                           for(j in 1:L){
                                             beta_PEESE[j]   ~ dnorm(0, 1/0.5^2)
                                             mu[j]           ~ dnorm(Pm[j], 1/Pe[j]^2)
                                             tau1[j]         <- 1/sigma^2}

                                           #I2 per level
                                           for(l in 1:L){
                                             I2[l] <- (1/tau1[level[l]])/(1/tau2[l]+1/tau1[level[l]])}}}}
                                       else{
                                         if(mod_data$npw>1){
                                           ##3.##FE model with npw>1
                                           meta_analysis <- function(){

                                                 ##likelihood
                                                 for (i in 1:N){
                                                   est[i]  ~  dnorm(mu1[i], tau1[i])
                                                   mu1[i]  <- mu[level[i]] +ifelse(PEESE==1,  beta_PEESE[level[i]]*log(1/se[i]^2), 0)
                                                   tau1[i] <- 1/(se[i]^2)}

                                                 ##priors
                                                 for(j in 1:L){
                                                   beta_PEESE[j]   ~ dnorm(0, 1/0.5^2)

                                                 for(k in 1:npw){
                                                   mu_M[j,k]       ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}
                                                   d[j]            ~ dcat(pw[1:npw])
                                                   mu[j]           <- inprod(mu_M[j, 1:npw], d[j] == 1:npw)}}}
                                         else{
                                           ##4.##FE model with npw=0
                                           meta_analysis <- function(){

                                                   ##likelihood
                                                   for (i in 1:N){
                                                     est[i]  ~  dnorm(mu1[i], tau1[i])
                                                     mu1[i]  <- mu[level[i]] + ifelse(PEESE==1, beta_PEESE[level[i]]*log(1/se[i]^2), 0)
                                                     tau1[i] <- 1/(se[i]^2)}

                                                   ##priors
                                                   for(j in 1:L){
                                                     beta_PEESE[j]   ~ dnorm(0, 1/0.5^2)
                                                     mu[j]           ~ dnorm(Pm[j], 1/Pe[j]^2)}}}}

                                         #Run the model
                                         model <- jags.parallel(data = mod_data,
                                                                model.file = meta_analysis,
                                                                parameters.to.save =  c("mu", "d", "I2", "sigma", "beta_PEESE"),
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
                                         mcmc_mu           <- extract_chain(model$BUGSoutput$sims.list$mu, mod_data)
                                         split_par         <- split(mcmc_mu, mcmc_mu$parameter)

                                         #summarize mu
                                         interval_level              <- 0.5+c(-1,1)*interval/2
                                         maxpost <- function(x){d <- density(x); d$x[which.max(d$y)]}
                                         basic_summary        <- aggregate(data=split_par$b1, estimate~., function(x) c(maxpost(x), mean(x), sd(x),
                                                                                                                        quantile(x, interval_level[1]), quantile(x, interval_level[2])))
                                         if(RE==T){
                                           mcmc_I2       <- extract_chain(model$BUGSoutput$sims.list$I2, mod_data)
                                           mcmc_I2_b1    <- mcmc_I2[mcmc_I2$parameter!="b0",]
                                           mcmc_sigma    <- model$BUGSoutput$sims.list$sigma

                                           i2_summary                    <- aggregate(data=mcmc_I2_b1, estimate~., mean)$estimate
                                           basic_summary                 <- cbind(basic_summary[c(1:4)], round(unlist(basic_summary$estimate),4), round(i2_summary, 4))
                                           colnames(basic_summary)[5:10] <- c("map", "mu", "se", "ll", "ul", "I2")}
                                         else{
                                           mcmc_I2 <- NULL
                                           mcmc_sigma <- NULL

                                           basic_summary                 <- cbind(basic_summary[c(1:4)], round(unlist(basic_summary$estimate),4))
                                           colnames(basic_summary)[5:9] <- c("map", "mu", "se", "ll", "ul")}

                                         #Extract chains for PEESE
                                         if(PEESE==T){
                                           mcmc_peese         <- extract_chain(model$BUGSoutput$sims.list$beta_PEESE, mod_data)}else{mcmc_peese <- NULL}

                                         #Extract chains for posterior weights
                                         if(mod_data$npw>1){
                                           mcmc_podd       <- extract_chain(model$BUGSoutput$sims.list$d, mod_data)}else{mcmc_podd <- NULL}

                                         return(list(Summary=basic_summary,
                                                     Estimates=split(mcmc_mu, mcmc_mu$parameter),
                                                     Chains_mu=mcmc_mu,
                                                     Chains_I2=mcmc_I2,
                                                     Chains_sigma=mcmc_sigma,
                                                     Chains_PEESE=mcmc_peese,
                                                     Chains_podd=mcmc_podd,
                                                     N_level=table(mod_data$level),
                                                     model=list(Prior_weight=mod_data$pw, JAGS_model=model,
                                                                Data=mod_data, Priors=data.frame(Levels=levels(mod_data$level),
                                                                                                 Prior_mu=mod_data$Pm,
                                                                                                 Prior_se=mod_data$Pe))))}}

