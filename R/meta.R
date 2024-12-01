#' Bayesian meta analysis function
#'
#' @param estimate Parameter b0 or b1 estimated from a LM or GLM
#' @param stderr Standard error belonging to b0 and b1
#' @param parameter A category (name) either b0 or b1
#' @param predictor  A predictor name (e.g., salinity) as multiple predictors can be handled
#' @param link_function The link function of the LM or GLM currently only 'identity', 'logit' and 'log' are supported and will only be used in the hop function
#' @param grouping A category (name) for the group as multiple groups can be handled
#' @param random An argument that needs to be a vector or  matrix of factors of the same number of rows as the estimate
#' @param method Indicates which adjustment performed 0 (='none'), 1 (='egger') or 2 (='peters')
#' @param RE An argument indicating if RE or FE should be used (default RE=TRUE)
#' @param Nsamp A vector with the number of samples for each estimate (only used when method = 2)
#' @param prior_mu Prior for the se which can be vector (Bayesian meta-analysis) or matrix (Bayesian meta-analysis with model averaging)
#' @param prior_mu_se Prior for the se which can be vector or matrix
#' @param prior_sigma_max Prior for sigma is and is a uniform prior starting at 0 restricted and given value (default=5) not used when RE=FALSE
#' @param interval Credibility intervals for the summary (default=0.9)
#' @param get_prior_only If it is unclear how many levels and how to formulate multiple priors for each level this argument will only return a data frame of priors so that formulating priors for levels is easier
#' @param n_chain Number of chains
#' @param n_thin Thinning interval of the chains
#' @param n_seed Seed
#' @param n_iter Number of iterations
#' @param n_burnin Burn-in period
#' @param Rhat_warn Warning level for Rhat
#' @param Eff_warn Warning level for Effective sample size
#' @param print_summary If TRUE it prints a summary
#'
#' @description
#' Full Bayesian meta-analytic method using Bayesian Model Averaging, with random-effect (RE), fixed-effect (FE) in combination
#' with Egger adjustment using inverse of the standard error (1/SE) (Moreno et al., 2009; Stanley  and Doucouliagos, 2013) or Peters adjustment (Moreno et al., 2009) using
#' the inverse of the sample size (1/N). The full use of the meta function is the combination with pdplot and hop function. Hence it  is used to analyse
#' the parameter estimates intercepts (b0) and regression coefficients (b1) on of LM or GLM models applied when the independent
#' variable is continues.
#'
#' This package was especially develop to analyse the  parameters when the independent variable (predictor variable) was log (natural log)
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
#' ecological senseful about the relation between predictor and target variable. This 'sense' is lost when looking at correlation coefficients, or standardized
#' effect sizes (Tukey 1969; Baguley 2009; Correl et al., 2020). Hence, what does Cohens' D tell us about the impact of oxygen depletion on EPT-taxa
#' When it is found to be 1: The mean difference between two groups divided by its pooled standard deviation. And how do we judge this to
#' be ecologically 'relevant'. Relatively little, considering there is a gradual relation between O2 and the number of EPT-taxa, which is only
#' retained within the units of the parameter e.g., -0.2 EPT-taxa/log(O2 mg/L). Additionally is our data so noisy that any deviation is
#' expected due to the noise alone. More importantly we cannot either asses our long-run estimations or used them in our prior or compare
#' our posterior. Thus transforming to standardized-effect sizes is an incredible waste of information if the underlying data not is
#' provided.
#'
#'@note
#'If Nsamp (NA) this is replaced by 1 since 1/N where N=1 is 1=1/1. If parameters b0 or b1 or link functions are not given still a meta analysis is performed. One
#'still needs to fill in the columns.
#'
#'@examples
#'data("example1")
#'
#'#Standard random effect (RE) meta-analysis
#'mod <-meta(estimate=example1$est, stderr=example1$se,
#'           parameter=example1$parameter, predictor=example1$predictor,
#'           link_function=example1$link, grouping=example1$group)
#'
#'#Standard random effect (RE) meta-analysis with egger's correction
#'mod <-meta(estimate=example1$est, stderr=example1$se,
#'           parameter=example1$parameter, predictor=example1$predictor,
#'           link_function=example1$link, grouping=example1$group, method=1)
#'
#'#Standard random effect (RE) meta-analysis with peter's correction
#'mod <-meta(estimate=example1$est, stderr=example1$se,
#'           parameter=example1$parameter, predictor=example1$predictor,
#'           link_function=example1$link, grouping=example1$group,
#'           Nsamp=example1$n, method=2)
#'
#' @importFrom stats sd
#' @importFrom stats aggregate
#' @importFrom stats density
#'
#' @export
meta <- function(estimate, stderr, parameter=NULL, predictor=NULL,
                link_function=NULL, grouping=NULL, random=NULL,
                method=0, RE=TRUE, Nsamp=NULL,
                prior_mu=0,
                prior_mu_se=10,
                prior_sigma_max=5,
                interval=0.9,
                get_prior_only = FALSE,
                n_chain = 2,
                n_thin = 1,
                n_seed = 666,
                n_iter=10000,
                n_burnin=1000,
                Rhat_warn = 1.01,
                Eff_warn = 1000,
                print_summary=FALSE){

  #If parameter is not given set to b1
  if(is.null(parameter)){parameter <- rep("b1", length(estimate))}

  #If predictor is not given set to none specified
  if(is.null(predictor)){predictor <- rep("none-specified", length(estimate))}

  #If link function is not given set to identity
  if(is.null(link_function)){link_function <- rep("identity", length(estimate))}

  #If grouping is not given set to non-specified
  if(is.null(grouping)){grouping <- rep("none-specified", length(estimate))}

  #Combine input to generate levels
  level <- paste(parameter, predictor, link_function, grouping, sep = "_")

  #Unique number of levels
  Ll <- length(unique(level))

  #If multiple priors and weights are given asses the number
  if(is.null(nrow(prior_mu)) && length(prior_mu) == 1 &&
     is.null(nrow(prior_mu_se)) && length(prior_mu_se) == 1){
    npw <- 1}else if(is.data.frame(prior_mu) && is.data.frame(prior_mu_se) &&
                     ncol(prior_mu) == ncol(prior_mu_se)){npw <- ncol(prior_mu)
    }else{stop(paste0("The the priors are not nummeric or the priors for prior_mu and prior_mu_se are not of the same dimensions."))}

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

  #a random effect data frame
  if(!is.null(random) && is.data.frame(random)){
    if(nrow(random)!=length(estimate)){stop("Random  effect  needs to be a data frame or matrix with
                                           factors of same length as the estimates")}}

  if(is.null(random)){RL <- 0; random <- NULL}else if(is.null(dim(random))){
  RL <- 1; random <- as.factor(random)}else{
  RL <- ncol(random); random[] <- lapply(random, as.factor); random  <- do.call(cbind, random)}

  if(all(method != c(0, 1, 2))){
    stop("method needs to be either 0 (='none'), 1 (='egger'), 2 (='peters')")}
  if(method == 2){
    if(is.null(Nsamp) | length(Nsamp) != length(estimate)){stop("If method is 'peters' then length of Nsamp needs to be of the same length as
                                             the estimates")}
    Nsamp[is.na(Nsamp)] <- 1}else{Nsamp <- rep(1, length(estimate))}

  #Place all given data in a list
  mod_data <- list(est=as.numeric(estimate),
                   se=as.numeric(stderr),
                   N=length(estimate),
                   level=factor(level, levels = unique(level)[order(unique(level))]),
                   L=Ll,
                   Nsamp=Nsamp,
                   Pm=prior_mu,
                   Pe=prior_mu_se,
                   Ps=prior_sigma_max,
                   npw=npw,
                   alpha_pw=rep(1, npw),
                   random=random,
                   R=RL,
                   method=method)

  #If only the prior is needed return
  if(get_prior_only){return(data.frame(Levels=levels(mod_data$level),
                                       Prior_mu=mod_data$Pm,
                                       Prior_se=mod_data$Pe))}else if(get_prior_only == F){

                                      if(RE==TRUE){
                                        if(mod_data$npw>1){
                                          ##1.##RE model with npw>1 and
                                          ##1.1 no random
                                          if(mod_data$R==0){
                                          meta_analysis <- function(){

                                            ##likelihood
                                            for (i in 1:N){
                                              est[i]  ~ dnorm(mu2[i], tau2[i])
                                              mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                              mu1[i]  <- mu[level[i]] + ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))
                                              tau2[i] <- 1/(se[i]^2)}

                                            resid <- est-mu2

                                            sigma  ~ dunif(0, Ps)

                                            ##Stochastic behavior for the weights
                                            pw[1:npw] ~ ddirch(alpha_pw[1:npw])

                                            ##priors
                                            for(j in 1:L){
                                              beta_adjust[j]~ dnorm(0, 1/10^2)
                                              tau1[j]      <- 1/sigma^2

                                              for(k in 1:npw){
                                                mu_M[j,k]    ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}
                                                d[j]         ~ dcat(pw[1:npw])
                                                mu[j]        <- inprod(mu_M[j, 1:npw], d[j] == 1:npw)}

                                            #I2 per level
                                            for(l in 1:L){
                                              I2[l] <- (1/tau1[level[l]])/(1/tau2[l]+1/tau1[level[l]])}}}
                                          ##1.2 only 1 random
                                          else if(mod_data$R==1){
                                          meta_analysis <- function(){

                                              ##likelihood
                                              for (i in 1:N){
                                                est[i]  ~ dnorm(mu2[i], tau2[i])
                                                mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                                mu1[i]  <- mu[level[i]] + ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2)) +
                                                           beta_random[level[i]]*random[i]
                                                tau2[i] <- 1/(se[i]^2)}

                                              resid <- est-mu2

                                              sigma  ~ dunif(0, Ps)

                                              ##Stochastic behavior for the weights
                                              pw[1:npw] ~ ddirch(alpha_pw[1:npw])

                                              ##priors
                                              for(j in 1:L){
                                                beta_random[j] ~ dnorm(0, 1/0.5^2)
                                                beta_adjust[j]  ~ dnorm(0, 1/10^2)
                                                tau1[j]        <- 1/sigma^2

                                                for(k in 1:npw){
                                                  mu_M[j,k]    ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}
                                                d[j]         ~ dcat(pw[1:npw])
                                                mu[j]        <- inprod(mu_M[j, 1:npw], d[j] == 1:npw)}

                                              #I2 per level
                                              for(l in 1:L){
                                                I2[l] <- (1/tau1[level[l]])/(1/tau2[l]+1/tau1[level[l]])}}}
                                          ##1.3 more than 1 random
                                          else{
                                          meta_analysis <- function(){

                                              ##likelihood
                                              for (i in 1:N){
                                                est[i]  ~ dnorm(mu2[i], tau2[i])
                                                mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                                mu1[i]  <- mu[level[i]] +ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))+
                                                  inprod(beta_random[level[i], ], random[i, ])
                                                tau2[i] <- 1/(se[i]^2)}

                                              resid <- est-mu2

                                              sigma  ~ dunif(0, Ps)

                                              ##Stochastic behavior for the weights
                                              pw[1:npw] ~ ddirch(alpha_pw[1:npw])

                                              ##priors
                                              for(j in 1:L){
                                                for(r in 1:R){
                                                beta_random[j, r]  ~ dnorm(0, 1/0.5^2)}
                                                beta_adjust[j]      ~ dnorm(0, 1/10^2)
                                                tau1[j]            <- 1/sigma^2

                                                for(k in 1:npw){
                                                  mu_M[j,k]    ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}
                                                d[j]         ~ dcat(pw[1:npw])
                                                mu[j]        <- inprod(mu_M[j, 1:npw], d[j] == 1:npw)}

                                              #I2 per level
                                              for(l in 1:L){
                                                I2[l] <- (1/tau1[level[l]])/(1/tau2[l]+1/tau1[level[l]])}}}}
                                        else{
                                          ##2.##RE model with npw=0
                                          ##1. no random
                                          if(mod_data$R==0){
                                          meta_analysis <- function(){

                                            ##likelihood
                                            for (i in 1:N){
                                              est[i]  ~ dnorm(mu2[i], tau2[i])
                                              mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                              mu1[i]  <- mu[level[i]] + ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))
                                              tau2[i] <- 1/(se[i]^2)}

                                            resid <- est-mu2

                                            sigma  ~ dunif(0, Ps)

                                            ##priors
                                            for(j in 1:L){
                                              beta_adjust[j]   ~ dnorm(0, 1/10^2)
                                              mu[j]           ~ dnorm(Pm[j], 1/Pe[j]^2)
                                              tau1[j]         <- 1/sigma^2}

                                            #I2 per level
                                            for(l in 1:L){
                                              I2[l] <- (1/tau1[level[l]])/(1/tau2[l]+1/tau1[level[l]])}}}
                                          ## only 1 random
                                          else if(mod_data$R==1){
                                          meta_analysis <- function(){

                                              ##likelihood
                                              for (i in 1:N){
                                                est[i]  ~ dnorm(mu2[i], tau2[i])
                                                mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                                mu1[i]  <- mu[level[i]] + ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))+
                                                           beta_random[level[i]]*random[i]
                                                tau2[i] <- 1/(se[i]^2)}

                                              resid <- est-mu2

                                              sigma  ~ dunif(0, Ps)

                                              ##priors
                                              for(j in 1:L){
                                                beta_random[j]  ~ dnorm(0, 1/0.5^2)
                                                beta_adjust[j]   ~ dnorm(0, 1/10^2)
                                                mu[j]           ~ dnorm(Pm[j], 1/Pe[j]^2)
                                                tau1[j]         <- 1/sigma^2}

                                              #I2 per level
                                              for(l in 1:L){
                                                I2[l] <- (1/tau1[level[l]])/(1/tau2[l]+1/tau1[level[l]])}}}
                                          ## more than 1 random
                                          else{
                                          meta_analysis <- function(){

                                              ##likelihood
                                              for (i in 1:N){
                                                est[i]  ~ dnorm(mu2[i], tau2[i])
                                                mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                                mu1[i]  <- mu[level[i]] +ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))+
                                                  inprod(beta_random[level[i], ], random[i, ])
                                                tau2[i] <- 1/(se[i]^2)}

                                              resid <- est-mu2

                                              sigma  ~ dunif(0, Ps)

                                              ##priors
                                              for(j in 1:L){
                                                for(r in 1:R){
                                                beta_random[j, r]  ~ dnorm(0, 1/0.5^2)}
                                                beta_adjust[j]      ~ dnorm(0, 1/10^2)
                                                mu[j]              ~ dnorm(Pm[j], 1/Pe[j]^2)
                                                tau1[j]            <- 1/sigma^2}

                                              #I2 per level
                                              for(l in 1:L){
                                                I2[l] <- (1/tau1[level[l]])/(1/tau2[l]+1/tau1[level[l]])}}}}}
                                      else{
                                         if(mod_data$npw>1){
                                             ##3.##FE model with npw>1 and
                                             ##3.1 no random
                                             if(mod_data$R==0){
                                             meta_analysis <- function(){

                                             ##likelihood
                                             for (i in 1:N){
                                               est[i]  ~  dnorm(mu1[i], tau1[i])
                                               mu1[i]  <- mu[level[i]] + ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))
                                               tau1[i] <- 1/(se[i]^2)}

                                             resid <- est-mu1

                                             ##Stochastic behavior for the weights
                                             pw[1:npw] ~ ddirch(alpha_pw[1:npw])

                                             ##priors
                                             for(j in 1:L){
                                               beta_adjust[j]   ~ dnorm(0, 1/10^2)

                                               for(k in 1:npw){
                                                 mu_M[j,k]       ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}
                                               d[j]            ~ dcat(pw[1:npw])
                                               mu[j]           <- inprod(mu_M[j, 1:npw], d[j] == 1:npw)}}}
                                             ##3.2 only 1 random
                                             else if(mod_data$R==1){
                                             meta_analysis <- function(){

                                                 ##likelihood
                                                 for (i in 1:N){
                                                   est[i]  ~  dnorm(mu1[i], tau1[i])
                                                   mu1[i]  <- mu[level[i]]+ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))+
                                                              beta_random[level[i]]*random[i]
                                                   tau1[i] <- 1/(se[i]^2)}

                                                 resid <- est-mu1

                                                 ##Stochastic behavior for the weights
                                                 pw[1:npw] ~ ddirch(alpha_pw[1:npw])

                                                 ##priors
                                                 for(j in 1:L){
                                                   beta_random[j]  ~ dnorm(0, 1/0.5^2)
                                                   beta_adjust[j]   ~ dnorm(0, 1/10^2)

                                                   for(k in 1:npw){
                                                     mu_M[j,k]       ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}
                                                   d[j]            ~ dcat(pw[1:npw])
                                                   mu[j]           <- inprod(mu_M[j, 1:npw], d[j] == 1:npw)}}}
                                             ##3.3 more than 1 random
                                             else{
                                             meta_analysis <- function(){

                                                 ##likelihood
                                                 for (i in 1:N){
                                                   est[i]  ~  dnorm(mu1[i], tau1[i])
                                                   mu1[i]  <- mu[level[i]]+ ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))+
                                                     inprod(beta_random[level[i], ], random[i, ])
                                                   tau1[i] <- 1/(se[i]^2)}

                                                 resid <- est-mu1

                                                 ##Stochastic behavior for the weights
                                                 pw[1:npw] ~ ddirch(alpha_pw[1:npw])

                                                 ##priors
                                                 for(j in 1:L){
                                                   for(r in 1:R){
                                                     beta_random[j, r]  ~ dnorm(0, 1/0.5^2)}
                                                     beta_adjust[j]      ~ dnorm(0, 1/10^2)

                                                   for(k in 1:npw){
                                                     mu_M[j,k]       ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}
                                                   d[j]            ~ dcat(pw[1:npw])
                                                   mu[j]           <- inprod(mu_M[j, 1:npw], d[j] == 1:npw)}}}}
                                         else{
                                             ##4.##FE model with npw=0 and
                                             ##4.1 no random
                                             if(mod_data$R==0){
                                             meta_analysis <- function(){

                                                 ##likelihood
                                                 for (i in 1:N){
                                                   est[i]  ~  dnorm(mu1[i], tau1[i])
                                                   mu1[i]  <- mu[level[i]] + ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))
                                                   tau1[i] <- 1/(se[i]^2)}

                                                 resid <- est-mu1

                                                 ##priors
                                                 for(j in 1:L){
                                                   beta_adjust[j]   ~ dnorm(0, 1/10^2)
                                                   mu[j]           ~ dnorm(Pm[j], 1/Pe[j]^2)}}}
                                             ##4.2 only 1 random
                                             else if(mod_data$R==1){
                                             meta_analysis <- function(){

                                               ##likelihood
                                               for (i in 1:N){
                                                 est[i]  ~  dnorm(mu1[i], tau1[i])
                                                 mu1[i]  <- mu[level[i]]+ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2)) +
                                                            beta_random[level[i]]*random[i]
                                                 tau1[i] <- 1/(se[i]^2)}

                                               resid <- est-mu1

                                               ##priors
                                               for(j in 1:L){
                                                 beta_random[j]  ~ dnorm(0, 1/0.5^2)
                                                 beta_adjust[j]   ~ dnorm(0, 1/10^2)
                                                 mu[j]           ~ dnorm(Pm[j], 1/Pe[j]^2)}}}
                                             ##4.3 more than 1 random
                                             else{
                                             meta_analysis <- function(){

                                                 ##likelihood
                                                 for (i in 1:N){
                                                   est[i]  ~  dnorm(mu1[i], tau1[i])
                                                   mu1[i]  <- mu[level[i]]+ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2)) +
                                                     inprod(beta_random[level[i], ], random[i, ])
                                                   tau1[i] <- 1/(se[i]^2)}

                                                 resid <- est-mu1

                                                 ##priors
                                                 for(j in 1:L){
                                                   for(r in 1:R){
                                                   beta_random[j, r]  ~ dnorm(0, 1/0.5^2)}
                                                   beta_adjust[j]      ~ dnorm(0, 1/10^2)
                                                   mu[j]              ~ dnorm(Pm[j], 1/Pe[j]^2)}}}}}

                                         #Run the model
                                         model <- R2jags::jags.parallel(data = mod_data,
                                                                model.file = meta_analysis,
                                                                parameters.to.save =  c("mu", "d", "I2", "sigma", "beta_random", "beta_adjust", "resid"),
                                                                n.chains = n_chain,
                                                                n.thin = n_thin,
                                                                jags.seed = n_seed,
                                                                n.iter = n_iter,
                                                                n.burnin = n_burnin)

                                         #Warning messages for bad mixing chains
                                         if(any(c(model$BUGSoutput$summary[-1,8]>Rhat_warn, model$BUGSoutput$summary[-1,9]<Eff_warn))){warning("The Rhat or/and effective sample size for some >",Rhat_warn," or/and <",Eff_warn,". Chains might not be mixing.")}

                                         #Function to extract chains
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

                                         #Extract posterior residual means of the model
                                         mcmc_mu_resid     <- colMeans(model$BUGSoutput$sims.list$resid)

                                         #Extract chains for means
                                         mcmc_mu           <- extract_chain(model$BUGSoutput$sims.list$mu, mod_data)
                                         split_par         <- split(mcmc_mu, mcmc_mu$parameter)

                                         #summarize mu
                                         #ETI intervals
                                         #interval_level      <- 0.5+c(-1,1)*interval/2
                                         #HDI intervals
                                         hdi_fun              <- function(x, level){

                                           orddata   <- sort(x)
                                           nord      <- length(x)
                                           infomass  <- ceiling(level*nord)
                                           outmass   <- nord-infomass

                                           min_width <- Inf
                                           ll        <- NA
                                           ul        <- NA

                                           for(i in 1:(outmass+1)){
                                             int_width <- orddata[i+infomass-1]-orddata[i]

                                             if(int_width < min_width){
                                               min_width <- int_width
                                               ll <- orddata[i]
                                               ul <- orddata[i+infomass-1]}}

                                           c(ll, ul)}
                                         maxpost              <- function(x){d <- density(x); d$x[which.max(d$y)]}
                                         basic_summary        <- aggregate(data=split_par$b1, estimate~., function(x) c(maxpost(x), mean(x), sd(x),
                                                                                                                        hdi_fun(x, level = interval)))
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

                                         #Extract chains for peese or peters
                                         if(method!=0){
                                           mcmc_adjust        <- extract_chain(model$BUGSoutput$sims.list$beta_adjust, mod_data)}else{mcmc_adjust <- NULL}

                                         #Extract chains for posterior weights
                                         if(mod_data$npw>1){
                                           mcmc_podd       <- extract_chain(model$BUGSoutput$sims.list$d, mod_data)}else{mcmc_podd <- NULL}

                                         if(print_summary==T){
                                         print(basic_summary)}

                                         return(invisible(list(Summary=basic_summary,
                                                     Estimates=split(mcmc_mu, mcmc_mu$parameter),
                                                     Chains_mu=mcmc_mu,
                                                     Chains_I2=mcmc_I2,
                                                     Chains_sigma=mcmc_sigma,
                                                     Chains_adjust=mcmc_adjust,
                                                     Chains_podd=mcmc_podd,
                                                     Residuals=mcmc_mu_resid,
                                                     N_level=table(mod_data$level),
                                                     model=list(JAGS_model=model,
                                                                Data=mod_data, Priors=data.frame(Levels=levels(mod_data$level),
                                                                                                 Prior_mu=mod_data$Pm,
                                                                                                 Prior_se=mod_data$Pe)))))}}

