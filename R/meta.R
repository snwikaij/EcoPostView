#' Bayesian Meta-Analysis
#'
#' @param estimate Parameter b0 or b1 estimated from a LM or GLM
#' @param stderr Standard error belonging to b0 and b1
#' @param parameter A category (name) either b0 or b1
#' @param predictor  A predictor name (e.g., salinity) as multiple predictors can be handled
#' @param link_function The link function of the LM or GLM currently only 'identity', 'logit' and 'log' are supported and will only be used in the hop function
#' @param grouping A category (name) for the group as multiple groups can be handled
#' @param random An argument that needs to be a vector or  matrix of factors of the same number of rows as the estimate
#' @param method Indicates which adjustment performed 0 (='none'), 1 (='egger') or 2 (='peters')
#' @param Nsamp A vector with the number of samples for each estimate (only used when method = 2)
#' @param prior_mu Prior for the mean which can be vector (Bayesian meta-analysis) or matrix (Bayesian meta-analysis with model averaging)
#' @param prior_mu_se Prior for the se which can be vector or matrix
#' @param prior_study_var Prior for the between study variance if RE=T. If the argument prior_fam_var indicates the family is 'unif1' or 'unif2', then a uniform prior is used between 0 and the maximum value indicated by prior_study_var (default=abs(max(estimate-mean(estimate)))*2) meaning Uniform(0, prior_study_var).
#' If prior_fam_var is 'exp' the the rate indicated bye prior_study_var. This is by default 0.001 meaning Exponential(prior_study_var).
#' @param prior_fam_var The distribution family used for the prior of the between study variance. Uniform1 = 'unif1' this uses the uniform distribution as a 'scaling parameter over all studies'. This can be beneficial is one does not want to make
#' assumptions about the between study variance and a large number of levels is present with sometimes small sample sizes. The between study variance is then estimated as the overall between study variance, but not for each level individual.
#' Uniform2 = 'unif2' estimates the between study variance for each individual level. This is useful when one no assumptions wants to make over the between study variance and a large number of studies is present for each level.
#' The Exponential = 'exp' estimates the between study variance for each individual level based on the exponential distribution (default='unif1').This is beneficial if only a small number of levels and samples are present. With unif1 or unif2 the intervals can become
#' very broad with might be beneficial if one wants to be conservative. However, it can be a disadvantage because that the variability of the parameter gets overly broad (default = 'unif1')
#' @param fixed_prior by default the prior weights are treated as stochastically generated from a Dirichlet distribution (fixed_prior=FALSE). It is possible to treat them as fixed (fixed_prior=TRUE). Under fixed conditions all priors are appointed equal prior weights as 1/number of priors.
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
#' Full Bayesian Random-Effect meta-analytic method with the possiblity of using Bayesian Model Averaging. It additionally has the possibility - internally - to adjust
#' for publication bias with Egger adjustment using the standard error (*SE^2*; Moreno et al., 2009; Stanley  and Doucouliagos, 2013) or Peters adjustment (Moreno et al., 2009) using
#' the inverse of the sample size (N). The full use of the meta function is the combination with pdplot and hop function. Hence it  is used to analyse
#' the parameter estimates intercepts (b0) and regression coefficients (b1) on of LM or GLM models applied when the independent
#' variable is continues.
#'
#' This package was first develop to analyse the  parameters when the independent variable (predictor variable) was log (natural log)
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
#'Which prior variance 'family' needs to be selected on the between study variance is still unclear to me.
#'I am not sure if there is an approximate optimal answer but would like to have one. It is possible to get unreasonable answers
#'due to the extremely large or small variance among studies and variability in sample sizes when using 'unif2' or 'exp'.
#'However extremely wide intervals across all levels are also possible when using 'unif1'. This means then that smaller pattern will not
#'come out very clear.
#'
#'The 'unif1' could prevents over fitting to small-group variances for smaller sample sizes. Therefore it a conservative
#'but stable approximate answer seems preferable at first glance. Consequently repetition could be performed in an second study
#'with a less conservative hierarchical structure where 'unif2' and 'exp' could be utilized.
#'
#'If Nsamp (NA) this is replaced by 1 since 1/N where N=1 is 1=1/1.
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
                 moderator=NULL, method=0, Nsamp=NULL,
                 prior_mu=0,
                 prior_mu_se=10,
                 prior_study_var=NULL,
                 prior_fam_var="unif1",
                 fixed_prior=FALSE,
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

  argument.call <- match.call()

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
    npw <- 1}else if(length(prior_mu) == length(prior_mu_se) &&
                     is.null(dim(prior_mu)) && is.null(dim(prior_mu_se))){npw <- 1
    }else if(is.data.frame(prior_mu) && is.data.frame(prior_mu_se) &&
             ncol(prior_mu) == ncol(prior_mu_se)){npw <- ncol(prior_mu)
    }else{stop(paste0("The the priors are not nummeric or the priors for prior_mu and prior_mu_se are not of the same dimensions."))}

  #If prior weights are fixed and not stochastic npw > 1
  if(fixed_prior==T && npw == 1){stop("The number of priors provided needs to be > 1.")}
  if(fixed_prior==T){FP <- c(0,1)}else{FP <- c(1,0)}

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

  #A moderator data frame
  if(!is.null(moderator) && is.data.frame(moderator)){
    if(nrow(moderator)!=length(estimate)){stop("Moderator needs to be a numeric data frame or matrix with the same
                                           length as the estimates")}}

  if(is.null(moderator)){ML <- 1; MP <- 0; moderator <- rep(0, length(estimate))}else if(is.null(dim(moderator))){
    ML <- 1; MP <- 1; moderator <- as.numeric(moderator)}else{
      ML <- ncol(moderator); moderator<- apply(moderator, 2, function(x) as.numeric(x))}

  #A random effect data frame
  if(!is.null(random) && is.data.frame(random)){
    if(nrow(random)!=length(estimate)){stop("Random  effect  needs to be a data frame or matrix with
                                           factors of same length as the estimates")}}

  if(is.null(random)){RL <- 1; RP <- 0; random <- as.factor(rep(0, length(estimate)))}else if(is.null(dim(random))){
    RL <- 1; RP <- 1; random <- as.factor(random)}else{
      RL <- ncol(random); RP <- 1; random <- apply(random, 2, function(x) as.numeric(factor(x)))}

  #Set correction method for bias
  if(all(method != c(0, 1, 2))){
    stop("method needs to be either 0 (='none'), 1 (='egger'), 2 (='peters')")}
  if(method == 0){AP <- 0}else{AP <- 1}
  if(method == 2){
    if(is.null(Nsamp) | length(Nsamp) != length(estimate)){stop("If method is 'peters' then length of Nsamp needs to be of the same length as
                                             the estimates")}
    Nsamp[is.na(Nsamp)] <- 1}else{Nsamp <- rep(1, length(estimate))}

  #Set the type of prior for between study variance 1 (='unif1'), 2 (='unif2'), 2 (='exp')
    if(prior_fam_var=="unif1"){prior_fam_var <- c(1, 0, 0)
    if(is.null(prior_study_var)){prior_study_var <- max(abs(estimate-mean(estimate)))*2}
    }else if(prior_fam_var=="unif2"){prior_fam_var <- c(0, 1, 0)
    if(is.null(prior_study_var)){prior_study_var <- max(abs(estimate-mean(estimate)))*2}
    }else if(prior_fam_var=="exp"){prior_fam_var <- c(0, 0, 1)
    if(is.null(prior_study_var)){prior_study_var <- 0.001}
    }else{stop("Not a correct family choosen for between study variance: 'unif1', 'unif2' or 'exp'")}

  #Place all given data in a list
  mod_data <- list(est=as.numeric(estimate),
                   se=as.numeric(stderr),
                   N=length(estimate),
                   level=factor(level, levels = unique(level)[order(unique(level))]),
                   L=Ll,
                   Nsamp=Nsamp,
                   Pm=prior_mu,
                   Pe=prior_mu_se,
                   Ps=prior_study_var,
                   npw=npw,
                   alpha_pw=rep(1, npw),
                   fix_pw=rep(1/npw, npw),
                   FP=FP,
                   random=random,
                   R=RL,
                   RP=RP,
                   moderator=moderator,
                   M=ML,
                   MP=MP,
                   prV=prior_fam_var,
                   method=method,
                   AP=AP)

  #If only the prior is needed return
  if(get_prior_only){return(data.frame(Levels=levels(mod_data$level),
                                       Prior_mu=mod_data$Pm,
                                       Prior_se=mod_data$Pe))}else if(get_prior_only == F){

                                        if(npw>1){
                                         if(RL<=1 && ML<=1){
                                           meta_analysis <- function(){

                                           for (i in 1:N) {
                                             est[i]  ~ dnorm(mu2[i], tau2[i])
                                             mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                             mu1[i]  <- mu[level[i]] + beta_random[level[i]] * random[i] + beta_moderator[level[i]] * moderator[i] +
                                               ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))
                                             tau2[i] <- 1/(se[i]^2)}

                                           resid[1:N] <- est[1:N] - mu2[1:N]

                                           #Prior weights/odds
                                           pw_stochastic[1:npw] ~ ddirch(alpha_pw[1:npw])
                                           pw_fixed            <- fix_pw
                                           pw                  <- FP[1]*pw_stochastic+FP[2]*pw_fixed

                                           #All priors
                                           sigma2 ~ dunif(0, Ps)

                                           for (j in 1:L) {

                                             #random effects
                                               beta_R[j] ~ dnorm(0, 1/1000^2)
                                               beta_random[j] <- beta_R[j] * RP

                                             #Moderator
                                               beta_M[j] ~ dnorm(0, 1/1000^2)
                                               beta_moderator[j] <- beta_M[j] * MP

                                             #Bias adjustement
                                             beta_A[j]      ~ dnorm(0, 1/1000^2)
                                             beta_adjust[j] <- beta_A[j] * AP

                                             sigma1[j] ~ dunif(0, Ps)
                                             tau1a[j] <- 1 / sigma1[j]^2
                                             tau1b[j] <- 1 / sigma2^2
                                             tau1c[j] ~ dexp(Ps)
                                             tau1[j] <- tau1a[j]*prV[1] + tau1b[j]*prV[2] + tau1c[j]*prV[3]

                                             for(k in 1:npw){
                                               mu_M[j,k] ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}

                                             d[j]  ~ dcat(pw[1:npw])
                                             mu[j] <- inprod(mu_M[j, 1:npw], equals(d[j], 1:npw))
                                             }

                                           #I2
                                           for(l in 1:L) {
                                             I2[l] <- (1/tau1[level[l]]) / (1/tau2[l] + 1/tau1[level[l]])}
                                           }
                                         }else if(RL>1 && ML>1){
                                           meta_analysis <- function(){

                                               for (i in 1:N) {
                                                 est[i]  ~ dnorm(mu2[i], tau2[i])
                                                 mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                                 mu1[i]  <- mu[level[i]] + inprod(beta_random[level[i], 1:R], random[i, 1:R]) + inprod(beta_moderator[level[i], 1:M], moderator[i, 1:M]) +
                                                   ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))
                                                 tau2[i] <- 1/(se[i]^2)}

                                               resid[1:N] <- est[1:N] - mu2[1:N]

                                               #Prior weights/odds
                                               pw_stochastic[1:npw] ~ ddirch(alpha_pw[1:npw])
                                               pw_fixed            <- fix_pw
                                               pw                  <- FP[1]*pw_stochastic+FP[2]*pw_fixed

                                               #All priors
                                               sigma2 ~ dunif(0, Ps)

                                               for (j in 1:L) {

                                                 #random effects
                                                 for(r in 1:R){
                                                   beta_R[j,r]  ~ dnorm(0, 1/1000^2)
                                                   beta_random[j,r] <- beta_R[j,r] * RP}

                                                 #Moderator
                                                 for(m in 1:M){
                                                   beta_M[j,m] ~ dnorm(0, 1/1000^2)
                                                   beta_moderator[j,m] <- beta_M[j,m] * MP}

                                                 #Bias adjustement
                                                 beta_A[j]      ~ dnorm(0, 1/1000^2)
                                                 beta_adjust[j] <- beta_A[j] * AP

                                                 sigma1[j] ~ dunif(0, Ps)
                                                 tau1a[j] <- 1 / sigma1[j]^2
                                                 tau1b[j] <- 1 / sigma2^2
                                                 tau1c[j] ~ dexp(Ps)
                                                 tau1[j] <- tau1a[j]*prV[1] + tau1b[j]*prV[2] + tau1c[j]*prV[3]

                                                 for(k in 1:npw){
                                                   mu_M[j,k] ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}

                                                 d[j]  ~ dcat(pw[1:npw])
                                                 mu[j] <- inprod(mu_M[j, 1:npw], equals(d[j], 1:npw))
                                               }

                                               #I2
                                               for(l in 1:L) {
                                                 I2[l] <- (1/tau1[level[l]]) / (1/tau2[l] + 1/tau1[level[l]])}
                                             }
                                         }else if(RL>1 && ML<=1){
                                           meta_analysis <- function(){

                                               for (i in 1:N) {
                                                 est[i]  ~ dnorm(mu2[i], tau2[i])
                                                 mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                                 mu1[i]  <- mu[level[i]] + inprod(beta_random[level[i], 1:R], random[i, 1:R]) + beta_moderator[level[i]] * moderator[i] +
                                                   ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))
                                                 tau2[i] <- 1/(se[i]^2)}

                                               resid[1:N] <- est[1:N] - mu2[1:N]

                                               #Prior weights/odds
                                               pw_stochastic[1:npw] ~ ddirch(alpha_pw[1:npw])
                                               pw_fixed            <- fix_pw
                                               pw                  <- FP[1]*pw_stochastic+FP[2]*pw_fixed

                                               #All priors
                                               sigma2 ~ dunif(0, Ps)

                                               for (j in 1:L) {

                                                 #random effects
                                                 for(r in 1:R){
                                                   beta_R[j,r]  ~ dnorm(0, 1/1000^2)
                                                   beta_random[j,r] <- beta_R[j,r] * RP}

                                                 #Moderator
                                                 beta_M[j] ~ dnorm(0, 1/1000^2)
                                                 beta_moderator[j] <- beta_M[j] * MP

                                                 #Bias adjustement
                                                 beta_A[j]      ~ dnorm(0, 1/1000^2)
                                                 beta_adjust[j] <- beta_A[j] * AP

                                                 sigma1[j] ~ dunif(0, Ps)
                                                 tau1a[j] <- 1 / sigma1[j]^2
                                                 tau1b[j] <- 1 / sigma2^2
                                                 tau1c[j] ~ dexp(Ps)
                                                 tau1[j] <- tau1a[j]*prV[1] + tau1b[j]*prV[2] + tau1c[j]*prV[3]

                                                 for(k in 1:npw){
                                                   mu_M[j,k] ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}

                                                 d[j]  ~ dcat(pw[1:npw])
                                                 mu[j] <- inprod(mu_M[j, 1:npw], equals(d[j], 1:npw))
                                               }

                                               #I2
                                               for(l in 1:L) {
                                                 I2[l] <- (1/tau1[level[l]]) / (1/tau2[l] + 1/tau1[level[l]])}
                                             }
                                         }else if(RL<=1 && ML>1){
                                           meta_analysis <- function(){

                                             for (i in 1:N) {
                                               est[i]  ~ dnorm(mu2[i], tau2[i])
                                               mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                               mu1[i]  <- mu[level[i]] + beta_random[level[i]] * random[i] +  inprod(beta_moderator[level[i], 1:M], moderator[i, 1:M]) +
                                                 ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))
                                               tau2[i] <- 1/(se[i]^2)}

                                             resid[1:N] <- est[1:N] - mu2[1:N]

                                             #Prior weights/odds
                                             pw_stochastic[1:npw] ~ ddirch(alpha_pw[1:npw])
                                             pw_fixed            <- fix_pw
                                             pw                  <- FP[1]*pw_stochastic+FP[2]*pw_fixed

                                             #All priors
                                             sigma2 ~ dunif(0, Ps)

                                             for (j in 1:L) {

                                               #random effects
                                               beta_R[j] ~ dnorm(0, 1/1000^2)
                                               beta_random[j] <- beta_R[j] * RP

                                               #Moderator
                                               for(m in 1:M){
                                                 beta_M[j,m] ~ dnorm(0, 1/1000^2)
                                                 beta_moderator[j,m] <- beta_M[j,m] * MP}

                                               #Bias adjustement
                                               beta_A[j]      ~ dnorm(0, 1/1000^2)
                                               beta_adjust[j] <- beta_A[j] * AP

                                               sigma1[j] ~ dunif(0, Ps)
                                               tau1a[j] <- 1 / sigma1[j]^2
                                               tau1b[j] <- 1 / sigma2^2
                                               tau1c[j] ~ dexp(Ps)
                                               tau1[j] <- tau1a[j]*prV[1] + tau1b[j]*prV[2] + tau1c[j]*prV[3]

                                               for(k in 1:npw){
                                                 mu_M[j,k] ~ dnorm(Pm[j,k], 1/Pe[j,k]^2)}

                                               d[j]  ~ dcat(pw[1:npw])
                                               mu[j] <- inprod(mu_M[j, 1:npw], equals(d[j], 1:npw))
                                             }

                                             #I2
                                             for(l in 1:L) {
                                               I2[l] <- (1/tau1[level[l]]) / (1/tau2[l] + 1/tau1[level[l]])}
                                           }
                                         }
                                        }else{#If npw=1
                                          if(RL<=1 && ML<=1){
                                           meta_analysis <- function(){

                                             #Likelihood
                                             for (i in 1:N) {
                                               est[i]  ~ dnorm(mu2[i], tau2[i])
                                               mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                               mu1[i]  <- mu[level[i]] + beta_random[level[i]] * random[i] + beta_moderator[level[i]] * moderator[i] +
                                                          ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))
                                               tau2[i] <- 1/(se[i]^2)}

                                             resid[1:N] <- est[1:N] - mu2[1:N]

                                             #All priors
                                             sigma2 ~ dunif(0, Ps)

                                             for (j in 1:L) {

                                               #random effects
                                                 beta_R[j] ~ dnorm(0, 1/1000^2)
                                                 beta_random[j] <- beta_R[j] * RP

                                               #Moderator
                                                 beta_M[j] ~ dnorm(0, 1/1000^2)
                                                 beta_moderator[j] <- beta_M[j] * MP

                                               #Bias adjustement
                                               beta_A[j]      ~ dnorm(0, 1/1000^2)
                                               beta_adjust[j] <- beta_A[j] * AP

                                               mu[j]           ~ dnorm(Pm[j], 1/Pe[j]^2)

                                               sigma1[j] ~ dunif(0, Ps)
                                               tau1a[j] <- 1 / sigma1[j]^2
                                               tau1b[j] <- 1 / sigma2^2
                                               tau1c[j] ~ dexp(Ps)
                                               tau1[j] <- tau1a[j]*prV[1] + tau1b[j]*prV[2] + tau1c[j]*prV[3]}

                                             #I2 heterogenity
                                             for(l in 1:L) {
                                               I2[l] <- (1/tau1[level[l]]) / (1/tau2[l] + 1/tau1[level[l]])
                                             }
                                           }
                                          }else if(RL>1 && ML>1){
                                           meta_analysis <- function(){

                                             #Likelihood
                                             for (i in 1:N) {
                                               est[i]  ~ dnorm(mu2[i], tau2[i])
                                               mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                               mu1[i]  <- mu[level[i]] + inprod(beta_random[level[i], 1:R], random[i, 1:R]) + inprod(beta_moderate[level[i], 1:M], moderate[i, 1:M]) +
                                                 ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))
                                               tau2[i] <- 1/(se[i]^2)}

                                             resid[1:N] <- est[1:N] - mu2[1:N]

                                             #All priors
                                             sigma2 ~ dunif(0, Ps)

                                             for (j in 1:L) {

                                               #random effects
                                               for(r in 1:R){
                                                 beta_R[j,r]  ~ dnorm(0, 1/1000^2)
                                                 beta_random[j,r] <- beta_R[j,r] * RP}

                                               #Moderator
                                               for(m in 1:M){
                                                 beta_M[j,m] ~ dnorm(0, 1/1000^2)
                                                 beta_moderator[j,m] <- beta_M[j,m] * MP}

                                               #Bias adjustement
                                               beta_A[j]      ~ dnorm(0, 1/1000^2)
                                               beta_adjust[j] <- beta_A[j] * AP

                                               mu[j]           ~ dnorm(Pm[j], 1/Pe[j]^2)

                                               sigma1[j] ~ dunif(0, Ps)
                                               tau1a[j] <- 1 / sigma1[j]^2
                                               tau1b[j] <- 1 / sigma2^2
                                               tau1c[j] ~ dexp(Ps)
                                               tau1[j] <- tau1a[j]*prV[1] + tau1b[j]*prV[2] + tau1c[j]*prV[3]}

                                             #I2 heterogenity
                                             for(l in 1:L) {
                                               I2[l] <- (1/tau1[level[l]]) / (1/tau2[l] + 1/tau1[level[l]])
                                             }
                                           }
                                          }else if(RL>1 && ML<=1){
                                           meta_analysis <- function(){

                                               #Likelihood
                                               for (i in 1:N) {
                                                 est[i]  ~ dnorm(mu2[i], tau2[i])
                                                 mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                                 mu1[i]  <- mu[level[i]] + inprod(beta_random[level[i], 1:R], random[i, 1:R]) +  beta_moderator[level[i]] * moderator[i] +
                                                   ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))
                                                 tau2[i] <- 1/(se[i]^2)}

                                               resid[1:N] <- est[1:N] - mu2[1:N]

                                               #All priors
                                               sigma2 ~ dunif(0, Ps)

                                               for (j in 1:L) {

                                                 #random effects
                                                 for(r in 1:R){
                                                   beta_R[j,r]  ~ dnorm(0, 1/1000^2)
                                                   beta_random[j,r] <- beta_R[j,r] * RP}

                                                 #Moderator
                                                 beta_M[j] ~ dnorm(0, 1/1000^2)
                                                 beta_moderator[j] <- beta_M[j] * MP

                                                 #Bias adjustement
                                                 beta_A[j]      ~ dnorm(0, 1/1000^2)
                                                 beta_adjust[j] <- beta_A[j] * AP

                                                 mu[j]           ~ dnorm(Pm[j], 1/Pe[j]^2)

                                                 sigma1[j] ~ dunif(0, Ps)
                                                 tau1a[j] <- 1 / sigma1[j]^2
                                                 tau1b[j] <- 1 / sigma2^2
                                                 tau1c[j] ~ dexp(Ps)
                                                 tau1[j] <- tau1a[j]*prV[1] + tau1b[j]*prV[2] + tau1c[j]*prV[3]}

                                               #I2 heterogenity
                                               for(l in 1:L) {
                                                 I2[l] <- (1/tau1[level[l]]) / (1/tau2[l] + 1/tau1[level[l]])
                                               }
                                             }
                                          }else if(RL<=1 && ML>1){
                                           meta_analysis <- function(){

                                              #Likelihood
                                              for (i in 1:N) {
                                                est[i]  ~ dnorm(mu2[i], tau2[i])
                                                mu2[i]  ~ dnorm(mu1[i], tau1[level[i]])
                                                mu1[i]  <- mu[level[i]] + beta_random[level[i]] * random[i] + inprod(beta_moderate[level[i], 1:M], moderate[i, 1:M]) +
                                                  ifelse(method==0, 0, ifelse(method==1, beta_adjust[level[i]]*se[i]^2, beta_adjust[level[i]]*Nsamp[i]^2))
                                                tau2[i] <- 1/(se[i]^2)}

                                              resid[1:N] <- est[1:N] - mu2[1:N]

                                              #All priors
                                              sigma2 ~ dunif(0, Ps)

                                              for (j in 1:L) {

                                                #random effects
                                                beta_R[j] ~ dnorm(0, 1/1000^2)
                                                beta_random[j] <- beta_R[j] * RP

                                                #Moderator
                                                for(m in 1:M){
                                                  beta_M[j,m] ~ dnorm(0, 1/1000^2)
                                                  beta_moderator[j,m] <- beta_M[j,m] * MP}

                                                #Bias adjustement
                                                beta_A[j]      ~ dnorm(0, 1/1000^2)
                                                beta_adjust[j] <- beta_A[j] * AP

                                                mu[j]           ~ dnorm(Pm[j], 1/Pe[j]^2)

                                                sigma1[j] ~ dunif(0, Ps)
                                                tau1a[j] <- 1 / sigma1[j]^2
                                                tau1b[j] <- 1 / sigma2^2
                                                tau1c[j] ~ dexp(Ps)
                                                tau1[j] <- tau1a[j]*prV[1] + tau1b[j]*prV[2] + tau1c[j]*prV[3]}

                                              #I2 heterogenity
                                              for(l in 1:L) {
                                                I2[l] <- (1/tau1[level[l]]) / (1/tau2[l] + 1/tau1[level[l]])
                                              }
                                            }
                                          }
                                         }

                                         #Run the model
                                         model <- R2jags::jags.parallel(data = mod_data,
                                                                        model.file = meta_analysis,
                                                                        parameters.to.save =  c("mu", "d", "I2", "sigma1", "sigma2", "beta_random", "beta_moderator", "beta_adjust", "resid"),
                                                                        n.chains = n_chain,
                                                                        n.thin = n_thin,
                                                                        jags.seed = n_seed,
                                                                        n.iter = n_iter,
                                                                        n.burnin = n_burnin)

                                         #Warning messages for bad mixing chains
                                         amr       <- c("beta_adjust", "beta_moderator", "beta_random")[c(1*mod_data$AP, 2*mod_data$MP, 3*mod_data$RP)]
                                         which_par <- c("mu", "sigma1", "sigma2", "I2", amr)

                                         if(any(c(model$BUGSoutput$summary[rownames(model$BUGSoutput$summary) %in% which_par,8]>Rhat_warn, model$BUGSoutput$summary[rownames(model$BUGSoutput$summary) %in% which_par,9]<Eff_warn))){warning("The Rhat or/and effective sample size for some >",Rhat_warn," or/and <",Eff_warn,". Chains might not be mixing.")}

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

                                         #Mode or MAP
                                         maxpost              <- function(x){d <- density(x); d$x[which.max(d$y)]}

                                         #Generate summary table
                                         basic_summary        <- aggregate(data=split_par$b1, estimate~., function(x) c(maxpost(x), mean(x), sd(x), hdi_fun(x, level = interval)))

                                         #Derive number of studies via factors
                                         b1_n          <- as.data.frame(table(mod_data$level))
                                         b1_n          <- setNames(cbind(do.call(rbind.data.frame, strsplit(as.character(b1_n[,1]), "_")), b1_n[,2]), c("parameter", "predictor", "link", "group"))

                                         mcmc_I2       <- extract_chain(model$BUGSoutput$sims.list$I2, mod_data)
                                         mcmc_I2_b1    <- mcmc_I2[mcmc_I2$parameter!="b0",]
                                         mcmc_sigma    <- model$BUGSoutput$sims.list$sigma

                                         i2_summary                    <- aggregate(data=mcmc_I2_b1, estimate~., mean)$estimate
                                         basic_summary                 <- cbind(basic_summary[c(1:4)], round(unlist(basic_summary$estimate),4), round(i2_summary, 4))
                                         basic_summary                 <- merge(basic_summary, b1_n, by=c("parameter", "predictor", "link", "group"))
                                         colnames(basic_summary)[5:11] <- c("map", "mu", "se", "ll", "ul", "I2", "n")

                                         #Extract chains for peese or peters
                                         if(AP==1){
                                           mcmc_adjust     <- extract_chain(model$BUGSoutput$sims.list$beta_adjust, mod_data)
                                           summary_adjust  <- aggregate(data=mcmc_adjust, estimate~., function(x) c(maxpost(x), mean(x), sd(x), hdi_fun(x, level = interval)))
                                           summary_adjust  <- cbind(summary_adjust[-5], summary_adjust$estimate)
                                           colnames(summary_adjust)[5:9] <- c("map", "mu", "se", "ll", "ul")
                                         }else{mcmc_adjust <- NULL; summary_adjust  <- NULL}

                                         #Extract chains for moderator
                                         if(MP==1){
                                           mcmc_moderator     <- extract_chain(model$BUGSoutput$sims.list$beta_moderator, mod_data)
                                           summary_moderator  <- aggregate(data=mcmc_moderator, estimate~., function(x) c(maxpost(x), mean(x), sd(x), hdi_fun(x, level = interval)))
                                           summary_moderator  <- cbind(summary_moderator[-5], summary_moderator$estimate)
                                           colnames(summary_moderator)[5:9] <- c("map", "mu", "se", "ll", "ul")
                                         }else{mcmc_moderator <- NULL; summary_moderator  <- NULL}

                                         #Extract chains for random effect
                                         if(RP==1){
                                           mcmc_random     <- extract_chain(model$BUGSoutput$sims.list$beta_random, mod_data)
                                           summary_random  <- aggregate(data=mcmc_random, estimate~., function(x) c(maxpost(x), mean(x), sd(x), hdi_fun(x, level = interval)))
                                           summary_random  <- cbind(summary_random[-5], summary_random$estimate)
                                           colnames(summary_random)[5:9] <- c("map", "mu", "se", "ll", "ul")
                                        }else{mcmc_random <- NULL; summary_random  <- NULL}

                                         #Extract chains for posterior weights
                                         if(mod_data$npw>1){
                                           mcmc_podd       <- extract_chain(model$BUGSoutput$sims.list$d, mod_data)}else{mcmc_podd <- NULL}

                                         #Print summary
                                         if(print_summary==T){
                                           cat("Call:\n")
                                           print(argument.call)
                                           cat("\n")
                                           cat("Summary:\n")
                                           print(basic_summary)}

                                         total_summary <- list(Main=basic_summary,
                                                               Adjust=summary_adjust,
                                                               Moderator=summary_moderator,
                                                               Random=summary_random)

                                         total_summary <- total_summary[!sapply(total_summary, is.null)]

                                         tbls <- lapply(names(total_summary), function(name) {
                                           x <- total_summary[[name]]
                                           if (!is.null(x)) {
                                             x$Name <- name
                                             return(x)
                                           } else {
                                             return(NULL)
                                           } })

                                         tbls <- Filter(Negate(is.null), tbls)

                                         all_cols <- unique(unlist(lapply(tbls, names)))

                                           tbls_aligned <- lapply(tbls, function(df) {
                                           missing <- setdiff(all_cols, names(df))
                                           if (length(missing) > 0) {
                                             for (m in missing) df[[m]] <- NA
                                           }

                                           df <- df[all_cols]
                                           df})

                                         total_summary <- do.call(rbind, tbls_aligned)

                                         total_summary <- total_summary[c("Name", setdiff(names(total_summary), "Name"))]

                                         return(invisible(list(Summary=total_summary,
                                                               Estimates=split(mcmc_mu, mcmc_mu$parameter),
                                                               Chains_mu=mcmc_mu,
                                                               Chains_I2=mcmc_I2,
                                                               Chains_sigma=mcmc_sigma,
                                                               Chains_adjust=mcmc_adjust,
                                                               Chains_moderator=mcmc_moderator,
                                                               Chains_random=mcmc_random,
                                                               Chains_podd=mcmc_podd,
                                                               Residuals=mcmc_mu_resid,
                                                               N_level=table(mod_data$level),
                                                               model=list(JAGS_model=model,
                                                                          Data=mod_data, Priors=data.frame(Levels=levels(mod_data$level),
                                                                                                           Prior_mu=mod_data$Pm,
                                                                                                           Prior_se=mod_data$Pe)))))}}
