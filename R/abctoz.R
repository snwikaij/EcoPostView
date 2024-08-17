#' Aproximated Bayesian Computation to derive Z (ABC-to-z)
#'
#' @param z A numeric vector of z values
#' @param interval Credibility level at 0.9
#' @param nsim Number of simulations
#' @param accept_max_dist Threshold for minimal acceptable distance
#' @param prior_mu Prior values for mu(Z) following a truncated (0) normal distribution, which needs to notated c(mean, sd)
#' @param prior_sd Prior values for sd(Z) following a gamma distribution, which needs to be notated c(shape, rate)
#' @param prior_cens Prior values for the fractions of censored samples following a beta distribution, which needs to be notated as c(alpha, beta)
#' @param n_dens Number of density lines simulated from the posterior values mu(Z), sd(Z) and censored fraction (f)
#' @param print_progress Print the progress of simulations
#' @param seed Set the seed
#'
abctoz <- function(z, interval=0.9,
                   nsim=50000,
                   accept_max_dist=0.4,
                   prior_mu=c(1.5, 0.35),
                   prior_sd=c(0.4, 0.4),
                   prior_cens=c(4, 2),
                   n_dens=100,
                   print_progress=T,
                   seed=123){

#number of samples and prior array
z                <- abs(z)
nz               <- length(z)
priors           <- as.data.frame(array(NA, dim=c(nsim, 6)))
colnames(priors) <- c("mu",  "sd", "cens", "n1", "n2", "dist")

#all priors
set.seed(seed)
priors[,1] <- truncnorm::rtruncnorm(nsim, a = 0, b=Inf, mean=prior_mu[1], sd=prior_mu[2])
priors[,2] <- rgamma(nsim, prior_sd[1], prior_sd[2])
priors[,3] <- rbeta(nsim, prior_cens[1], prior_cens[2])
priors[,4] <- round(nz*priors[,3])
priors[,5] <- nz-priors[,4]

mod_dens <- function(x){d <- density(x); d$x[which.max(d$y)]}

#parameters of the data
data_mu  <- mean(z)
data_mod <- mod_dens(z)
data_sd  <- sd(z)
data_q   <- quantile(z, c(.25, .75))

#create some free space (is it at all needed?)
gc()

#set seed
set.seed(seed)

#ABC-rejection algorithm
for(i in 1:nrow(priors)){

if(print_progress == T) print(i)

sim       <- c(truncnorm::rtruncnorm(priors$n1[i], a=1.96, b=Inf, mean=priors$mu[i], sd=priors$sd[i]),
               truncnorm::rtruncnorm(priors$n2[i], a=0, b=Inf, mean=priors$mu[i], sd=priors$sd[i]))

priors$dist[i]    <- sum(abs(c(mean(sim), mod_dens(sim), sd(sim), quantile(sim, c(.25, .75)))-c(data_mu, data_mod, data_sd, data_q)))/5}

#Extract the posterior
priors           <- priors[order(priors$dist),]
posterior        <- priors[priors$dist<accept_max_dist,]

#local (generalized) linear correction function
loc_glm_cor <- function(par, dist, fam="Gamma"){

dfloc <- data.frame(par=par, dist=dist)

if(fam=="Gamma"){
wt  <- 1/abs(par-mean(par))
raw <- predict(glmmTMB(par~dist, family=Gamma(link="identity"), data=posterior), type = "response")
wtp <- predict(glmmTMB(par~dist, family=Gamma(link="identity"), data=posterior, weights = wt), type = "response")

}else if(fam=="Beta"){
wt  <- 1/abs(par-mean(par))
raw <- predict(glmmTMB(par~dist, family=beta_family(link="logit"), data=posterior), type = "response")
wtp <- predict(glmmTMB(par~dist,  family=beta_family(link="logit"), data=posterior, weights = wt), type = "response")}

par+(wtp-raw)}

#apply loc glm correction to mu, sd and censor
posterior$mu_adj   <- loc_glm_cor(posterior$mu, posterior$dist)
posterior$sd_adj   <- loc_glm_cor(posterior$sd, posterior$dist)
posterior$cens_adj <- loc_glm_cor(posterior$cens, posterior$dist, fam = "Beta")

#Plot the curvature of the number of iterations and log(distance)
distcurve <- ggplot(data.frame(x=log(priors$dist), y=1:nrow(priors)), aes(x, y))+
  geom_line()+xlab("log(average absolute distance)")+ylab("iteration")+
  geom_vline(xintercept = log(accept_max_dist), col="tomato3", lwd=0.6, lty=2)+
  theme_bw()

#Create a list for the density curves
sim_dens <- list()

#Simulate density curves
for(i in sample(1:nrow(posterior), size=100, replace = T)){

sim_dens[[i]]    <- c(truncnorm::rtruncnorm(posterior$n1[i], a=1.96, b=Inf, mean=posterior$mu_adj[i], sd=posterior$sd_adj[i]),
               truncnorm::rtruncnorm(posterior$n2[i], a=0, b=Inf, mean=posterior$mu_adj[i], sd=posterior$sd_adj[i]))}

#Unlist the simulated values
sim_z <- unlist(sim_dens)

#create n_dens density lines
densline <- list()
for(i in 1:100){
  x <- sample(sim_z, nz)
  sim_dens <- density(x, bw=0.1)
  dens_df  <- data.frame(i=i, x=sim_dens$x, y=sim_dens$y)
  densline[[i]] <- dens_df}

#Create a long format
densline <- do.call(rbind, densline)
densline <- densline[densline$x>0,]

#Create a density curve for the data
zdens <- data.frame(x=density(z, bw=0.2)$x, y=density(z, bw=0.2)$y)
zdens <- zdens[zdens$x>0,]

#Plot the density lines
pldens <- ggplot(densline, aes(x, y, group=as.factor(i)))+
  geom_line(alpha=0.5, col="grey80")+xlab("z-value")+
  ylab("Density")+xlim(0, 10)+
  geom_line(data=zdens, aes(x, y), inherit.aes = F)+
  #geom_vline(xintercept = qnorm(1 - sig_level/2), lty=2, col="tomato3")+
  theme_classic()+scale_alpha(guide = 'none')

#Use selected interval level
interval_level              <- 0.5+c(-1,1)*interval/2

#Generate a simple summary
output <- data.frame(
  Statistic = c("c", "mu(z)", "sd(z)"),
  Mean = round(c(mean(posterior$cens_adj), mean(posterior$mu_adj), mean(posterior$sd_adj)), 4),
  SE = round(c(sd(posterior$cens_adj), sd(posterior$mu_adj), sd(posterior$sd_adj)), 4),
  ll = round(c(quantile(posterior$cens_adj, interval_level[1]),
               quantile(posterior$mu_adj, interval_level[1]),
               quantile(posterior$sd_adj, interval_level[1])), 4),
  ul = round(c(quantile(posterior$cens_adj, interval_level[2]),
               quantile(posterior$mu_adj, interval_level[2]),
               quantile(posterior$sd_adj, interval_level[2])), 4))

#Select the relevant values
posterior_subset <- posterior[c("mu_adj", "sd_adj", "cens_adj")]
posterior_long   <- tidyr::gather(posterior_subset)

#Create a df for the point and intervals
table_hist <- output[c(1,2,4,5)]
colnames(table_hist)[1] <- "key"
table_hist$key <- c('cens_adj', 'mu_adj', 'sd_adj')

#Create text for labels in facet wraps
create_labs <- function(mu, ll, ul, text){
paste0(text,mu," (ll=",ll,", ","ul=",ul,")")}

#Change the names of the facet wraps
facet_labels <- c('cens_adj' = create_labs(round(output$Mean[1],2), round(output$ll[1],2), round(output$ul[1],2), "Censored fraction="),
                  'mu_adj' = create_labs(round(output$Mean[2],2), round(output$ll[2],2), round(output$ul[2],2), "Mean(Z)"),
                  'sd_adj' = create_labs(round(output$Mean[3],2), round(output$ll[3],2), round(output$ul[3],2), "SD(Z)"))

#Plot a histogram of the posterior iterations
#Histogram seems easier than a density plot
#due to the sometimes sparse results during
#optimization of fitting different priors and
#assessing the fit
plhist <- ggplot(posterior_long, aes(value))+
  geom_histogram(fill="grey80", alpha=0.2, col="black",
                 position="identity",
                 boundary = 0, closed = "left",
                 bins=20)+
  geom_point(data=table_hist, aes(x=Mean, y=0), size=4, inherit.aes = F)+
  geom_errorbarh(data=table_hist, aes(xmin=ll, xmax=ul, y=0),
                 height=0, lwd=1.2, inherit.aes = F)+
  facet_wrap(key~., scales = "free", labeller = labeller(key = facet_labels))+
  theme_classic()+ylab("Count")+
  theme(axis.title.x = element_blank())

return(list(summary=output,
            iterationcurve=distcurve,
            densitycurve=pldens,
            posterior=plhist,
            posterior_and_priors=list(priors=priors,
            posterior=posterior)))}
