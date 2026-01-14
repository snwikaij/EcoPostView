pred_hyp     <- function(object, variable=NULL, input_values,
                         family="normal", interval_type="CI",
                         prior_mu1=NULL, prior_sd1=NULL,
                         prior_mu0=NULL, prior_sd0=NULL,
                         level=0.9, xlim=NULL){

  if(!interval_type %in% c("PI", "CI")){stop("interval_type needs to be `PI` or `CI`")}

  #Check if all priors are used
  if(any(sapply(list(prior_mu1, prior_mu0, prior_sd1, prior_sd0), is.null))){stop("mu and sd for both priors need to be appointed.")}

  #Check if all roots are named
  if(!all(names(input_values) %in% object$Roots)){stop("Not the correct names as input values.")}

  pred <- predict_ppmn(object, new_data=input_values, interval_type = interval_type)

  #Check if there is only but only variable.
  if(is.null(variable)){stop("Select one variable.")}else if(
    length(variable)>1){stop("Only one variable allowed.")}

  #Check if the variable is in there
  pred_var <- names(pred$Variance)
  if(any(pred_var == variable)==F){stop("Variable name not found.")}

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

  posterior <- as.numeric(unlist(pred$Variance[[variable]]))
  posterior <- na.omit(posterior)
  post_pred <- c(mean(posterior, na.rm = T), sd(posterior, na.rm = T))

  priors_mu <- c(prior_mu1, prior_mu0)
  priors_sd <- c(prior_sd1, prior_sd0)

  if(family == "normal"){

    grid_range <- lapply(1:3, function(i){
      mu <- c(post_pred[1], priors_mu)[i]
      sd <- c(post_pred[2], priors_sd)[i]

      qnorm(c(0.001,0.999), mu, sd)})

    seq_min    <-  min(unlist(grid_range), na.rm = T)
    seq_max    <-  max(unlist(grid_range), na.rm = T)

    grid       <- seq(seq_min, seq_max, length.out=300000)
    dens       <- density(posterior)

    likelihood <- approx(dens$x, dens$y, rule = 2, xout = grid)$y
    prior1     <- dnorm(grid, priors_mu[1], priors_sd[1])
    prior2     <- dnorm(grid, priors_mu[2], priors_sd[2])

    m1 <- sum(likelihood * prior1)*(grid[2]-grid[1])
    m2 <- sum(likelihood * prior2)*(grid[2]-grid[1])

    results <-  c(m1/m2, "map"=maxpost(posterior), hdi_fun(posterior, level=level))
    xlims   <-  c(seq_min, seq_max)

  }else if(family=="gamma"){

    grid_range <- lapply(1:3, function(i){
      mu <- c(post_pred[1], priors_mu)[i]
      sd <- c(post_pred[2], priors_sd)[i]

      qgamma(c(0.001,0.999), mu^2/sd^2, mu/sd^2)})

    seq_min    <-  min(unlist(grid_range), na.rm = T)
    seq_max    <-  max(unlist(grid_range), na.rm = T)

    grid       <- seq(seq_min, seq_max, length.out=300000)
    dens       <- density(posterior)

    likelihood <- approx(dens$x, dens$y, rule = 2, xout = grid)$y
    prior1     <- dgamma(grid, priors_mu[1]^2/priors_sd[1]^2, priors_mu[1]/priors_sd[1]^2)
    prior2     <- dgamma(grid, priors_mu[2]^2/priors_sd[2]^2, priors_mu[2]/priors_sd[2]^2)

    m1 <- sum(likelihood * prior1)*(grid[2]-grid[1])
    m2 <- sum(likelihood * prior2)*(grid[2]-grid[1])

    results <-  c(m1/m2, "map"=maxpost(posterior), hdi_fun(posterior, level=level))
    xlims   <-  c(seq_min, seq_max)

  }else if(family=="beta"){

    if(any(priors_mu[1] <= 0 | priors_mu[1] >= 1)){stop("When using the familiy beta the priors for mu1 should be >0 or <1.")}
    if(any(priors_mu[2] <= 0 | priors_mu[2] >= 1)){stop("When using the familiy beta the priors for mu0 should be >0 or <1.")}

    if(any(priors_sd[1] <= 0 | priors_sd[1] >= 1)){stop("When using the familiy beta the priors for sd1 should be >0 or <1.")}
    if(any(priors_sd[2] <= 0 | priors_sd[2] >= 1)){stop("When using the familiy beta the priors for sd0 should be >0 or <1.")}

    if (priors_sd[1]^2 > priors_mu[1] * (1 - priors_mu[1])){stop("the sd^2 for prior 1 should be <mu*(1−mu)")}
    if (priors_sd[2]^2 > priors_mu[2] * (1 - priors_mu[2])){stop("the sd^2 for prior 0 should be <mu*(1−mu)")}

    grid_range <- lapply(1:3, function(i){
      mu <- c(post_pred[1], priors_mu)[i]
      sd <- c(post_pred[2], priors_sd)[i]

      phi   <- (mu * (1 - mu) / sd^2) - 1
      alpha <- mu * phi
      beta  <- (1 - mu) * phi

      qbeta(c(0.001,0.999), alpha, beta)})

    seq_min    <-  min(unlist(grid_range))
    seq_max    <-  max(unlist(grid_range))

    grid       <- seq(seq_min, seq_max, length.out=300000)
    dens       <- density(posterior)

    likelihood <- approx(dens$x, dens$y, rule = 2, xout = grid)$y
    phi1   <- (priors_mu[1] * (1 - priors_mu[1]) / priors_sd[1]^2) - 1
    phi2   <- (priors_mu[2] * (1 - priors_mu[2]) / priors_sd[2]^2) - 1

    prior1     <- dbeta(grid, priors_mu[1]*phi1, (1-priors_mu[1])*phi1)
    prior2     <- dbeta(grid, priors_mu[2]*phi2, (1-priors_mu[2])*phi2)

    m1 <- sum(likelihood * prior1)*(grid[2]-grid[1])
    m2 <- sum(likelihood * prior2)*(grid[2]-grid[1])

    results <- c(m1/m2, "map"=maxpost(posterior), hdi_fun(posterior, level=level))
    xlims   <- c(seq_min, seq_max)

  }else{stop("Not the correct family selected.")}

  results        <- c(plogis(log(results[1])), results)
  names(results) <- c("P(M1>M0|Data)", "BF10", "predictive map", "ll", "ul")
  results        <- round(results, 3)

  df <-  data.frame(g=rep(c("Likelihood", "Prior 1", "Prior 0"), each=length(grid)),
                    x=rep(grid,3),
                    y=c(likelihood/max(likelihood), prior1/max(prior1), prior2/max(prior2)))

  if(is.null(xlim)){xlims <- xlims}else{xlims <- xlim}

  like_fig <- suppressWarnings(ggplot(df[df$g != "Posterior",], aes(x,y,group=g,colour = g))+
                                 geom_line()+xlim(xlims)+xlab(variable)+labs(colour = NULL)+
                                 scale_color_manual(breaks = c("Likelihood", "Prior 1", "Prior 0"),
                                                    values = c("tomato3", "dodgerblue", "green2"))+
                                 theme_classic()+ylab("Scaled probablility density")+
                                 theme(legend.position = "bottom",
                                       axis.text.y = element_blank(),
                                       axis.ticks.y = element_blank()))

  return(list(summary=results, likelihood=like_fig))}
