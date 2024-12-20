#' Extract ABC values
#'
#' @param obj An object from the abctoz function
#' @param dist_threshold A threshold to accept the number of iterations
#' @param interval Interval level
#' @param n_dens Number of density curves to generate
#' @param xpos Values for position of text along x-axis
#' @param ypos_lim A value in fraction of the maximum density that indicates limits of the y-axis
#' @param alpha_dens Transparency of the density lines
#'
#' @importFrom stats predict
#' @importFrom stats glm
#' @importFrom stats gaussian
#' @importFrom stats lm
#' @importFrom stats resid
#' @importFrom ggplot2 geom_smooth
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 xlim
#' @importFrom ggplot2 ylim
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 ggplot_build
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 scale_alpha
#'
#' @export
extrabc <- function(obj, dist_threshold=NULL,
                    interval=0.9, n_dens=100,
                    xpos=6, ypos_lim=0.99, alpha_dens=0.3){

  #Extract the posterior by threshold distance
  priors           <- obj$iterations

  seq_val <- seq(0.001, 0.05, 0.001)
  seq_df  <- array(NA, dim=c(length(seq_val),2))

  for(i in 1:nrow(seq_df)){
    seq_set <- priors[priors$distance<seq_val[i],]
    if(nrow(seq_set)==0){
      acceptance <- 0}else{
        acceptance <- nrow(seq_set)/nrow(priors)}
  seq_df[i,]       <- c(seq_val[i],acceptance)}

  seq_df           <- setNames(as.data.frame(seq_df), c("dist", "acc"))
  seq_df           <- na.omit(seq_df)

  #Accepted simulations
  if(is.null(dist_threshold)){
  dy  <- diff(seq_df$acc)/diff(seq_df$dist)
  d2y <- diff(dy)/diff(seq_df$dist[-1])

  max_curve      <- which.max(abs(d2y))
  dist_threshold <- seq_df$dist[max_curve + 1]

  accepted         <- which(abs(priors$distance)<dist_threshold)}else{
  accepted         <- which(abs(priors$distance)<dist_threshold)}

  #Posterior accepted values
  posterior        <- priors[accepted,]

  #local linear correction function
  loc_lm_cor <- function(par, dist){

    dfloc <- data.frame(par=par, dist=dist)

    p_org   <- predict(mod1 <- lm(par~dist, data=dfloc), type = "response")
    wt      <- 1/resid(lm(par~1, data=dfloc))^2
    p_wt    <- predict(mod2 <- lm(par~dist, data=dfloc, weights = wt), type = "response")

    adjust  <- dfloc$par+p_wt-p_org

    org_adj_plot <- rbind(
      data.frame(type="Adjusted", x=dfloc$dist, y=adjust),
      data.frame(type="Original", x=dfloc$dist, y=dfloc$par))

    l <- geom_smooth(method = "glm", formula="y~x", se=F, col='tomato2', lty=2, lwd=1.2)

    plot <- ggplot(org_adj_plot, aes(x, y))+
      geom_point()+l+
      facet_wrap(.~type)+xlab("Distance")+ylab("Posterior iterations")+
      theme_classic()

    return(list(adjust=adjust, plot=plot))}

  #apply loc lm correction to mu, sd and censor
  #mu
  mu_extracted   <- loc_lm_cor(posterior$mu, posterior$distance)
  mu_adj         <- mu_extracted$adjust
  mu_plot        <- mu_extracted$plot

  #sd
  sd_extracted   <- loc_lm_cor(posterior$sd, posterior$distance)
  sd_adj         <- sd_extracted$adjust
  sd_plot        <- sd_extracted$plot

  #cens
  cens_extracted <- loc_lm_cor(posterior$cens, posterior$distance)
  cens_adj       <- cens_extracted$adjust
  cens_plot      <- cens_extracted$plot

  #create n_dens density lines
  densline <- vector("list", nrow(posterior))

  if(!is.list(obj$data)){
    #Create a density curve for the data
    zdens <- data.frame(x=density(obj$data, bw=0.2)$x, y=density(obj$data, bw=0.2)$y)
    zdens <- zdens[zdens$x>0,]}else{mu_data_dens <- vector("list", nrow(posterior))}

  #vec for estimated values
  est_array <- array(NA, dim=c(nrow(posterior), 5))

  for(i in 1:nrow(posterior)){

    x        <- obj$simulations[[accepted[i]]]
    sim_dens <- density(x, bw=0.1)
    dens_df  <- data.frame(i=i, x=sim_dens$x, y=sim_dens$y)

    if(is.list(obj$data)){
      zd        <- obj$data[[accepted[i]]]
      sim_densz <- density(zd, bw=0.1)
      densz_df  <- data.frame(i=i, x=sim_densz$x, y=sim_densz$y)

      #set a length  for approximation of density curve
      xlen     <- seq(0, 10, length.out=100)
      data_val <- approx(densz_df$x, densz_df$y, xout = xlen)$y

      mu_data_dens[[i]] <- cbind(xlen, approx(densz_df$x, densz_df$y, xout = xlen)$y)
    }else{

      #set a length  for approximation of density curve
      xlen <- seq(max(c(min(zdens$x), min(dens_df$x))), min(c(max(zdens$x), max(dens_df$x))), length.out=100)
      data_val <- approx(zdens$x, zdens$y, xout = xlen)$y}

    sim_val  <- approx(dens_df$x, dens_df$y, xout = xlen)$y

    est_array[i,]   <- c(cor(data_val, sim_val, use = 'pairwise.complete.obs')^2,
                         mu_adj[i],
                         sd_adj[i],
                         cens_adj[i],
                         posterior$threshold[i])

    densline[[i]] <- dens_df}

  #Select n_dens
  if(length(densline)<n_dens) n_dens <- length(densline)

  #Create a long format
  densline <- do.call(rbind, densline[sample(1:length(densline), size=n_dens, replace = F)])
  densline <- densline[densline$x>0,]

  #create label
  mu_stat <- round(colMeans(est_array, na.rm = T),2)
  se_stat <- round(apply(est_array, 2, function(x) sd(x, na.rm = T)),2)
  lab     <- paste0(c("R-squared=", "mean(Z)=", "sd(Z)=", "Censored=", "Threshold="), mu_stat, " (SE=", se_stat, ")\n", collapse = "")

  mean_max <- mean(aggregate(data=densline, y~i, max)[,2])
  sd_max   <- sd(aggregate(data=densline, y~i, max)[,2])

  #Plot the histogram
  plhist <- ggplot(data.frame(z=obj$raw_data), aes(z))+xlim(0, 10)+
    xlab("z-value")+
    geom_vline(xintercept = 1.96, lty=2, lwd=0.6)+
    geom_histogram(col="black", fill="grey70",
                   alpha=0.2, position="identity",
                   binwidth = 0.5,
                   boundary = 0, closed = "left")+
    xlab("z-value") +
    ylab("Count") +
    theme_classic() +
    scale_alpha(guide = 'none')

  #maximum value histogram
  max_hist   <- max(ggplot_build(plhist)$data[[2]]$ymax)

  #maximum range for histogram
  ymax_hist <- max_hist+max_hist*sd_max

  #adjustment factor
  adj_factor <- max_hist/mean_max

  #final figures
  plhist <- plhist+geom_line(data=densline, aes(x = x, y = y*adj_factor, group = as.factor(i)), alpha = alpha_dens, color = "grey30", inherit.aes = F)+
    ylim(0, ymax_hist)+annotate("text", x = xpos, y = ymax_hist/2, label = lab)

  if(is.list(obj$data)){
    zmu <- do.call(rbind, mu_data_dens)
    zmu <- aggregate(data=zmu, zmu[,2]~zmu[,1], function(x) mean(x, na.rm=T))
    zmu <- setNames(zmu, c("x", "y"))}else{zmu <- zdens}

  ylimit <- quantile(densline$y,ypos_lim)
  ypos <- quantile(densline$y,.95)
  #Plot the histogram
  pldens <- ggplot(densline, aes(x, y, group = as.factor(i)))+xlim(0, 10)+
    xlab("z-value")+annotate("text", x = xpos, y = ypos, label = lab)+ylim(0,ylimit)+
    geom_vline(xintercept = 1.96, lty=2, lwd=0.6)+
    geom_line(data=densline, aes(x = x, y = y,
                                 group = as.factor(i)),
              alpha = alpha_dens, color = "grey80", inherit.aes = F)+
    geom_line(data=zmu, aes(x, y))+
    xlab("z-value") +
    ylab("Density") +
    theme_classic() +
    scale_alpha(guide = 'none')

  #Use selected ETI interval level (I specifically choose ETI) due to its simplicity
  #Also, I do not know how HDI behaves exactly with limited posterior samples
  interval_level              <- 0.5+c(-1,1)*interval/2

  posth0h1 <- table(posterior$`H0/H1`)
  if(!length(posth0h1) == 2){
  if(names(posth0h1) == "0"){posth0h1 <- cbind("0"=posth0h1, "1"=0)}
  if(names(posth0h1) == "1"){posth0h1 <- cbind("1"=0, "1"=posth0h1)}}

  #summarize support
  supstat <- round(c(setNames(posth0h1,c("Uncensored", "Censored")),
                  setNames(posth0h1[2]/posth0h1[1],c("BayesFactor(H1/H0)")),
                  setNames(plogis(log(posth0h1[2]/posth0h1[1])),c("Probability(H1/H0)")),
                  tolerance=dist_threshold),2)

  #Generate a simple summary
  output <- data.frame(
    Statistic = c("c", "mu(z)", "sd(z)"),
    Mean = round(c(mean(cens_adj), mean(mu_adj), mean(sd_adj)), 4),
    SE = round(c(sd(cens_adj), sd(mu_adj), sd(sd_adj)), 4),
    ll = round(c(quantile(cens_adj, interval_level[1]),
                 quantile(mu_adj, interval_level[1]),
                 quantile(sd_adj, interval_level[1])), 4),
    ul = round(c(quantile(cens_adj, interval_level[2]),
                 quantile(mu_adj, interval_level[2]),
                 quantile(sd_adj, interval_level[2])), 4))

  print(list('Supportive statistics'=supstat,
            'Summary statistics'=output))

  return(invisible(list(summary=output,
              distance=seq_df,
              density=pldens,
              hist=plhist,
              posterior=posterior,
              check=list(mu_plot, sd_plot, cens_plot))))}
