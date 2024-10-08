#' Extract ABC values
#'
#' @param obj An object from the abctoz function
#' @param dist_threshold A threshold to accept the number of iterations
#' @param interval Interval level
#' @param n_dens Number of density curves to generate
#' @param xpos Values for position of text along x-axis
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
#'
#' @export
extrabc <- function(obj, dist_threshold=0.3,
                    interval=0.9, n_dens=100,
                    xpos=6, alpha_dens=0.3){

  #Extract the posterior by threshold distance
  priors           <- obj$iterations

  #Accepted simulations
  accepted <- which(abs(priors$dist)<dist_threshold)

  #Posterior accepted values
  posterior        <- priors[accepted,]

  #local linear correction function
  loc_lm_cor <- function(par, dist, link="log"){

    dfloc <- data.frame(par=par, dist=dist)

    if(link=="log"){
      p_org   <- predict(mod1 <- glm(par~dist, data=dfloc, family = gaussian(link="log")), type = "response")
      wt      <- 1/resid(lm(par~1, data=dfloc))^2
      p_wt    <-  predict(mod2 <- glm(par~dist, data=dfloc, family = gaussian(link="log"), weights = wt), type = "response")}
    else if(link=="identity"){
      p_org   <- predict(mod1 <- glm(par~dist, data=dfloc, family = gaussian(link="identity")), type = "response")
      wt      <- 1/resid(lm(par~1, data=dfloc))^2
      p_wt    <-  predict(mod2 <- glm(par~dist, data=dfloc, family = gaussian(link="identity"), weights = wt), type = "response")}

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
  mu_extracted   <- loc_lm_cor(posterior$mu, posterior$dist)
  mu_adj         <- mu_extracted$adjust
  mu_plot        <- mu_extracted$plot

  #sd
  sd_extracted   <- loc_lm_cor(posterior$sd, posterior$dist)
  sd_adj         <- sd_extracted$adjust
  sd_plot        <- sd_extracted$plot

  #cens
  cens_extracted <- loc_lm_cor(posterior$cens, posterior$dist, link = "logit")
  cens_adj       <- cens_extracted$adjust
  cens_plot      <- cens_extracted$plot

  #create n_dens density lines
  densline <- list()

  #Create a density curve for the data
  zdens <- data.frame(x=density(obj$data, bw=0.2)$x, y=density(obj$data, bw=0.2)$y)
  zdens <- zdens[zdens$x>0,]

  #vec for estimated values
  est_array <- array(NA, dim=c(nrow(posterior), 5))

  for(i in 1:nrow(posterior)){

    x        <- obj$simulations[[accepted[i]]]
    sim_dens <- density(x, bw=0.1)
    dens_df  <- data.frame(i=i, x=sim_dens$x, y=sim_dens$y)

    #set a length  for approximation of density curve
    xlen <- seq(max(c(min(zdens$x), min(dens_df$x))), min(c(max(zdens$x), max(dens_df$x))), length.out=100)

    data_val <- approx(zdens$x, zdens$y, xout = xlen)$y
    sim_val  <- approx(dens_df$x, dens_df$y, xout = xlen)$y

    est_array[i,]   <- c(cor(data_val, sim_val)^2,
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
  plhist <- ggplot(data.frame(z=obj$data), aes(z))+xlim(0, 10)+
    xlab("z-value")+
    geom_histogram(col="black", fill="grey70",
                   alpha=0.2, position="identity",
                   binwidth = 0.5,
                   boundary = 0, closed = "left")+
    xlab("z-value") +
    ylab("Density") +
    xlim(0, 10) +
    theme_classic() +
    scale_alpha(guide = 'none')

  #maximum value histogram
  max_hist   <- max(ggplot_build(plhist)$data[[1]]$ymax)

  #maximum range for histogram
  ymax_hist <- max_hist+max_hist*sd_max

  #adjustment factor
  adj_factor <- max_hist/mean_max

  #final figures
  plhist <- plhist+geom_line(data=densline, aes(x = x, y = y*adj_factor, group = as.factor(i)), alpha = alpha_dens, color = "grey30", inherit.aes = F)+
    ylim(0, ymax_hist)+annotate("text", x = xpos, y = ymax_hist/2, label = lab)

  #Plot the histogram
  pldens <- ggplot(densline, aes(x, y, group = as.factor(i)))+xlim(0, 10)+
    xlab("z-value")+annotate("text", x = xpos, y = mean_max, label = lab)+
    geom_line(data=densline, aes(x = x, y = y,
                                 group = as.factor(i)),
              alpha = alpha_dens, color = "grey80", inherit.aes = F)+
    geom_line(data=zdens, aes(x, y))+
    xlab("z-value") +
    ylab("Density") +
    xlim(0, 10) +
    theme_classic() +
    scale_alpha(guide = 'none')

  #Use selected interval level
  interval_level              <- 0.5+c(-1,1)*interval/2

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

  return(list(summary=output,
              density=pldens,
              hist=plhist,
              posterior=posterior,
              check=list(mu_plot, sd_plot, cens_plot)))}
