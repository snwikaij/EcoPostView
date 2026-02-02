select_root_target <- function(object, target_value=NULL, plot_posterior=NULL,
                               epsilon=NULL, xlab="Variable", level=0.9, warning=50){

  object <- object$simulations

  if(is.null(target_value)){stop("target value is not selected")}

  if(is.null(plot_posterior)){stop("a variable for the posterior is not selected")}

  if(is.null(epsilon)){warning(paste("epsilon is not set manual therefore epsilon =" , round(sd(object$yhat)/2000,3), "make sure this is reasonable"))
    epsilon <- sd(object$yhat)/2000}

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

  #kernel for priors
  kern           <- object[,-1]
  kernels        <- lapply(kern, density)

  #select kernel
  kern_plot   <- data.frame(type="Prior", x=kernels[[plot_posterior]]$x, y=kernels[[plot_posterior]]$y)
  kern_plot$y <- kern_plot$y/max(kern_plot$y)

  #Distance/tolerance for abc selection
  object$epsilon <- object$yhat-target_value
  object         <- object[abs(object$epsilon)<epsilon,]

  #select posterior
  post_plot   <- density(post_vals <- object[,plot_posterior])
  post_plot   <- data.frame(type="Posterior", x=post_plot$x, y=post_plot$y)
  post_plot$y <- post_plot$y/max(post_plot$y)
  postprior   <- rbind(kern_plot, post_plot)

  results        <- c(maxpost(post_vals), mean(post_vals), sd(post_vals),
                      hdi_fun(post_vals, level=level),
                      length(post_vals))
  names(results) <- c("map", "mu", "se", "ll", "ul", "accepted")
  results        <- round(results, 3)

  if(length(post_vals)<warning){warning(paste0("the number of accepted prior values is <", warning, ", the posterior might be unstable"))}

  #plot prior and posterior
  pplot <- ggplot(postprior, aes(x, y, col=type))+
    geom_line()+xlab(xlab)+ylab("Scaled probability density")+
    geom_point(data=as.data.frame(t(results)), aes(x=map, y=0), inherit.aes = F, size=3)+
    geom_errorbarh(data=as.data.frame(t(results)), aes(y=0, xmin=ll, xmax=ul), height=0, inherit.aes = F)+
    scale_color_manual(breaks=c("Posterior", "Prior"),
                       values=c("tomato3", "dodgerblue"))+
    theme_classic()+labs(colour = NULL)+
    theme(legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  return(list(summary=results, figure=pplot, abc=object))}
