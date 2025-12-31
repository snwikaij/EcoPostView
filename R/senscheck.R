#' Sensitivity check
#'
#' @param mod1 First model preferable the model fitted using clearly defined priors
#' @param mod0 An empty model only using a vague prior
#' @param interval What credibility level to display
#' @param order_predictor The order of the predictors
#' @param order_group The order of the groups
#' @param xlab The x-label text
#'
#' @importFrom ggplot2 coord_flip
#' @importFrom ggplot2 geom_density
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 scale_fill_gradientn
#' @importFrom ggplot2 geom_text
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 position_dodge
#'
#' @export
senscheck <- function(mod1, mod0, interval=0.9,
                      order_predictor=NULL, order_group=NULL,
                      xlab="Log(P(M1|Data, Info)/P(M0|Data, Info))"){

  if(!all(mod1$Chains_mu[-5]$parameter ==  mod0$Chains_mu[-5]$parameter &
          mod1$Chains_mu[-5]$predictor == mod0$Chains_mu[-5]$predictor &
          mod1$Chains_mu[-5]$link == mod0$Chains_mu[-5]$link &
          mod1$Chains_mu[-5]$group == mod0$Chains_mu[-5]$group)){stop("Not the same model for mod1 as for mod0")}
  if(!any(names(mod1$Estimates)%in%c("b0", "b1") & names(mod0$Estimates) %in%c("b0", "b1"))){stop("Currently only possible when b0 and b1 are indicated.")}

  maxpost       <- function(x){d <- density(x); d$x[which.max(d$y)]}
  hdi_fun       <- function(x, level){

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

  overlaydf     <- rbind(cbind(Model="M1", mod1$Chains_mu[-5], Est=mod1$Chains_mu$estimate), cbind(Model="M0", mod0$Chains_mu[-5], Est=mod0$Chains_mu$estimate))
  overlaydf     <- overlaydf[overlaydf$parameter != "b0",]

  sens_df_odds  <- suppressWarnings(cbind(mod1$Chains_mu[-5], odds=log(mod1$Chains_mu$estimate/mod0$Chains_mu$estimate)))
  aggr_df_odds  <- aggregate(data=sens_df_odds, odds~parameter+predictor+link+group, function(x) c(maxpost(x), mean(x, na.rm=T), sd(x, na.rm = T), hdi_fun(x, interval)))
  sens_analysis <- setNames(cbind(code=paste(aggr_df_odds[,1], aggr_df_odds[,2],
                                             aggr_df_odds[,3], aggr_df_odds[,4], sep="_"),
                                  aggr_df_odds[-5], unlist(aggr_df_odds$odds)),c("code", "parameter", "predictor", "link", "group", "map", "mu", "se", "ll", "ul"))

  sens_df_sign   <- cbind(mod1$Chains_mu, mod0$Chains_mu[,5])
  colnames(sens_df_sign)[5:6] <- c("mod1", "mod0")
  shift_function <- function(x) plogis(log(sum(x$mod1>0)/sum(x$mod0>0)))
  code           <- paste(sens_df_sign[,1], sens_df_sign[,2], sens_df_sign[,3], sens_df_sign[,4], sep="_")

  pos_df <- as.data.frame(unlist(lapply(split(sens_df_sign, code), shift_function)))
  pos_df <- data.frame(code=rownames(pos_df), post=pos_df[,1])
  sens_analysis <- merge(sens_analysis, pos_df)
  split_odds    <- split(sens_analysis, sens_analysis$parameter)

  if(!is.null(order_predictor)){
    if(length(unique(split_odds$b1$predictor)) == length(order_predictor)){
      split_odds$b1$predictor <- factor(split_odds$b1$predictor, levels = order_predictor)
      overlaydf$predictor     <- factor(overlaydf$predictor, levels = order_predictor)}else{stop("Number of predictor names to order by is not of the same length as the number of predictors in the model.")}}

  if(!is.null(order_group)){
    if(!is.null(order_group) && length(unique(split_odds$b1$group)) == length(order_group)){
      split_odds$b1$group     <- factor(split_odds$b1$group, levels = rev(order_group))
      overlaydf$group         <- factor(overlaydf$group, levels =  order_group)}else{stop("Number of group names to order by is not of the same length as the number of groups in the model.")}}

  overlaydf$Model <- factor(overlaydf$Model, levels =  c("M1", "M0"))
  overlaydf       <- split(overlaydf, overlaydf$link)

  pl1 <- lapply(overlaydf, function(x) ggplot(x, aes(Est, col=Model))+
                  ggplot2::geom_density()+xlab("Estimate")+ylab("Posterior density")+
                  ggplot2::scale_color_manual(breaks=c("M1", "M0"),
                                     values = c("grey10",  "grey70"))+
                  ggplot2::geom_vline(xintercept = 0, lty=2, lwd=0.6, col="tomato3")+
                  ggplot2::facet_wrap(.~group+predictor, scales = "free")+
                  theme_classic()+
                  theme(legend.position = "bottom",
                        axis.text.x = element_text(size = 6),
                        axis.text.y = element_blank(),
                        axis.ticks.y = element_blank()))

  pl2 <- ggplot(split_odds$b1, aes(x=group, y=map, col=link))+
    ggplot2::coord_flip()+
    ggplot2::scale_color_manual(breaks=c("log", "logit", "identity"),
                       values = c("grey10", "grey50", "grey90"))+
    ggplot2::facet_wrap(.~predictor)+
    ggplot2::geom_errorbar(data=split_odds$b1, aes(ymin=as.numeric(ll), ymax=as.numeric(ul), y=0),
                  width=0, lwd=0.6, position = ggplot2::position_dodge(width = 0.7))+
    ylab(xlab)+labs(col="Link")+
    ggplot2::geom_hline(yintercept = 0)+
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.7))+
    theme_classic()+
    theme(axis.title.y = element_blank())

  if(is.null(order_group)){lord_group <- length(unique(split_odds$b1$group))}else{lord_group <- length(order_group)}
  if(is.null(order_predictor)){lord_pred <- length(unique(split_odds$b1$predictor))}else{lord_pred <- length(order_predictor)}

  pl3 <- ggplot(split_odds$b1 , aes(x= group, y = predictor, fill = post)) +
    ggplot2::geom_tile()+labs(fill="Shift in posterior positive direction")+
    ggplot2::facet_wrap(.~link, nrow=lord_group, ncol=lord_pred)+
    ggplot2::geom_text(aes(label = round(post, 2)), color = "black")+
    ggplot2::scale_fill_gradientn(colors = c("tomato3", "white", "dodgerblue3"),
                                  breaks=seq(0,1,0.2),
                                  limits=c(0,1)) +
    theme_classic()+
    theme(axis.title = element_blank(),
          legend.position = "bottom")

  return(invisible(list(overlay=pl1, posterior_odds=pl2, posterior_shift=pl3, table=sens_analysis)))}
