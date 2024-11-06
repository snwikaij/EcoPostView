#' Title
#'
#' @param mod1 First model preferable the model fitted using clearly defined priors
#' @param mod0 An empty model only using a vague prior
#' @param interval What credibility level to display
#' @param xlab The x-label text
#'
#' @export
senscheck <- function(mod1, mod0,
                      interval=0.9,
                      xlab="Log(P(M1|Data, Info)/P(M0|Data, Info))"){

  if(!all(mod1$Chains_mu[-5]$parameter ==  mod0$Chains_mu[-5]$parameter &
          mod1$Chains_mu[-5]$predictor == mod0$Chains_mu[-5]$predictor &
          mod1$Chains_mu[-5]$link == mod0$Chains_mu[-5]$link &
          mod1$Chains_mu[-5]$group == mod0$Chains_mu[-5]$group)){stop("Not the same model for mod1 as for mod0")}
  if(!any(names(mod1$Estimates)%in%c("b0", "b1") & names(mod0$Estimates)%in%c("b0", "b1"))){stop("Currently only possible when b0 and b1 are indicated.")}

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

sens_df_odds  <- cbind(mod1$Chains_mu[-5], odds=log(mod1$Chains_mu$estimate/mod0$Chains_mu$estimate))
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

pl1 <- ggplot(split_odds$b1, aes(x=reorder(group, -map), y=map, col=link))+
  coord_flip()+
  scale_color_manual(breaks=c("log", "logit", "identity"),
                     values = c("grey10", "grey50", "grey90"))+
  facet_wrap(.~predictor)+
  geom_errorbar(data=split_odds$b1, aes(ymin=as.numeric(ll), ymax=as.numeric(ul), y=0),
                width=0, lwd=0.6, position = position_dodge(width = 0.7))+
  ylab(xlab)+labs(col="Link")+
  geom_hline(yintercept = 0)+
  geom_point(position = position_dodge(width = 0.7))+
  theme_classic()+
  theme(axis.title.y = element_blank())

pl2 <- ggplot(split_odds$b1 , aes(x= group, y = predictor, fill = post)) +
  ggplot2::geom_tile()+labs(fill="Shift in posterior positive direction")+
  facet_wrap(.~link)+
  geom_text(aes(label = round(post, 2)), color = "black")+
  ggplot2::scale_fill_gradientn(colors = c("tomato3", "white", "dodgerblue3"),
                                breaks=seq(0,1,0.2),
                                limits=c(0,1)) +
  theme_classic()+
  theme(axis.title = element_blank(),
        legend.position = "bottom")

return(invisible(list(posterior_odds=pl1, posterior_shift=pl2, table=sens_analysis)))}
