#' P-to-z-values transformation
#'
#' @param estimate The Effect-size
#' @param stderr Standard error
#' @param zi z-value
#' @param pi p-value
#' @param names Additional names
#' @param alpha Alpha level
#' @param absolute Display estimates on the absolute scale
#' @param xlim Limits of the x-axis
#' @param ylim Limits of the y-axis
#'
#' @importFrom stats rgamma
#' @importFrom stats rbeta
#' @importFrom stats density
#' @importFrom stats sd
#' @importFrom stats qnorm
#' @importFrom stats pnorm
#' @importFrom stats quantile
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 scale_color_manual
#' @importFrom ggplot2 theme
#'
#' @export
ptoz <- function(estimate=NULL, stderr=NULL, zi=NULL, pi=NULL, names=NULL, alpha=0.05,
                 absolute=F, xlim=NULL, ylim=NULL){

  zint <- abs(qnorm(alpha/2))

  if(!is.null(estimate)){

    if(!is.null(stderr)){
      se  <- stderr
      z   <- estimate/se
      p   <- (1-pnorm(abs(z)))*2
      ll  <- estimate-zint*se
      ul  <- estimate+zint*se}

    if(!is.null(zi)){
      z   <- zi
      se  <- abs(estimate)/abs(z)
      p   <- (1-pnorm(abs(z)))*2
      ll  <- estimate-zint*se
      ul  <- estimate+zint*se}

    if(!is.null(pi)){
      p   <- pi
      z   <- qnorm(1 - p/2)
      se  <- abs(estimate)/abs(z)
      ll  <- estimate-zint*se
      ul  <- estimate+zint*se}

    if(!is.null(names)){
      if(length(estimate)!=length(names)){stop("Names has not the same lenght as the data")}else{
        df <- data.frame(names=names, estimate=estimate, se=se, z=z, p=p, ll=ll, ul=ul)

        df2 <- data.frame(estimate=df$estimate,
                          z=df$estimate/df$se,
                          pre=1/df$se,
                          se=df$se,
                          s=names)

        df2$sig <- ifelse(abs(df2$z)>abs(zint),"sig","nonsig")

        df2     <- df2[!(is.infinite(df2$estimate) | is.infinite(df2$pre)),]
        spldf   <- split(df2, df2$s)

        rcoef             <- quantile(abs(df2$estimate), c(.005, .995))
        yl_pos            <- seq(0, ifelse(!is.null(ylim), ylim[2], rcoef[2]), length.out=50)
        yl_neg            <- seq(0, ifelse(!is.null(ylim), ylim[1], -1*rcoef[2]), length.out=50)
        sigline_pos       <- data.frame(coef=yl_pos, se=abs(yl_pos/1.96))
        sigline_neg       <- data.frame(coef=yl_neg, se=abs(yl_neg/1.96))
        sigline_pos$linvse<- 1/sigline_pos$se
        sigline_neg$linvse<- 1/sigline_neg$se

        if(!is.null(ylim)){ylim_plot <- ylim(ylim)}else{ylim_plot <- NULL}
        if(!is.null(xlim)){xlim_plot <- xlim(xlim)}else{xlim_plot <- NULL}

        if(absolute==F){
        plot1 <- ggplot(df2, aes(x=pre, y=coef, col=sig))+
          geom_point(alpha=0.4)+xlab("1/se")+ylab("Est.")+ylim_plot+xlim_plot+
          geom_line(data=sigline_pos, aes(linvse, coef), lty=2, col="tomato4")+
          geom_line(data=sigline_neg, aes(linvse, coef), lty=2, col="tomato4")+
          geom_hline(yintercept = 0, lty=3, col="tomato4")+
          scale_color_manual(breaks = c("sig", "nonsig"),
                             values = c("tomato3", "dodgerblue3"))+
          theme_classic()+
          facet_wrap(.~s, scales = "free")+
          theme(legend.position = "none")}else{

        plot1 <- ggplot(df2, aes(x=pre, y=abs(coef), col=sig))+
          geom_point(alpha=0.4)+xlab("1/se")+ylab("|Est.|")+ylim_plot+xlim_plot+
          geom_line(data=sigline_pos, aes(linvse, coef), lty=2, col="tomato4")+
          geom_hline(yintercept = 0, lty=3, col="tomato4")+
          scale_color_manual(breaks = c("sig", "nonsig"),
                             values = c("tomato3", "dodgerblue3"))+
          theme_classic()+
          facet_wrap(.~s, scales = "free")+
          theme(legend.position = "none")}

        }
    }else{df    <- data.frame(names=1:length(estimate), estimate=estimate, se=se, z=z, p=p, ll=ll, ul=ul)
          plot1 <- NULL}

  }else{

    if(!is.null(zi)){
      z   <- zi
      p   <- (1-pnorm(abs(z)))*2}

    if(!is.null(pi)){
      p   <- pi
      z   <- qnorm(1 - p/2)}

    if(!is.null(names)){
      if(length(p)!=length(names)){stop("Names has not the same lenght as the data")}else{
        df       <- data.frame(names=names, z=z, p=p)
        plot1    <- NULL
      }
    }else{df       <- data.frame(names=1:length(p), z=z, p=p)
    plot1          <- NULL}

    if(!is.null(zi) | !is.null(pi)){

    }else{stop("Either z- or p-values need to be provided")}}

  cninfnan <- function(x){!any(is.na(x), is.infinite(x), is.nan(x))}

  vec      <- rep(NA, nrow(df))
  for(i in 1:nrow(df)){
    vec[i] <- cninfnan(df$se[i])}

  df                 <- df[vec,]
  df$z[is.nan(df$z)] <- 0
  return(list(data=df, funnel_plot=plot1))}
