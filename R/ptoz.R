#' Title
#'
#' @param coef The Effect-size
#' @param sei Standard error
#' @param zi z-value
#' @param pi p-value
#' @param names Additional names
#' @param alpha alpha level
#'
ptoz <- function(coef=NULL, sei=NULL, zi=NULL, pi=NULL, names=NULL, alpha=0.05){

  zint <- abs(qnorm(alpha/2))

  if(!is.null(coef)){

    if(!is.null(sei)){
      se  <- sei
      z   <- coef/se
      p   <- (1-pnorm(abs(z)))*2
      ll  <- coef-zint*se
      ul  <- coef+zint*se}

    if(!is.null(zi)){
      z   <- zi
      se  <- abs(coef)/abs(z)
      p   <- (1-pnorm(abs(z)))*2
      ll  <- coef-zint*se
      ul  <- coef+zint*se}

    if(!is.null(pi)){
      p   <- pi
      z   <- qnorm(1 - p/2)
      se  <- abs(coef)/abs(z)
      ll  <- coef-zint*se
      ul  <- coef+zint*se}

    if(!is.null(names)){
      if(length(coef)!=length(names)){stop("Names has not the same lenght as the data")}else{
        df <- data.frame(names=names, coef=coef, se=se, z=z, p=p, ll=ll, ul=ul)

        df2 <- data.frame(coef=df$coef,
                          lae=log(abs(df$coef)),
                          logaz=log(abs(df$coef/df$se)),
                          pre=1/df$se,
                          logpre=log(1/df$se),
                          se=df$se,
                          s=names)

        df2$sig <- ifelse(df2$logaz>log(zint),"sig","nonsig")

        df2     <- df2[!(is.infinite(df2$lae) | is.infinite(df2$pre)),]
        spldf   <- split(df2, df2$s)

        lint <- lapply(spldf, function(x) c(lm(lae ~ I(se^2), weight = pre, data = x)$coef[1],
                                            lm(lae~logpre, data = x)$coef[1],
                                            lm(lae~logpre, data = x)$coef[2]))

        rlogpre           <- quantile(df2$logpre, c(.005, .995))
        xlp               <- seq(rlogpre[1], rlogpre[2], length.out=10)
        b_lines           <- cbind.data.frame(do.call(rbind, lapply(lint, function(x) cbind(xlp, x[2]+x[3]*xlp))), rep(names(lint), each=50))
        colnames(b_lines) <- c("logpre", "lae", "s")

        adjline           <- setNames(cbind.data.frame(do.call(rbind, lint))[-c(2:3)], "intercept")
        adjline$slope     <- 0
        adjline$s         <- rownames(adjline)

        rlae              <- quantile(df2$lae, c(.005, .995))
        ylp               <- seq(exp(rlae[1]), exp(rlae[2]), length.out=10)
        sigline           <- data.frame(lae=log(ylp), se=ylp/1.96)
        sigline$linvse    <- log(1/sigline$se)

        plot1 <- ggplot(df2, aes(x=logpre, y=lae, col=sig))+
          geom_point(alpha=0.4)+xlab("Ln(1/se)")+ylab("Ln(|Est.|)")+
          geom_line(data=sigline, aes(linvse, lae), lty=2, col="tomato4")+
          #geom_line(data=b_lines, aes(logpre, lae), lty=1, col="tomato4")+
          geom_abline(data=adjline, aes(intercept = intercept, slope = slope), col="blue4") +
          scale_color_manual(breaks = c("sig", "nonsig"),
                             values = c("tomato3", "dodgerblue3"))+
          theme_classic()+
          facet_wrap(.~s, scales = "free")+
          theme(legend.position = "none")

        }
    }else{df    <- data.frame(names=1:length(coef), coef=coef, se=se, z=z, p=p, ll=ll, ul=ul)
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
