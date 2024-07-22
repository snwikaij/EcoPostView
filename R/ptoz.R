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
      se  <- coef/z
      p   <- (1-pnorm(abs(z)))*2
      ll  <- coef-zint*se
      ul  <- coef+zint*se}

    if(!is.null(pi)){
      p   <- pi
      z   <- qnorm(1 - p/2)
      se  <- coef/z
      ll  <- coef-zint*(coef/z)
      ul  <- coef+zint*(coef/z)}

    if(!is.null(names)){
      if(length(coef)!=length(names)){stop("Names has not the same lenght as the data")}else{
        df <- data.frame(names=names, coef=coef, se=se, z=z, p=p, ll=ll, ul=ul)
      }
    }else{df <- data.frame(names=1:length(coef), coef=coef, se=se, z=z, p=p, ll=ll, ul=ul)}

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
      }
    }else{df       <- data.frame(names=1:length(p), z=z, p=p)}

    if(!is.null(zi) | !is.null(pi)){

    }else{stop("Either z- or p-values need to be provided")}}

  cninfnan <- function(x){!any(is.na(x), is.infinite(x), is.nan(x))}

  vec      <- rep(NA, nrow(df))
  for(i in 1:nrow(df)){
    vec[i] <- cninfnan(df$se[i])}

  df                 <- df[vec,]
  df$z[is.nan(df$z)] <- 0
  df}
