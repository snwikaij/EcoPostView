#' Randomized p-values to transform to z-values
#'
#' @param p A vector of p-values from literature
#' @param operator A corresponding vector if this was "smaller" (e.g., p<.05) or "bigger" (p>.1)
#' @param specific If specific is TRUE all for p<.05 a random values is simulated as p~Uniform(0.01, 0.05) for <.01 Uniform(0.001, 0.01) etc.
#' If specific is False for every operator that is "smaller' p~Uniform(0, 0.05) often a proper approach.
#' @param seed Fixed seed set by default to 123
#'
randomizep <- function(p, operator, specific=T, seed=123){

  df  <- data.frame(p=p, operator=operator)

  if(specific == F)df$p[df$operator == "smaller"] <- 0.05

  r_p <- function(p, operator){

    if(is.na(operator)){
      p <- p}else{
        if(operator == "smaller"){
          if(p == 0.05){p <- runif(1, 0.01, 0.05)}
          else if(p == 0.01){p <- runif(1, 0.001, 0.01)}
          else if(p == 0.001){p <- runif(1, 0, 0.001)}
          else if(p == 0.0001){p <- runif(1, 0, 0.0001)}
          else(p <- runif(1, 0,  p))
        }else if(operator == "bigger"){
          p <- runif(1, p, 1)}}

    return(p)}

  set.seed(seed)

  dfpz <- data.frame(p=p, z=qnorm(1 - mapply(r_p, df$p, df$operator)/2))
  return(dfpz)}
