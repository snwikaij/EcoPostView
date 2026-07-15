#' A generic function that defines and predicts with different functions
#'
#' @param func Name of the function
#' @param beta Parameters of the functions
#' @param trans Transformation type applied over all independent variables
#' @param x matrix of independent variables
#'
#' @export
predict_edge <- function(func, beta, trans, x){

  #transformation functions under estimation
  if (!is.na(trans)){
    if (trans == "log") {
      x[x <= 0] <- NA_real_
      x <- log(x)}
    if (trans == "sqrt") {
      x[x < 0] <- NA_real_
      x <- sqrt(x)}}

  #link function of standard glms
  if (func %in% c("identity", "log", "logit")){eta <- beta[1] + drop(x %*% beta[-1])}
  if (func == "identity"){return(eta)}
  if (func == "log"){eta <- pmax(pmin(eta, 30), -30); return(exp(eta))}
  if (func == "logit"){return(plogis(eta))}

  #non standard models that are useful
  if (func == "gompertz"){
    b0 <- beta[1]
    b1 <- beta[2]
    b2 <- beta[3]

    x <- as.matrix(x)
    x1 <- x[, 1]

    # sum(bj * xj) for j >= 3 corresponds to x columns 2..p
    if (ncol(x) > 1) {
      gamma <- drop(x[, -1, drop = FALSE] %*% beta[-(1:3)])
    } else {
      gamma <- 0
    }

    eta <- b1 + gamma - b2 * x1
    return(b0 * exp(-exp(eta)))}
  if (func == "sigmoidal"){
    b0     <- beta[1]
    b1     <- beta[2]
    b2     <- beta[3]

    x <- as.matrix(x)
    x1 <- x[, 1]

    if (ncol(x) > 1) {gamma <- drop(x[, -1, drop = FALSE] %*% beta[-(1:3)])}else{gamma <- 0}
    eta <- (x1 - b1 + gamma) / b2
    return(b0 / (1 + exp(-eta)))}
  if (func == "gaussian"){return(beta[1] * exp(-0.5 * ((x - beta[2]) / beta[3])^2))}
  if (func == "class"){p <- exp(dnorm(x, beta[1], beta[2], log = TRUE) - dnorm(x, beta[3], beta[4], log = TRUE))
  return(as.numeric((p * beta[5]) / ((p * beta[5]) + (1 - beta[5]))))}}
