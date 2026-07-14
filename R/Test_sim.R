sim_once <- function(n = 20, nsim = 3000, max_lambda = 1000,
                     lambda_grid = seq(0, 1000, length.out = 200)){

  eps <- .Machine$double.eps

  empty_out <- c(
    n = n, lambda = NA_real_, score = NA_real_,
    b0_glm = NA_real_, b0_post = NA_real_,
    b1_glm = NA_real_, b1_post = NA_real_,
    se_b0_glm = NA_real_, se_b0_post = NA_real_,
    se_b1_glm = NA_real_, se_b1_post = NA_real_
  )

  # true model
  b0_true <- 4
  b1_true <- -0.5
  sigma_y <- 40

  # training data
  x  <- runif(n, -2, 2)
  mu <- exp(b0_true + b1_true * x)

  shape <- mu^2 / sigma_y^2
  rate  <- mu / sigma_y^2
  y     <- rgamma(n, shape = shape, rate = rate)
  y     <- pmax(y, eps)

  # prior particles
  b0 <- rnorm(nsim, 4, 3)
  b1 <- rnorm(nsim, -0.5, 2)

  # particle predictions
  eta_hat <- sapply(seq_len(nsim), function(s) b0[s] + b1[s] * x)
  eta_hat <- pmin(pmax(eta_hat, -30), 30)
  yhat    <- exp(eta_hat)

  data_resid <- sweep(yhat, 1, y, "-")

  # global scaling
  scale_data     <- mad(as.vector(data_resid), na.rm = TRUE)
  data_resid     <- data_resid/max(scale_data, eps)

  # particle loss
  loss  <- colMeans(log(cosh(data_resid)), na.rm = T)

  ok <- is.finite(loss) & is.finite(b0) & is.finite(b1)

  if(sum(ok) < 2) return(empty_out)

  loss <- loss[ok]
  b0   <- b0[ok]
  b1   <- b1[ok]
  yhat <- yhat[, ok, drop = FALSE]
  eta_hat <- eta_hat[, ok, drop = FALSE]

  # SafeBayes-style predictive lambda selection
  calibration_score <- function(lambda){

    a <- -lambda * loss
    a <- a - max(a)

    w <- exp(a)
    w <- w / sum(w)

    y_pred <- as.numeric(yhat %*% w)

    se_data <- sd(y - y_pred) / sqrt(n)

    sd_post <- mean(
      apply(yhat, 1, function(v){
        mu <- sum(w * v)
        sqrt(sum(w * (v - mu)^2))
      })
    )

    ESS <- 1 / sum(w^2)

    se_post <- sd_post / sqrt(ESS)

    abs(se_post - se_data)
  }

  scores <- sapply(lambda_grid, calibration_score)
  scores[!is.finite(scores)] <- Inf

  if(all(!is.finite(scores))){return(empty_out)}

  lambda <- lambda_grid[which.min(scores)]
  score  <- min(scores, na.rm = TRUE)

  # posterior weights
  a <- -lambda * loss
  a <- a - max(a, na.rm = TRUE)

  w <- exp(a)
  w <- w / sum(w)

  # updated parameters
  b0_post <- sum(w * b0)
  b1_post <- sum(w * b1)

  # posterior uncertainty
  se_b0_post <- sqrt(sum(w * (b0 - b0_post)^2))
  se_b1_post <- sqrt(sum(w * (b1 - b1_post)^2))

  # GLM baseline
  fit_glm <- tryCatch(
    glm(
      y ~ x,
      data = data.frame(x = x, y = y),
      family = Gamma(link = "log"),
      start = c(log(mean(y)), 0),
      control = glm.control(maxit = 100, epsilon = 1e-8)
    ),
    error = function(e) NULL
  )

  if(is.null(fit_glm) || !fit_glm$converged) return(empty_out)

  glm_se <- summary(fit_glm)$coef[, 2]

  c(
    n = n,
    lambda = lambda,
    score = score,

    b0_glm = unname(coef(fit_glm)[1]),
    b0_post = b0_post,
    b1_glm = unname(coef(fit_glm)[2]),
    b1_post = b1_post,

    se_b0_glm = unname(glm_se[1]),
    se_b0_post = se_b0_post,
    se_b1_glm = unname(glm_se[2]),
    se_b1_post = se_b1_post
  )
}

ns <- c(10,30,60,120,300)

out <- do.call(rbind, lapply(ns, function(n){
  t(replicate(40, sim_once(n = n), simplify = "matrix"))}))

colnames(out)

testout <- aggregate(
  out[, c("b0_post", "b0_glm", "b1_post",
          "b1_glm","se_b0_post", "se_b0_glm",
          "se_b1_post", "se_b1_glm","lambda","score")],
  by = list(n = out[, "n"]),
  mean)

testout

pl1 <- ggplot(testout, aes(b1_post, b1_glm, size=c(log10(n))))+
  geom_vline(xintercept = -0.5, col="tomato3", lwd=.6, lty=2)+
  geom_hline(yintercept = -0.5, col="tomato3", lwd=.6, lty=2)+
  geom_point()+theme_classic()+xlab("b1 Gibbs posterior")+labs(size=NULL)+
  ylab("b1 from GLM fit")+
  theme(legend.position = "none")

pl2 <- ggplot(testout, aes(se_b1_post, se_b1_glm, size=c(log10(n))))+
  geom_abline(intercept = 0, slope = 1, col="tomato3", lwd=0.6, lty=2)+
  xlim(0, max(testout$se_b1_post))+
  ylim(0, max(testout$se_b1_glm))+
  geom_point()+
  theme_classic()+xlab("SE(b1) Gibbs posterior")+
  ylab("SE(b1) from GLM fit")+labs(size="Sample size")

pl12 <- cowplot::plot_grid(pl1, pl2)

xl      <- seq(-2, 2, 0.1)
yl_glm  <- exp(testout$b1_glm[1]*xl+testout$b0_glm[1])
yl_post <- exp(testout$b1_post[1]*xl+testout$b0_post[1])

x  <- runif(300, -2, 2)
mu <- exp(-0.4*x+4)

y     <- rgamma(length(x), mu^2/40^2, mu/40^2)
y     <- pmax(y, .Machine$double.eps)
line_df <- data.frame(xl, yl_glm, yl_post)

pl3 <- ggplot(data = data.frame(x, y), aes(x, y))+
  geom_point(pch=19, size=3, alpha=0.3)+theme_classic()+
  geom_line(data=line_df, aes(xl, yl_glm), lwd=1.2, col="dodgerblue3")+
  geom_line(data=line_df, aes(xl, yl_post), lwd=1.2, col="tomato3")

cowplot::plot_grid(pl12, pl3, ncol=1, rel_heights = c(0.3, 0.7))
testout
