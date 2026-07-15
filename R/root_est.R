#' Root estimation function
#'
#' @param object A foundational or updated PPMN
#' @param child The name of the child vertex observed in the network
#' @param target_value The value of the child vertex observed in the network
#' @param prior A data frame of value describing the prior values in the
#' root vertices
#' @param level Credibility level (default level = 0.9)
#' @param plot_lab Plots the point estimate (median) and standard error (se)
#' in the figure
#' @param lab_x Position of the label in percentage of the x-axis
#' @param lab_y Position of the label in percentage of the y-axis
#' @param round_lab Rounding the numbers in the label
#' @param lab_size Size of the label
#' @param kde_method The type of kernel density method (default kde_method = Gaussian)
#' can be Gaussian or Gamma
#' @param kde_adjust Adjustment value of the kernel bandwidth (default kde_adjust = 1)
#' @param kl_frac The strength of the update can be derived from the Kullenback-Leiber
#' divergence applied over the number of simulations (default kl_frac = 0.9) smaller
#' values mean less precision and broader intervals
#'
#' @export
root_est <- function(object, child = NULL, target_value = NULL, prior = NULL, level = 0.9,
                     plot_lab = F, lab_x=0.75, lab_y=0.75, round_lab=1, lab_size=3.5,
                     kde_method = "Gaussian", kde_adjust = 1, kl_frac = 0.9){

  if(level>0.99999){stop("level cannot be larger than .99999")}
  if(level<0.00001){stop("level cannot be smaller than .00001")}
  if(kl_frac>0.99){stop("update fraction cannot be larger than .99")}
  if(kl_frac<0.01){stop("update fraction cannot be smaller than .01")}
  if(is.null(child)) stop("child is not provided")
  if(is.null(prior)) stop("prior data frame is not provided")
  if(is.null(target_value)) stop("target_value is not provided")
  if(!child %in% igraph::V(object$Structure$visual_dag$graph)$name) {stop("child name not found in graph")}
  if(!all(object$Roots %in% colnames(prior))){stop("not all root nodes are provided in prior")}

  #update kl_target
  kl_target <- log(1/(1-kl_frac))

  #constant
  eps <- .Machine$double.eps

  #median
  w_median <- function(x, w) {
    ok <- is.finite(x) & is.finite(w)
    x <- x[ok]; w <- w[ok]
    o <- order(x)
    x <- x[o]; w <- w[o]
    cw <- cumsum(w / sum(w))
    x[which(cw >= 0.5)[1]]}

  #eti intervals
  w_quantile <- function(x, w, probs) {
    ok <- is.finite(x) & is.finite(w)
    x <- x[ok]; w <- w[ok]
    o <- order(x)
    x <- x[o]; w <- w[o]
    w <- w / sum(w)
    cw <- cumsum(w)
    sapply(probs, function(p) x[which(cw >= p)[1]])}

  #hdi intervals
  w_hdi <- function(x, w, level) {

    ok <- is.finite(x) & is.finite(w)
    x <- x[ok]; w <- w[ok]

    # order by x
    o <- order(x)
    x <- x[o]; w <- w[o]

    w <- w / sum(w)
    cw <- cumsum(w)

    target_mass <- level

    n <- length(x)
    best_width <- Inf
    ll <- NA_real_
    ul <- NA_real_

    for (i in seq_len(n)) {
      # find the smallest j such that mass >= target
      start_mass <- if (i == 1) 0 else cw[i - 1]
      j <- which(cw >= start_mass + target_mass)[1]
      if (is.na(j)) break

      width <- x[j] - x[i]
      if (width < best_width) {
        best_width <- width
        ll <- x[i]
        ul <- x[j]
      }
    }

    c(ll = ll, ul = ul)
  }

  #density curves
  w_density <- function(x, w, adjust = kde_adjust, eps,
                        kde_method = c("Gaussian", "Gamma"),
                        n_grid = 512) {

    ok <- is.finite(x) & is.finite(w)
    x  <- x[ok]
    w  <- w[ok]
    w  <- w / sum(w)

    if (kde_method == "Gamma") {

      if (any(x <= 0)) {
        warning("Gamma KDE requires positive x. Falling back to Gaussian KDE.")
        kde_method <- "gaussian"
      } else {

        h <- stats::bw.nrd0(x) * adjust
        h <- max(h, eps)

        grid <- seq(max(eps, min(x) * 0.5), max(x) * 1.2, length.out = n_grid)

        y <- sapply(grid, function(g) {
          sum(w * stats::dgamma(
            g,
            shape = x / h + 1,
            scale = h
          ))
        })

        return(data.frame(
          x = grid,
          y = y / max(y)
        ))
      }
    } else {

      d <- suppressWarnings(stats::density(
        x,
        weights = w,
        adjust = adjust
      ))

      return(data.frame(x = d$x, y = d$y / max(d$y)))
    }
  }

  #KL optimizer function
  kl_optim <- function(lambda, loss, kl_target){

    S       <- length(loss)
    eps     <- .Machine$double.eps
    a       <- -lambda*loss
    a       <- a-max(a)
    w1      <- exp(a)
    weights <- w1/sum(w1)

    sum(weights*log(pmax(weights, eps)))+log(S)-kl_target}

  #KL stop function
  lambda_kl <- function(loss, kl_target, upper = 1000){

    f0    <- kl_optim(0, loss, kl_target)
    f_up  <- kl_optim(upper, loss, kl_target)

    if (!is.finite(f_up) || f_up<0){return(upper)}

    uniroot(kl_optim, interval=c(0, upper), loss=loss, kl_target=kl_target)$root}

  #Predictions
  sims <- do.call(rbind, lapply(seq_len(nrow(prior)), function(i) {
    out <- predict_ppmn(object, nsim = 1, new_data = prior[i, , drop = FALSE])
    cbind(yhat = unlist(out$Variance[[child]]), prior[i, , drop = FALSE])}))

  #loss functions
  resid                    <- sims$yhat - target_value
  resid[!is.finite(resid)] <- NA_real_

  #calculated MAD
  MAD <- max(mad(resid, center = median(resid, na.rm = T), na.rm = T), eps)

  #standardize residuals
  z    <- resid/MAD
  loss <- z^2

  #safety cap
  loss[!is.finite(loss)] <- max(loss[is.finite(loss)], na.rm = TRUE)

  #lambda for particles
  lambda_hat <- lambda_kl(loss, kl_target)

  #derive weights
  w <- exp(-lambda_hat * loss)
  w <- w/sum(w)

  #plot and summarise per root
  roots <- colnames(sims)[colnames(sims) != "yhat"]

  #create storage
  plots   <- list()
  summary <- list()

  for (r in roots) {

    post_vals <- sims[[r]]

    mu   <- sum(w * post_vals)
    med  <- w_median(post_vals, w)
    se   <- sqrt(sum(w * (post_vals - mu)^2))
    llul <- w_hdi(post_vals, w, level)

    prior_den      <- w_density(post_vals, rep(1 / length(w), length(w)), kde_method=kde_method, adjust = kde_adjust, eps=eps)
    prior_den$type <- "Prior"

    post_den       <- w_density(post_vals, w,  kde_method=kde_method, adjust = kde_adjust, eps=eps)
    post_den$type  <- "Posterior"

    summary[[r]]   <- data.frame(root = r, mu = mu, med = med, se = se, ll= llul[1], ul = llul[2])

    lab     <- paste0("Median = ", round(med, round_lab), " SE = ", round(se, round_lab))

    if(plot_lab == T){
      plot_lab_txt <- ggplot2::annotate("text", x=quantile(plot_df$x, lab_x), y=lab_y, size=lab_size, label=lab)
    } else {
      plot_lab_txt <- NULL
    }

    plot_df <- rbind(prior_den, post_den)

    plots[[r]] <-ggplot2::ggplot(plot_df,
                                 ggplot2::aes(x, y, col = type)) +
      plot_lab_txt+
      ggplot2::geom_line(linewidth = 0.9) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = r, y = "Scaled probability density", col = NULL) +
      ggplot2::scale_color_manual(values = c("Posterior" = "tomato3", "Prior" = "dodgerblue")) +
      ggplot2::theme(
        legend.position = "bottom",
        axis.title.y = ggplot2::element_text(size = 12),
        axis.title.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank())+
      ggplot2::geom_point(
        data = summary[[r]],
        ggplot2::aes(x = mu, y = 0),
        inherit.aes = FALSE,
        size = 2)+
      ggplot2::geom_errorbarh(
        data = summary[[r]],
        ggplot2::aes(y = 0, xmin = ll, xmax = ul),
        inherit.aes = FALSE,
        width = 0)}

  list(
    summary = do.call(rbind, summary),
    plots   = plots,
    sims    = sims,
    weights = w)}
