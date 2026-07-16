#' Plots the marginal expectations of each edge-function with credibility intervals
#'
#' @param object A foundational or updated PPMN
#' @param response The dependent (response) variable to display
#' @param data Data that contains x and y coordinates (points) to plot
#' @param level Credibility level (default level = 0.9)
#' @param nsim Number of simulations
#' @param n_grid Grid size of the added variable plot
#' @param ylim Limit of the y-axis
#' @param ylab Title of the y-axis
#' @param plot_lab Plot the performance of marginal function on the introduced data
#' @param lab_x Position of the label on the x-axis
#' @param lab_y Position of the label on the y-axis
#' @param round_lab Rounding the numbers in the label
#' @param lab_size Size of the label
#' @param lwd_exp Line width of the expected value
#' @param lwd_col Colour of the expected value
#' @param size_title Title size of the y-axis label
#' @param size_text Text size of the title
#' @param int_col Colour of the credibility intervals around the expected value
#' @param int_alpha Transparency of the credibility intervals around the expected value
#' @param pt_alpha Transparency of the points in the plot
#'
#' @description
#' This is essentially an added variable plot. It displays the marginal relation of the
#' dependent variable and all independent variables for a single edge-function. It can also
#' display a label with performance metrics on the provided data. They only serve a check
#' after updating, and do not serve as actual metrics, the default is therefore FALSE. To
#' display only the line and intervals when no data is provided a dataset can be simulate of
#' the dependent and independent variables. The pt_alpha can consequently be set to 0 so the
#' points are not displayed.
#'
#' @export
marginal_ppmn <- function(object, response, data,
                          level = 0.9, nsim = 1000, n_grid = 100,
                          ylim=NULL, ylab=NULL, plot_lab = F, lab_x=0.6, lab_y=0.85, round_lab=1,
                          lab_size=3.5, lwd_exp = 0.8, lwd_col="dodgerblue4", size_title=8, size_text=8,
                          int_col = "dodgerblue3", int_alpha=0.3, pt_alpha = 0.25){

  g          <- object$Structure$exec_dag$graph

  # draw parameters once
  mu_global        <- unlist(lapply(object$Parameters, function(e) sapply(e, `[[`, "mu")))
  names(mu_global) <- rownames(object$Sigma)

  if (all(is.infinite(unlist(lapply(object$Parameters, function(e) sapply(e, `[[`, "a")))))) {
    theta <- MASS::mvrnorm(nsim, mu_global, object$Sigma)
  } else {
    lower <- unlist(lapply(object$Parameters, function(e) sapply(e, `[[`, "a")))
    upper <- unlist(lapply(object$Parameters, function(e) sapply(e, `[[`, "b")))
    theta <- tmvtnorm::rtmvnorm(nsim, mu_global, object$Sigma, lower, upper)}

  theta <- matrix(theta,
                  nrow = nsim,
                  ncol = length(mu_global),
                  dimnames = list(seq_len(nsim), names(mu_global)))

  edge_name <- grep(paste0("^", response, "~"), names(object$Parameters), value = T)
  if(length(edge_name) == 0L){stop("No edge-function found for response: ", response)}
  if(length(edge_name) > 1L){stop("Multiple edge-functions found for response ",response,": ", paste(edge_name, collapse = ", "))}
  parents   <- strsplit(strsplit(edge_name, "~")[[1]][2], "\\+")[[1]]
  edge_pars <- object$Parameters[[edge_name]]

  if (length(parents) == 0){stop("response has no parents")}
  stopifnot(setequal(parents,igraph::V(g)$name[igraph::neighbors(g, response, mode = "in")]))

  param_idx <- paste0(edge_name, "_", names(edge_pars))
  param_idx <- param_idx[param_idx %in% colnames(theta)]

  if (length(param_idx) == 0){stop(paste("No parameters found for response:", response))}

  func  <- object$Structure$predict_dag$functions[[edge_name]]
  trans <- object$Structure$predict_dag$transformations[[edge_name]]

  plots    <- vector("list", length(parents))
  pred_exp <- predict_ppmn(object, data, nsim=1)$Expected

  get_parent_values <- function(p, data, pred_exp, trans = NA) {

    #choose source
    if (!p %in% colnames(data) || !is.numeric(data[[p]])) {
      x <- pred_exp[[p]]
    } else {
      x <- data[[p]]
    }

    #apply transformation
    if (!is.na(trans) && trans == "log") {
      x[x <= 0] <- NA_real_
      x <- log(x)}

    if (!is.na(trans) && trans == "sqrt") {
      x[x < 0] <- NA_real_
      x <- sqrt(x)}

    x[!is.finite(x)] <- NA_real_
    x }

  for(i in seq_along(parents)){

    vary   <- parents[i]
    others <- setdiff(parents, vary)

    if(!is.na(trans) && trans == "log"){

      if(!any(colnames(data) %in% vary)){
        grid <- seq(
          min(pred_exp[[vary]][pred_exp[[vary]] > 0], na.rm = TRUE),
          max(pred_exp[[vary]], na.rm = TRUE),
          length.out = n_grid)
      }else{
        grid <- seq(
          min(data[[vary]][data[[vary]] > 0], na.rm = TRUE),
          max(data[[vary]], na.rm = TRUE),
          length.out = n_grid)}
    }else{
      if(!any(colnames(data) %in% vary)){
        grid <- seq(
          min(pred_exp[[vary]], na.rm = TRUE),
          max(pred_exp[[vary]], na.rm = TRUE),
          length.out = n_grid)
      }else{
        grid <- seq(
          min(data[[vary]], na.rm = TRUE),
          max(data[[vary]], na.rm = TRUE),
          length.out = n_grid)}}

    X <- do.call(cbind, lapply(parents, function(p) {

      if(p == vary){
        grid
      }else{
        xc <- get_parent_values(p, data, pred_exp, trans)

        if(!is.na(trans) && trans == "log"){
          rep(exp(median(xc, na.rm = TRUE)), n_grid)
        }else if (!is.na(trans) && trans == "sqrt"){
          rep((median(xc, na.rm = TRUE))^2, n_grid)
        }else {
          rep(median(xc, na.rm = TRUE), n_grid)
        }}}))

    Y <- matrix(NA_real_, n_grid, nsim)

    for (s in seq_len(nsim)) {
      beta <- theta[s, param_idx]
      y <- predict_edge(func, beta, trans, X)
      y[!is.finite(y)] <- NA_real_
      #if (!is.null(object$Sigma_y) && response %in% names(object$Sigma_y)){
      #y <- y + rnorm(length(y), 0, object$Sigma_y[[response]])}
      Y[, s] <- y
    }

    df <- data.frame(
      x  = grid,
      mu = apply(Y, 1, median, na.rm = TRUE),
      ll = apply(Y, 1, quantile, probs = (1 - level)/2, na.rm = TRUE),
      ul = apply(Y, 1, quantile, probs = 1 - (1 - level)/2, na.rm = TRUE))

    if(!is.null(ylim)){
      df$ll     <- ifelse(ylim[1] > df$ll, ylim[1], df$ll)
      df$ul     <- ifelse(ylim[2] < df$ul, ylim[2], df$ul)
      ylim_func <- ylim(ylim)}else{ylim_func <- NULL}

    if(!any(colnames(data) %in% vary)){
      x_obs <- pred_exp[[vary]]
    }else{
      x_obs <- data[[vary]]}
    y_obs <- data[[response]]

    #buildmatrix at observed x, others fixed at mean
    X_obs <- do.call(cbind, lapply(parents, function(p) {

      if (p == vary) {
        x_obs
      } else {
        xc <- get_parent_values(p, data, pred_exp, trans)
        n  <- length(x_obs)

        if (!is.na(trans) && trans == "log") {
          rep(exp(median(xc, na.rm = TRUE)), n)
        } else if (!is.na(trans) && trans == "sqrt") {
          rep((median(xc, na.rm = TRUE))^2, n)
        } else {
          rep(median(xc, na.rm = TRUE), n)
        }
      }
    }))


    #posterior predictive mean at observed x
    y_hat <- rowMeans(
      sapply(seq_len(nsim), function(s) {
        beta <- theta[s, param_idx]
        yh <- predict_edge(func, beta, trans, X_obs)
        yh[!is.finite(yh)] <- NA_real_
        #if (!is.null(object$Sigma_y) && response %in% names(object$Sigma_y)){
        #  yh <- yh + rnorm(length(yh), 0, object$Sigma_y[[response]])}
        yh
      }),
      na.rm = TRUE
    )

    res <- y_obs - y_hat
    res <- res[is.finite(res)]

    MAE <- mean(abs(res))
    MD  <- median(res)
    ME  <- mean(res)

    label <- paste0(
      "Marginal residual \n",
      "MAE = ", round(MAE, round_lab), "\n",
      "Median = ", round(MD, round_lab), "\n",
      "Mean = ", round(ME, round_lab) )

    p <- ggplot(df, aes(x, mu)) +
      geom_line(lwd = lwd_exp, col = lwd_col) + ylim_func +
      geom_ribbon(aes(ymin = ll, ymax = ul), fill = int_col, alpha = int_alpha) +
      xlab(vary)+ylab(response) +
      theme_classic()+
      theme(axis.title = element_text(size=size_title),
            axis.text = element_text(size=size_text))+
      if(!any(colnames(data) %in% vary)){NULL}else{
        geom_point(data = data,
                   aes_string(x = vary, y = response),
                   inherit.aes = FALSE, pch = 19, alpha = pt_alpha)}

    ggplot_build(p)

    if (plot_lab) {

      pb <- ggplot_build(p)

      x_rng <- pb$layout$panel_params[[1]]$x.range
      y_rng <- pb$layout$panel_params[[1]]$y.range

      lab_x_pos <- x_rng[1] + diff(x_rng) * lab_x
      lab_y_pos <- y_rng[1] + diff(y_rng) * lab_y

      p <- p +
        if(!any(colnames(data) %in% vary)){NULL}else{
          annotate(
            "text",
            x = lab_x_pos,
            y = lab_y_pos,
            label = label,
            hjust = 0,
            vjust = 1,
            size = lab_size
          )}
    }

    plots[[i]] <- p
  }

  plots
}
