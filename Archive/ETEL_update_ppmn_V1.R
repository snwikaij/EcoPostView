update_ppmn <- function(object, new_data, nsim = 1000, level = 0.9,
                        focus = "local", center="median", penalize=NULL,
                        mad_constant=1.4826, covar=F, seed = 123, print_progress=F){

  #constant
  eps <- .Machine$double.eps

  #ETI intervals of posterior
  wquantile <- function(x, w, probs) {
    ok <- is.finite(x) & is.finite(w)
    x <- x[ok]; w <- w[ok]
    if (length(x) == 0) return(rep(NA_real_, length(probs)))
    o <- order(x)
    x <- x[o]; w <- w[o]
    w <- w / sum(w)
    cw <- cumsum(w)
    sapply(probs, function(p) x[which(cw >= p)[1]])}

  #structure data based on execution order
  new_data <- new_data[, colnames(new_data) %in% object$Structure$exec_dag$order, drop = FALSE]

  #keep rows with at least one finite entry
  keep     <- apply(new_data, 1, function(r) any(is.finite(as.numeric(r))))
  new_data <- new_data[keep, , drop = FALSE]

  #drop empty columns
  new_data <- new_data[,colSums(new_data, na.rm=T)>0]

  if (nrow(new_data) == 0) {stop("after filtering, new_data has 0 rows with finite values")}

  #from the new data create a matrix
  x_obs <- as.matrix(new_data)

  #Make predictions
  preds <- predict_ppmn(object, new_data, nsim = nsim)

  #select predicted nodes
  pred_nodes_all <- setdiff(names(preds$Variance), preds$Roots)

  #select observed nodes/variables
  obs_nodes <- intersect(colnames(x_obs), pred_nodes_all)

  if (length(obs_nodes) == 0) {stop("No observed (non-root) nodes in new_data match model predicted nodes")}

  ###########################
  #standardize over all nsim#
  ###########################
  #pre-extract observed vectors for speed
  obs_mat <- as.matrix(x_obs[, obs_nodes, drop = FALSE])

  #expected (median, mean or mode) value per node
  center_method <- center

  #different center possibilities
  center_fun <- switch(
    center_method,
    median = function(x){median(x, na.rm = TRUE)},
    mean = function(x){mean(x, na.rm = TRUE)},
    mode = function(x){
      x <- x[is.finite(x)]
      if (length(x) == 0) return(NA_real_)
      d <- density(x)
      d$x[which.max(d$y)]}, stop("unknown center method, needs to be, median, mean or mode"))

  #stack all residuals to derive center
  res_stack <- do.call(rbind,
                       lapply(seq_len(nsim), function(i) {
                         x_hat <- do.call(cbind, lapply(preds$Variance[obs_nodes], `[[`, i))
                         obs_mat - x_hat}))

  #derive mad from residuals
  MAD <- apply(res_stack, 2, function(x) {
    mad(x,center = median(x, na.rm = T), constant = mad_constant,na.rm = TRUE)})

  #if scale is basically 0 set some small value
  MAD[MAD==0] <- eps

  #set penalty per node
  node_penalty <- rep(1, length(obs_nodes))

  #appoint penalty to node
  node_penalty <- rep(1, length(obs_nodes))
  names(node_penalty) <- obs_nodes
  if (!is.null(penalize)){node_penalty[names(penalize)] <- penalize[names(penalize) %in% obs_nodes]}

  #detect binairy variables
  detect_binary <- function(x){ux <- unique(na.omit(x)); length(ux) <= 2 && all(ux %in% c(0, 1))}

  #detect count variables
  detect_count <- function(x){x <- na.omit(x); length(x) > 0 && all(x >= 0) && all(abs(x - round(x)) < eps)}

  tilt_fun <- function(r, sigma2 = 1, eps = 1e-12) {

    r <- as.numeric(r)
    ok <- is.finite(r)
    r  <- r[ok]

    if (length(r) == 0)
      return(rep(NA_real_, length(ok)))

    # Softmax weights with 2 tilt parameters
    softmax_w <- function(lambda) {
      a <- lambda[1] * r + lambda[2] * (r^2 - sigma2)
      a <- a - max(a)                 # stabilize
      w <- exp(a)
      w / (sum(w) + eps)
    }

    # Objective: squared moment conditions
    obj <- function(lambda) {
      w <- softmax_w(lambda)
      m1 <- sum(w * r)
      m2 <- sum(w * (r^2 - sigma2))
      m1^2 + m2^2
    }

    # Solve
    sol <- optim(c(0, 0), obj, method = "BFGS")
    lambda_star <- sol$par

    w <- softmax_w(lambda_star)

    # map back to full length
    out <- rep(0, length(ok))
    out[ok] <- w

    out / sum(out)
  }

  ###########
  #MERR dist#
  ###########
  residuals <- matrix(NA_real_, nrow = nsim, ncol = ncol(obs_mat))
  colnames(residuals) <- colnames(obs_mat)

  for (i in seq_len(nsim)) {
    if (print_progress) print(i)

    #select simulation i
    x_hat           <- do.call(cbind,lapply(preds$Variance[obs_nodes], `[[`, i))
    x_hat           <- as.matrix(x_hat)
    colnames(x_hat) <- obs_nodes

    #create empty matrix
    loss <- matrix(NA_real_, nrow = nrow(obs_mat), ncol = ncol(obs_mat))
    colnames(loss) <- obs_nodes

    #for column j
    for (j in seq_len(ncol(obs_mat))){

      y_raw  <- obs_mat[,j]
      y_hat  <- x_hat[,j]

      #detect bin
      if (detect_binary(y_raw)){

        loss[, j] <- y_raw-y_hat

        bad <- !is.finite(loss[,j])

        if (any(bad)) {
          finite_vals <- loss[is.finite(loss[,j])]
          cap         <- if (length(finite_vals)){max(finite_vals)}else{1}
          loss[bad,j] <- cap
        }

        #detect count
      }else if(detect_count(y_raw)){
        loss[, j]  <- sqrt(y_raw)-sqrt(y_hat)

        bad <- !is.finite(loss[,j])
        #set max otherwise 1
        if (any(bad)) {
          finite_vals  <- loss[is.finite(loss[, j]), j]
          cap          <- if(length(finite_vals)){max(finite_vals)}else{1}
          loss[bad,j]  <- cap}

        #detect cont
      }else{
        loss[, j]  <- (y_hat-y_raw)/MAD[j]

        bad <- !is.finite(loss[,j])

        #set to max
        if (any(bad)) {
          finite_vals <- loss[is.finite(loss[,j])]
          cap         <- if(length(finite_vals)){max(finite_vals)}else{1}
          loss[bad,j] <- cap}
      }
    }

    #penalty
    loss <- sweep(loss, 2, node_penalty, "*")

    residuals[i, ] <- colMeans(loss, na.rm = TRUE)}

  residuals <- apply(residuals, 2, function(col) {
    bad <- !is.finite(col)
    if (any(bad)) col[bad] <- max(col[is.finite(col)], na.rm = TRUE)
    col})

  colnames(residuals) <- obs_nodes
  #####################
  #derive node weights#
  #####################

  etel_w <- lapply(1:ncol(residuals), function(i){
    dist <- residuals[,i]
    tilt_fun(dist)})

  w_node           <- do.call(cbind.data.frame, etel_w)
  colnames(w_node) <- obs_nodes

  ################
  #global weights#
  ################

  logw <- rep(0, nsim)
  for (nm in names(etel_w)) {
    logw <- logw + log(etel_w[[nm]] + eps)}

  logw <- logw - max(logw)

  #prevent issues
  w_global <- exp(logw)

  #global weights
  w_global <- w_global/(sum(w_global) + eps)

  #######################################
  #summarize all the nodes using weights#
  #######################################
  node_summaries     <- list()
  node_params_weight <- list()

  for(node in names(preds$Parameters)){

    draws <- preds$Parameters[[node]]
    if(node %in% object$Roots) next

    #weight per node
    if(focus == "global"){
      w <- w_global
    }else if(focus == "local"){
      w <- w_node[[node]]
      if(is.null(w)){
        w <- rep(1/nsim, nsim)}
    }else{stop("focus needs to be global or local")}

    #appoint weights to parameter draws that are not inf
    weights   <- w/sum(w)

    #combine all parameter draws
    combi_params <- do.call(rbind, draws)
    combi_params <- as.matrix(combi_params)
    if (is.null(dim(combi_params))) next

    #provide column names
    if (is.null(colnames(combi_params))) {
      colnames(combi_params) <- paste0("b", seq_len(ncol(combi_params)) - 1)}

    edge_name <- names(object$Parameters)[ sub("~.*", "", names(object$Parameters)) == node][1]
    if (is.na(edge_name)) next

    #colnames(combi_params) <- paste0(edge_name, "_", colnames(combi_params))

    #parameter summary mu and se
    mu     <- colSums(combi_params * weights)
    diffmu <- sweep(combi_params, 2, mu, "-")
    se     <- sqrt(colSums((diffmu^2) * weights))

    #parameter intervals
    llul <- apply(combi_params,2, function(x) wquantile(x, weights, c((1-level)/2, 1-(1-level)/2)))

    node_summaries[[node]]     <- list(mu=mu, se=se, ll=llul[1, ], ul=llul[2, ])
    node_params_weight[[node]] <- combi_params}

  #set names of the summary
  node_summaries <- Filter(Negate(is.null), node_summaries)

  ###################
  #derive cor matrix#
  ###################
  #parameters for global weighting

  if(covar==T){
    params <- do.call(cbind, node_params_weight)
    params <- as.matrix(params)
    params[!is.finite(params)] <- 0

    #calculate mean
    mu_par      <- colSums(params * w_global)
    diff_mu_par <- sweep(params, 2, mu_par, "-")
    diff_mu_par[!is.finite(diff_mu_par)] <- 0

    #some jitter to prevent 0 var
    sigma <- t(diff_mu_par) %*% (diff_mu_par *  w_global)
    sigma <- sigma + diag(eps, ncol(sigma))

    rownames(sigma) <- rownames(object$Sigma)
    colnames(sigma) <- colnames(object$Sigma)

    warning("estimating the covariance over all parameter can result in degeneracy issues.")

    object$Sigma            <- sigma}

  ###################################
  #Move all back to object to update#
  ###################################
  old_params <- object$Parameters
  edge_names <- names(old_params)
  dep_of_edge <- sub("~.*", "", edge_names)

  for (e in seq_along(edge_names)) {
    dep <- dep_of_edge[e]
    if (!dep %in% names(node_summaries)) next

    summ    <- node_summaries[[dep]]
    sub_old <- old_params[[e]]

    for (b in names(sub_old)) {
      if (!b %in% names(summ$mu)) next
      sub_old[[b]]["mu"] <- as.numeric(summ$mu[[b]])
      sub_old[[b]]["se"] <- as.numeric(summ$se[[b]])
      sub_old[[b]]["ll"] <- as.numeric(summ$ll[[b]])
      sub_old[[b]]["ul"] <- as.numeric(summ$ul[[b]])
    }

    old_params[[e]] <- sub_old
  }

  object$Update           <- object$Update + 1
  object$`Old parameters` <- object$Parameters
  object$Parameters       <- old_params
  KL_div                  <- sum(w_global * log((w_global + eps) / (1/nsim)))
  object$Info             <- list(KL_div = KL_div, KL_exp = exp(-KL_div), rel_KL_div = KL_div/log(nsim))

  #############################
  #Update residual diagnostics#
  #############################
  #Make predictions with updated mod
  preds_update <- predict_ppmn(object, new_data, nsim = nsim)

  #select predicted nodes updated mod
  pred_nodes_update <- setdiff(names(preds_update$Variance), object$Roots)

  #updated model
  residuals_list <- list()

  for (i in seq_len(nsim)) {
    if (print_progress){print(i)}

    # predicted values of variance i
    x_hat <- do.call(cbind,lapply(preds_update$Variance[obs_nodes], `[[`, i))
    x_hat <- as.matrix(x_hat)
    colnames(x_hat) <- obs_nodes

    #residuals
    residuals_list[[i]] <- obs_mat - x_hat}

  #set names back to residuals
  resid_names <- colnames(residuals_list[[1]])

  #organize
  out <- lapply(resid_names, function(col) {
    do.call(cbind.data.frame, lapply(residuals_list, function(x) x[, col, drop = FALSE]))})

  #give names to mat
  names(out) <- resid_names

  #residual center
  center_resids <- lapply(out, function(x) apply(x, 1, function(x) center_fun(x)))

  #stack all residuals to derive center
  res_stack <- do.call(rbind,lapply(seq_len(nsim), function(i) {
    x_hat <- do.call(cbind, lapply(preds_update$Variance[obs_nodes], `[[`, i))
    obs_mat - x_hat}))

  #derive mad from residuals
  MAD <- apply(res_stack, 2, function(x){mad(x,center = median(x, na.rm = T), constant = mad_constant,na.rm = TRUE)})

  #extract and standardize median residuals
  full_pred_resid <- do.call(cbind, center_resids)
  full_resid      <- cbind(new_data[object$Roots], do.call(cbind, center_resids))
  stand_resid     <- sweep(full_pred_resid,  2, MAD, "/")

  #create long list for the residuals
  df_long <- data.frame(
    Row = rep(seq_len(nrow(stand_resid)), times = ncol(stand_resid)),
    Variable = rep(colnames(stand_resid), each = nrow(stand_resid)),
    Value = as.vector(stand_resid))

  #quantile limits of plot
  qy <- quantile(df_long$Value, c(.025,.975), na.rm = T)

  #residual boxplot
  rpl <- ggplot(df_long, aes(Variable, Value))+ylim(qy)+
    geom_hline(yintercept = 0, col="tomato3", lty=2)+ylab("Residuals")+
    geom_boxplot(outlier.shape = NA)+geom_jitter(alpha=0.2, pch=19, width = 0.2)+
    theme_classic()

  object$Residuals$list          <- residuals_list
  object$Residuals$residual_plot <- rpl

  return(object)
}
