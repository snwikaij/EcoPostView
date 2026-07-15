#' Update PPMN function
#'
#' @param object A foundational or updated PPMN
#' @param new_data New dataset with all root nodes and least one child node present
#' @param nsim Number of simulations needed for updating
#' @param level Credibility level (default level = 0.9)
#' @param kl_frac The strength of the update can be derived from the Kullenback-Leiber divergence applied over the number of simulations (default kl_frac = NULL)
#' @param mad_constant the constant that scale the Median Absolute Deviation (default mda_constant = 1)
#' @param seed Seed value 123
#' @param max_lambda Maximum Lambda that can be set for Generalized Bayesian Updating (default max_lambda = 1000)
#' @param covar Covariance matrix that is calculated between model parameters (default covar = FALSE)
#' @param print_progress Progress report during updating
#'
#' @export
update_ppmn <- function(object, new_data, nsim = 3000, level = 0.9,
                        kl_frac = NULL, mad_constant = 1, seed = 123,
                        max_lambda=1000, covar = F, print_progress = F){

  if(nsim<1){stop("number of simulations cannot be smaller than 1")}
  if(level>0.99999){stop("level cannot be larger than .99999")}
  if(level<0.00001){stop("level cannot be smaller than .00001")}
  if(!is.null(kl_frac)){
    if(kl_frac>1){stop("update fraction cannot be larger than 1")}
    if(kl_frac<0){stop("update fraction cannot be smaller than 0")}}
  if(mad_constant<0.01){stop("mad_constant cannot be smaller than 0.01")}

  #constant
  eps <- .Machine$double.eps

  #ETI intervals of posterior
  wquantile <- function(x, w, probs){

    ok <- is.finite(x) & is.finite(w)
    x <- x[ok]; w <- w[ok]

    if(length(x) == 0){return(rep(NA_real_, length(probs)))}

    o <- order(x)
    x <- x[o]; w <-w[o]
    w <- w/sum(w)
    cw <- cumsum(w)

    sapply(probs, function(p) x[which(cw >= p)[1]])}

  #structure data based on execution order
  new_data <- new_data[, colnames(new_data) %in% object$Structure$exec_dag$order, drop = FALSE]

  #keep rows with at least one finite entry
  keep     <- apply(new_data, 1, function(r) any(is.finite(as.numeric(r))))
  new_data <- new_data[keep, , drop = FALSE]

  #drop empty columns
  new_data <- new_data[, apply(new_data, 2, function(x) any(is.finite(as.numeric(x)))), drop=FALSE]

  if(nrow(new_data) == 0){stop("after filtering, new_data has 0 rows with finite values")}

  #from the new data create a matrix
  x_obs <- as.matrix(new_data)

  #Make predictions
  set.seed(seed)
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

  #stack all residuals to derive center
  res_stack <- do.call(rbind,
                       lapply(seq_len(nsim), function(i) {
                         x_hat <- do.call(cbind, lapply(preds$Variance[obs_nodes], `[[`, i))
                         obs_mat - x_hat}))

  #derive mad from residuals
  MAD <- apply(res_stack, 2, function(x){mad(x,center = median(x, na.rm = T), constant = mad_constant,na.rm = TRUE)})

  #if scale is basically 0 set some small value
  MAD[MAD==0] <- eps

  ###################################
  #Detection fun (P.S. not real fun)#
  ###################################

  #binary variables
  detect_binary <- function(x){ux <- unique(na.omit(x)); length(ux) <= 2 && all(ux %in% c(0, 1))}

  #count variables
  detect_count <- function(x){x <- na.omit(x); length(x) > 0 && all(x >= 0) && all(abs(x - round(x)) < eps)}

  #proportional data
  detect_beta <- function(x, eps){x <- na.omit(x); if (length(x) == 0){return(FALSE)}; in_range <- all(x >= 0 & x <= 1); many_unique <- length(unique(x)) > 2; in_range && many_unique}

  #transformation function
  transform <- function(y_raw, y_hat, eps, MAD=NULL){

    if(!is.null(MAD)){

      ##############
      #Some log fun#
      ##############

      #binary log loss
      logistic_loss <- function(y_raw, y_hat, eps){y_hat <- pmin(pmax(y_hat, eps), 1-eps); -(y_raw*log(y_hat)+(1-y_raw)*log(1-y_hat))}

      ######################
      #loss and resid funcs#
      ######################

      if(detect_binary(y_raw) || detect_beta(y_raw)){
        p_hat    <- pmin(pmax(y_hat, eps), 1-eps)
        r        <- (y_raw-p_hat)/sqrt(p_hat*(1-p_hat)+eps)
        l        <- logistic_loss(y_raw, y_hat, eps)
      }else if(detect_count(y_raw)){
        r <- sqrt(y_raw+3/8)-sqrt(pmax(y_hat, eps)+3/8)
        l <- abs(r)
      }else{
        if(is.null(MAD)){"why is MAD  null?"}
        r <- (y_raw-y_hat)/MAD
        l <- r^2}

      bad <- !is.finite(l)

      #set to max
      if(any(bad)){
        finite_loss_vals  <- l[is.finite(l)]

        cap_loss          <- if(length(finite_loss_vals)){max(finite_loss_vals)}else{1}
        cap_resid         <- if(length(finite_loss_vals)){r[which.max(abs(r))]}else{1}

        l[bad]            <- cap_loss
        r[bad]            <- cap_resid}

      return(list(loss=l, residuals=r))

    }else{

      if(detect_binary(y_raw) || detect_beta(y_raw)){
        p_hat    <- pmin(pmax(y_hat, eps), 1-eps)
        r        <- (y_raw-p_hat)/sqrt(p_hat*(1-p_hat)+eps)
      }else if(detect_count(y_raw)){
        r <- sqrt(y_raw+3/8)-sqrt(pmax(y_hat, eps)+3/8)
      }else{
        if(is.null(MAD)){"why is MAD  null?"}
        r <- (y_raw-y_hat)/MAD}

      return(r)}}

  #KL optimizer function
  kl_optim  <- function(lambda, loss, kl_target, eps){

    S       <- length(loss)
    a       <- -lambda*loss
    a       <- a-max(a)
    w1      <- exp(a)
    weights <- w1/sum(w1)

    sum(weights*log(pmax(weights, eps)))+log(S)-kl_target}

  #KL stop function
  lambda_kl <- function(loss, kl_target, eps, upper = max_lambda){

    f0    <- kl_optim(0, loss, kl_target, eps)
    f_up  <- kl_optim(upper, loss, kl_target, eps)

    if (!is.finite(f_up) || f_up<0){return(upper)}

    uniroot(kl_optim, interval=c(0, upper), loss=loss, kl_target=kl_target, eps=eps)$root}

  #stable weights to prevent underflow
  stable_weights <- function(loss, lambda, eps){

    a <- -lambda * loss
    a <- a-max(a, na.rm = TRUE)
    w <- exp(a)
    w/sum(w)}

  #Calibrator stuff
  lambda_calibrate <- function(loss, y_obs, y_particles, eps, max_lambda, family){

    ok_obs <- is.finite(y_obs)
    y_obs <- y_obs[ok_obs]
    y_particles <- y_particles[ok_obs, , drop = FALSE]

    ok_part <- is.finite(loss) & apply(y_particles, 2, function(x) all(is.finite(x)))
    loss <- loss[ok_part]
    y_particles <- y_particles[, ok_part, drop = FALSE]

    score <- function(lambda){

      w <- stable_weights(loss, lambda, eps)
      if(any(!is.finite(w)) || sum(w) <= 0) return(.Machine$double.xmax)

      y_pred <- as.numeric(y_particles %*% w)

      if(family %in% c("binary", "beta")){
        p <- pmin(pmax(y_pred, eps), 1 - eps)
        r_data <- (y_obs - p) / sqrt(p * (1 - p) + eps)

        r_part <- apply(y_particles, 2, function(pj){
          pj <- pmin(pmax(pj, eps), 1 - eps)
          mean((y_obs - pj) / sqrt(pj * (1 - pj) + eps), na.rm = TRUE)
        })

      } else if(family == "count"){
        r_data <- sqrt(y_obs + 3/8) - sqrt(pmax(y_pred, eps) + 3/8)

        r_part <- apply(y_particles, 2, function(mu){
          sqrt(y_obs + 3/8) - sqrt(pmax(mu, eps) + 3/8)
        })
        r_part <- colMeans(r_part, na.rm = TRUE)

      } else {
        r_data <- y_obs - y_pred
        r_part <- apply(y_particles, 2, function(mu) mean(y_obs - mu, na.rm = TRUE))
      }

      se_data <- sd(r_data, na.rm = TRUE) / sqrt(length(r_data))
      ESS <- 1 / sum(w^2)

      mu_r <- sum(w * r_part)
      sd_post <- sqrt(sum(w * (r_part - mu_r)^2))
      se_post <- sd_post / sqrt(ESS)

      out <- abs(se_post - se_data)
      if(!is.finite(out)) .Machine$double.xmax else out
    }

    optimize(score, interval = c(0, max_lambda))$minimum
  }

  ############
  #null model#
  ############

  xhat0_array <- array(NA, dim = c(nrow(obs_mat), length(obs_nodes), nsim),
                       dimnames = list(NULL, obs_nodes, NULL))

  for(i in seq_len(nsim)){xhat0_array[,,i] <- do.call(cbind, lapply(preds$Variance[obs_nodes], `[[`, i))}

  y_hat0 <- apply(xhat0_array, c(1, 2), mean, na.rm = TRUE)

  ################
  #loss functions#
  ################

  #create templates for loss and residuals per simulation
  loss0_sim <- loss1_sim <- residual0_sim <- residual1_sim <- matrix(NA_real_, nrow = nsim, ncol = ncol(obs_mat),
                                                                     dimnames = list(NULL, colnames(obs_mat)))

  for (i in seq_len(nsim)){
    if (print_progress) print(i)

    #select simulation i
    x_hat           <- do.call(cbind,lapply(preds$Variance[obs_nodes], `[[`, i))
    x_hat           <- as.matrix(x_hat)
    colnames(x_hat) <- obs_nodes

    #create empty matrix
    loss0 <- loss1  <- residual0 <- residual1  <- matrix(NA_real_, nrow = nrow(obs_mat), ncol = ncol(obs_mat))
    colnames(residual0) <- colnames(residual1) <- colnames(loss0) <- colnames(loss1) <- obs_nodes

    #for column j
    for (j in seq_len(ncol(obs_mat))){

      y_raw  <- obs_mat[,j]
      y_hat  <- x_hat[,j]

      list_transform1 <- transform(y_raw, y_hat, eps, MAD[j])
      list_transform0 <- transform(y_hat, y_hat0[,j], eps, MAD[j])

      loss1_sim[i,j]     <- median(list_transform1$loss)
      residual1_sim[i,j] <- median(list_transform1$residuals)

      loss0_sim[i,j]     <- median(list_transform0$loss)
      residual0_sim[i,j] <- median(list_transform0$residuals)
    }}

  #set inf or NA to max col
  mxr <- function(r){r <- apply(r, 2, function(col) {
    bad <- !is.finite(col)
    if(any(bad)) col[bad] <- max(col[is.finite(col)], na.rm = TRUE)
    col}); r}

  loss0_sim <- mxr(loss0_sim)
  loss1_sim <- mxr(loss1_sim)

  residual0_sim <- mxr(residual0_sim)
  residual1_sim <- mxr(residual1_sim)

  ###########################################
  #Set lambda GBU and derive weights for GBU#
  ###########################################

  node_lambda <- lapply(seq_along(obs_nodes), function(i){

    node  <- obs_nodes[i]

    y_obs <- obs_mat[,i]

    y_particles <- do.call(
      cbind,
      preds$Variance[[node]]
    )

    family <- if(detect_binary(y_obs)) {
      "binary"
    } else if(detect_beta(y_obs, eps)) {
      "beta"
    } else if(detect_count(y_obs)) {
      "count"
    } else {
      "continuous"
    }

    lambda_calibrate(
      loss = loss1_sim[, i],
      y_obs = y_obs,
      y_particles = y_particles,
      eps = eps,
      max_lambda = max_lambda,
      family = family
    )
  })

  #update kl_target
  if(!is.null(kl_frac)){kl_target <- rep(log(1/(1-kl_frac)), ncol(residual1_sim))

  node_lambda <- lapply(1:ncol(loss1_sim), function(i){
    lambda_kl(loss1_sim[,i], kl_target[i], eps)})

  }

  w_node <- lapply(1:ncol(loss1_sim), function(k){
    w <- stable_weights(loss1_sim[,k], node_lambda[[k]], eps)})

  #add names
  names(w_node) <- obs_nodes

  ################
  #global weights#
  ################

  #node_lambda to numeric vector
  node_lambda   <- setNames(as.numeric(unlist(node_lambda)), colnames(loss1_sim))
  lambda_vec    <- node_lambda
  loss1_sim     <- loss1_sim[, names(lambda_vec), drop = FALSE]

  #residual cols in same order
  loss1_sim <- loss1_sim[,names(lambda_vec), drop = FALSE]

  #global log weights
  logw <- -rowSums(sweep(loss1_sim, 2, lambda_vec, "*"), na.rm = TRUE)

  #stabalize prevent overflow
  logw <- logw - max(logw, na.rm = TRUE)

  #prevent issues
  w_global <- exp(logw)

  #global weights
  w_global <- w_global/(sum(w_global)+eps)

  #######################################
  #summarize all the nodes using weights#
  #######################################
  node_summaries     <- list()
  node_params_weight <- list()

  for(node in names(preds$Parameters)){

    draws <- preds$Parameters[[node]]
    if(node %in% object$Roots) next

    w <- w_node[[node]]
    if(is.null(w)){w <- rep(1/nsim, nsim)}

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
    mu_par      <- colSums(sweep(params, 1, w_global, "*"))
    diff_mu_par <- sweep(params, 2, mu_par, "-")
    diff_mu_par[!is.finite(diff_mu_par)] <- 0

    #some jitter to prevent 0 var
    sigma_weighted <- t(diff_mu_par) %*% sweep(diff_mu_par, 1, w_global, "*")
    sigma_weighted <- sigma_weighted + diag(eps, ncol(sigma_weighted))

    rownames(sigma_weighted) <- rownames(object$Sigma)
    colnames(sigma_weighted) <- colnames(object$Sigma)

    warning("estimating the covariance over all parameter can result in degeneracy issues.")}

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
  KL_g                    <- sum(w_global * log((w_global + eps) / (1/nsim)))
  KL_node                 <- lapply(w_node, function(x) sum(x* log((x + eps) / (1/nsim))))
  object$Info             <- list(KL_global = KL_g, rel_KL_div = KL_g/log(nsim), KL_local=KL_node, weights=w_node)

  #extract names
  theta_names <- unlist(Map(function(edge, pars) paste0(edge, "_", names(pars)), names(object$Parameters), object$Parameters))

  #extract se params
  theta_var <- unlist(lapply(object$Parameters, function(e) {sapply(e, function(p) as.numeric(p["se"])^2)}))

  Sigma_new           <- diag(theta_var)
  dimnames(Sigma_new) <- list(theta_names, theta_names)

  if (covar && !is.null(sigma_weighted)) {
    rownames(sigma_weighted) <- theta_names
    colnames(sigma_weighted) <- theta_names

    sigma_weighted[!is.finite(sigma_weighted)] <- 0
    diag(sigma_weighted) <- pmax(diag(sigma_weighted), eps)

    R <- cov2cor(sigma_weighted)
    R[!is.finite(R)] <- 0
    diag(R) <- 1

    sd_new <- sqrt(pmax(theta_var, eps))
    object$Sigma <- diag(sd_new) %*% R %*% diag(sd_new)
    object$Sigma[!is.finite(object$Sigma)] <- 0
    diag(object$Sigma) <- pmax(diag(object$Sigma), eps)
    dimnames(object$Sigma) <- list(theta_names, theta_names)

    #if covar is flase
  }else{
    object$Sigma <- Sigma_new}

  #############################
  #Update residual diagnostics#
  #############################
  #Make predictions with updated mod
  preds_update <- predict_ppmn(object, new_data, nsim = 1)

  #select predicted nodes updated mod
  pred_nodes_update <- setdiff(names(preds_update$Expected), object$Roots)

  x_hat           <- do.call(cbind, preds_update$Expected[obs_nodes])
  x_hat           <- as.matrix(x_hat)
  colnames(x_hat) <- obs_nodes

  mat           <- matrix(NA_real_, nrow = nrow(obs_mat), ncol = ncol(obs_mat))
  colnames(mat) <- obs_nodes

  for (j in seq_len(ncol(obs_mat))){

    y_raw <- obs_mat[, j]
    y_hat <- x_hat[, j]

    if(detect_binary(y_raw)){
      p_hat    <- pmin(pmax(y_hat, eps), 1-eps)
      mat[, j] <- (y_raw-p_hat)/sqrt(p_hat*(1-p_hat)+eps)

    }else if(detect_beta(y_raw, eps)){
      p_hat    <- pmin(pmax(y_hat, eps), 1-eps)
      mat[, j] <- (y_raw-p_hat)/sqrt(p_hat*(1-p_hat)+eps)

    }else if(detect_count(y_raw)){
      mat[, j] <- sqrt(y_raw+3/8)-sqrt(pmax(y_hat, eps)+3/8)

    }else{
      mat[, j] <- (y_raw-y_hat)
    }
  }

  #calculate mad
  MAD <- apply(mat, 2, function(x)
    mad(x, center = median(x, na.rm = TRUE),
        constant = mad_constant, na.rm = TRUE))

  #set mad = 0 to small nr
  MAD[MAD == 0] <- eps

  for(j in seq_len(ncol(mat))) {

    y_raw <- obs_mat[, j]

    if (!detect_binary(y_raw) && !detect_count(y_raw)) {
      mat[, j] <- mat[, j] / MAD[j]
    }
  }

  #kl-divergence function
  kl_div_fun <- function(w){kl <- sum(w * log((w + eps) / (1/nsim))); ifelse(kl<0,0,kl)}

  kl_node   <- apply(do.call(cbind, w_node), 2, kl_div_fun)
  res_medsd <- apply(mat, 2, function(x) c("median"=median(x, na.rm = T),
                                           "mu"=mean(x, na.rm = T)))
  #create labels
  forlabel <- round(rbind(res_medsd, kl_node),2)
  labels   <- apply(forlabel, 2, function(x) paste0("\nmedian=", x[1], "\nmu=", x[2], "\nKL=", x[3]))
  labels   <- paste0(obs_nodes, labels)

  #create long df
  list_long         <- lapply(seq_len(ncol(mat)), function(x) data.frame(labels[x], na.omit(mat[,x])))
  df_long           <- do.call(rbind, list_long)
  colnames(df_long) <- c("Variable", "Value")

  #quantile limits of plot
  qy <- quantile(df_long$Value, c(.025,.975), na.rm = T)

  #residual plot
  rpl <- ggplot(df_long, aes(x=Variable, y=Value))+ylim(qy)+
    geom_hline(yintercept = 0, col="tomato3", lty=2)+ylab("Standardized residuals")+
    geom_boxplot(outlier.shape = NA)+geom_jitter(alpha=0.2, pch=19, width = 0.2)+
    theme_classic()

  object$Residuals$list          <- df_long
  object$Residuals$residual_plot <- rpl

  return(object)
}
