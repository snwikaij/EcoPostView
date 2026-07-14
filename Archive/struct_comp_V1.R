struct_comp <- function(object, new_data, nsim = 3000, level=0.9,
                        kl_frac=NULL, mad_constant=1.4826, seed=123,
                        max_lambda=1000, print_progress=F, rotate_x=0, hjust=0){

  if(nsim<1){stop("number of simulations cannot be smaller than 1")}
  if(level>0.99999){stop("level cannot be larger than .99999")}
  if(level<0.00001){stop("level cannot be smaller than .00001")}
  if(!is.null(kl_frac)){
    if(kl_frac>0.99){stop("update fraction cannot be larger than .99")}
    if(kl_frac<0.01){stop("update fraction cannot be smaller than .01")}}
  if(mad_constant<0.01){stop("mad_constant cannot be smaller than 0.01")}

  #update kl_target
  kl_target <- log(1/(1-kl_frac))

  #constant
  eps <- .Machine$double.eps

  #create the null model
  make_M0 <- function(object, new_data, eps = .Machine$double.eps){

    obj0  <- object
    funcs <- obj0$Structure$predict_dag$functions
    edges <- names(obj0$Parameters)

    for (edge in edges) {

      dep  <- sub("~.*", "", edge)
      func <- funcs[[edge]]
      resp <- new_data[[dep]]

      if(is.null(resp)) next

      resp <- resp[is.finite(resp)]

      if(length(resp) == 0) next

      y_ident <- mean(resp, na.rm = TRUE)

      y_pos <- resp[resp > 0 & is.finite(resp)]
      y_log <- if (length(y_pos)) mean(log(y_pos), na.rm = TRUE) else log(eps)

      y_prop <- resp[resp > 0 & resp < 1 & is.finite(resp)]
      y_logit <- if (length(y_prop)){
        mean(qlogis(pmin(pmax(y_prop, eps), 1 - eps)), na.rm = TRUE)
      }else {0}

      bnames <- names(obj0$Parameters[[edge]])

      for (b in bnames) {
        obj0$Parameters[[edge]][[b]][["mu"]] <- 0}

      if(func == "identity") {
        if("b0" %in% bnames){obj0$Parameters[[edge]][["b0"]][["mu"]] <- y_ident}
      }else if(func == "log"){
        if("b0" %in% bnames){obj0$Parameters[[edge]][["b0"]][["mu"]] <- y_log}
      }else if (func == "logit") {
        if("b0" %in% bnames){obj0$Parameters[[edge]][["b0"]][["mu"]] <- y_logit}
      }else if (func == "sigmoidal"){
        if("b0" %in% bnames){obj0$Parameters[[edge]][["b0"]][["mu"]] <- y_ident}
        if("b1" %in% bnames){obj0$Parameters[[edge]][["b1"]][["mu"]] <- 0}
        if("b2" %in% bnames){obj0$Parameters[[edge]][["b2"]][["mu"]] <- 1e6}

      }else if (func == "gompertz"){
        if ("b0" %in% bnames){obj0$Parameters[[edge]][["b0"]][["mu"]] <- object$Parameters[[edge]][["b0"]][["mu"]]}
        if ("b1" %in% bnames){obj0$Parameters[[edge]][["b1"]][["mu"]] <- 0}
        if ("b2" %in% bnames){obj0$Parameters[[edge]][["b2"]][["mu"]] <- 0}

      }else if (func == "gaussian"){
        if ("b0" %in% bnames){obj0$Parameters[[edge]][["b0"]][["mu"]] <- object$Parameters[[edge]][["b0"]][["mu"]]}
        if ("b1" %in% bnames){obj0$Parameters[[edge]][["b1"]][["mu"]] <- 0}
        if ("b2" %in% bnames){obj0$Parameters[[edge]][["b2"]][["mu"]] <- 1e6}

      }else if (func == "class"){
        mu0 <- (object$Parameters[[edge]][["b0"]][["mu"]]+object$Parameters[[edge]][["b2"]][["mu"]])/2
        sd0 <- sqrt((object$Parameters[[edge]][["b1"]][["mu"]]^2+object$Parameters[[edge]][["b3"]][["mu"]]^2)/2)

        if("b0" %in% bnames){obj0$Parameters[[edge]][["b0"]][["mu"]] <- mu0}
        if("b2" %in% bnames){obj0$Parameters[[edge]][["b2"]][["mu"]] <- mu0}

        if("b1" %in% bnames){obj0$Parameters[[edge]][["b1"]][["mu"]] <- sd0}
        if("b3" %in% bnames){obj0$Parameters[[edge]][["b3"]][["mu"]] <- sd0}

        if("b4" %in% bnames){obj0$Parameters[[edge]][["b4"]][["mu"]] <- 0.5}
      }
    }

    obj0
  }

  #structure data based on execution order
  new_data <- new_data[, colnames(new_data) %in% object$Structure$exec_dag$order, drop = FALSE]

  #keep rows with at least one finite entry
  keep     <- apply(new_data, 1, function(r) any(is.finite(as.numeric(r))))
  new_data <- new_data[keep, , drop = FALSE]

  #drop empty columns
  new_data <- new_data[, apply(new_data, 2, function(x) any(is.finite(as.numeric(x)))), drop=FALSE]

  if (nrow(new_data) == 0){stop("after filtering, new_data has 0 rows with finite values")}

  #from the new data create a matrix
  x_obs <- as.matrix(new_data)

  #Create M1 and M0
  objectM1 <- object
  objectM0 <- make_M0(object, new_data)

  #Make predictions
  set.seed(seed)
  preds1 <- predict_ppmn(objectM1, new_data, nsim = nsim)
  preds0 <- predict_ppmn(objectM0, new_data, nsim = nsim)

  #select predicted nodes
  pred_nodes_all <- setdiff(names(preds1$Variance), preds1$Roots)

  #select observed nodes/variables
  obs_nodes <- intersect(colnames(x_obs), pred_nodes_all)

  if (length(obs_nodes) == 0) {stop("No observed (non-root) nodes in new_data match model predicted nodes")}

  ###########################
  #standardize over all nsim#
  ###########################
  #pre-extract observed vectors for speed
  obs_mat <- as.matrix(x_obs[, obs_nodes, drop = FALSE])

  #stack all residuals to derive center
  stack_resids <- function(predictions, observations, nsim){
    do.call(rbind,
            lapply(seq_len(nsim), function(i) {
              x_hat <- do.call(cbind, lapply(predictions$Variance[obs_nodes], `[[`, i))
              obs_mat - x_hat}))}

  stacked_M1 <- stack_resids(preds1, obs_mat, nsim)
  stacked_M0 <- stack_resids(preds0, obs_mat, nsim)
  sres       <- rbind(stacked_M1, stacked_M0)

  #derive mad from residuals
  MAD <- apply(sres, 2, function(x) {mad(x,center = median(x, na.rm = T), constant = mad_constant,na.rm = TRUE)})
  MAD[MAD==0] <- eps

  #log cosh fun
  log_cosh <- function(z){log(cosh(z))}

  #sqrt loss func
  sqrt_loss <- function(y, pred, eps){abs(sqrt(y) - sqrt(pmax(pred, eps)))}

  #binary log loss
  logistic_loss <- function(y, p, eps){y <- pmin(pmax(y, eps), 1 - eps); p <- pmin(pmax(p, eps), 1 - eps); - (y * log(p) + (1 - y) * log(1 - p))}

  #detect binary variables
  detect_binary <- function(x){ux <- unique(na.omit(x)); length(ux) <= 2 && all(ux %in% c(0, 1))}

  #detect count variables
  detect_count <- function(x){x <- na.omit(x); length(x) > 0 && all(x >= 0) && all(abs(x - round(x)) < eps)}

  #detect proportional data
  detect_beta <- function(x, eps){x <- na.omit(x); if (length(x) == 0){return(FALSE)}; in_range <- all(x >= 0 & x <= 1); many_unique <- length(unique(x)) > 2; in_range && many_unique}

  #KL optimizer function
  kl_optim <- function(lambda, loss, kl_target, eps){

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

  lambda_calibrate <- function(loss, y_obs, y_particles, eps, max_lambda, family = "continues"){

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

  ###########
  #MERR dist#
  ###########

  loss_fun <- function(predictions, obs_mat, obs_nodes, eps, MAD){

    residuals           <- matrix(NA_real_, nrow = nsim, ncol = ncol(obs_mat))
    colnames(residuals) <- colnames(obs_mat)

    for (i in seq_len(nsim)) {
      if (print_progress) print(i)

      #select simulation i
      x_hat           <- do.call(cbind,lapply(predictions$Variance[obs_nodes], `[[`, i))
      x_hat           <- as.matrix(x_hat)
      colnames(x_hat) <- obs_nodes

      #create empty matrix
      loss           <- matrix(NA_real_, nrow = nrow(obs_mat), ncol = ncol(obs_mat))
      colnames(loss) <- obs_nodes

      #for column j
      for (j in seq_len(ncol(obs_mat))){

        y_raw  <- obs_mat[,j]
        y_hat  <- x_hat[,j]

        #detect bin or prop
        if (detect_binary(y_raw) || detect_beta(y_raw)){

          loss[, j] <- logistic_loss(y_raw, y_hat, eps)
          bad       <- !is.finite(loss[,j])

          if(any(bad)) {
            finite_vals <- loss[is.finite(loss[, j]), j]
            cap         <- if (length(finite_vals)){max(finite_vals)}else{1}
            loss[bad,j] <- cap
          }

          #detect count values
        }else if(detect_count(y_raw)){
          loss[, j] <- sqrt_loss(y_raw, y_hat, eps)
          bad <- !is.finite(loss[,j])

          #set max otherwise 1
          if(any(bad)) {
            finite_vals <- loss[is.finite(loss[, j]), j]
            cap          <- if(length(finite_vals)){max(finite_vals)}else{1}
            loss[bad,j]  <- cap}

          #detect continues values
        }else{
          loss[, j] <- log_cosh((y_hat-y_raw)/MAD[j])
          bad <- !is.finite(loss[,j])

          #set to max
          if(any(bad)) {
            finite_vals <- loss[is.finite(loss[, j]), j]
            cap         <- if(length(finite_vals)){max(finite_vals)}else{1}
            loss[bad,j] <- cap}
        }
      }

      residuals[i, ] <- colMeans(loss, na.rm = TRUE)}

    residuals <- apply(residuals, 2, function(col) {
      bad <- !is.finite(col)
      if (any(bad)) col[bad] <- max(col[is.finite(col)], na.rm = TRUE)
      col})

    colnames(residuals) <- obs_nodes
    residuals}

  residM1 <- loss_fun(preds1, obs_mat, obs_nodes, eps, MAD)
  residM0 <- loss_fun(preds0, obs_mat, obs_nodes, eps, MAD)

  ###############################
  #derive weights for GBU Lambda#
  ###############################

  node_lambda <- lapply(seq_along(obs_nodes), function(i){

    node  <- obs_nodes[i]

    y_obs <- obs_mat[,i]

    y_particles <- do.call(
      cbind,
      preds1$Variance[[node]]
    )

    lambda_calibrate(
      loss = residM1[, i],
      y_obs = y_obs,
      y_particles = y_particles,
      eps = eps,
      max_lambda = max_lambda
    )
  })

  #update kl_target
  if(!is.null(kl_frac)){kl_target <- rep(log(1/(1-kl_frac)), ncol(residM1))

  node_lambda <- lapply(1:ncol(residM1), function(i){
    lambda_kl(residM1[,i], kl_target[i], eps)})

  }

  #weight of m1
  w_node1 <- lapply(1:ncol(residM1), function(k){
    w <- stable_weights(residM1[,k], node_lambda[[k]], eps)})

  #weight of m0
  w_node0 <- lapply(1:ncol(residM0), function(k){
    w <- stable_weights(residM0[,k], node_lambda[[k]], eps)})

  #add names
  names(w_node1) <- names(w_node0) <- obs_nodes

  #KL divergence
  kl_div <- unlist(lapply(1:ncol(residM1), function(i){
    p  <- w_node1[[i]]+eps; p  <- p/sum(p)
    q  <- w_node0[[i]]+eps; q  <- q/sum(q)
    sum(p*log(p/q))}))

  #add names again
  names(kl_div) <- obs_nodes

  #place into a dataset
  df              <- data.frame(Variables=obs_nodes,
                                Value=kl_div)
  rownames(df)    <- NULL

  colvals <- 1-level

  ggplot(df, aes(x=Variables, y=Value))+
    geom_bar(stat = "identity", col="black", fill="grey70")+
    theme_classic()+ylab("KL-Divergence (nats)")+
    theme(axis.text.x = element_text(angle=rotate_x, hjust = hjust),
          axis.title.x = element_blank())}
