ppmn_fit_metric <- function(object, data, level=0.9, round=2){

  #check roots
  miss_roots <- setdiff(object$Roots, colnames(data))
  if(length(miss_roots)>0){stop("missing root nodes: ", paste(miss_roots, collapse = ", "))}

  #hdi function
  hdi_fun              <- function(x, level){

    orddata   <- sort(x)
    nord      <- length(x)
    infomass  <- ceiling(level*nord)
    outmass   <- nord-infomass

    min_width <- Inf
    ll        <- NA
    ul        <- NA

    for(i in 1:(outmass+1)){
      int_width <- orddata[i+infomass-1]-orddata[i]

      if(int_width < min_width){
        min_width <- int_width
        ll <- orddata[i]
        ul <- orddata[i+infomass-1]}}

    c(ll, ul)}

  #select predicted variance
  pred_nodes <- setdiff(names(object$Variance), object$Roots)
  obs_nodes  <- intersect(colnames(data), pred_nodes)

  if(length(obs_nodes)==0){return(NA)}

  obs_mat <- as.matrix(data[, obs_nodes, drop = FALSE])
  obs_mat <- apply(obs_mat, 2, as.numeric)

  mae     <-  matrix(NA, nrow=object$nsim, ncol=length(obs_nodes))
  rsq     <-  matrix(NA, nrow=object$nsim, ncol=length(obs_nodes))

  for(s in seq_len(object$nsim)){
    pred_mat <- do.call(cbind, lapply(object$Variance[obs_nodes], `[[`, s))

    pred_mat           <- as.matrix(pred_mat)
    colnames(pred_mat) <- obs_nodes

    ae                 <- abs(obs_mat-pred_mat)
    ae[!is.finite(ae)] <- NA


    if(all(is.na(ae))){next}

    mae[s,] <- apply(ae, 2, function(x) median(x, na.rm = T))
    rsq[s,] <- unlist(lapply(seq_len(ncol(obs_mat)), function(x) cor(obs_mat[,x], pred_mat[,x], use = "pairwise.complete.obs")))^2}

  mae_df <- round(apply(mae, 2, function(x) c(median(x, na.rm = T), hdi_fun(x, level))),round)
  rsq_df <- round(apply(rsq, 2, function(x) c(median(x, na.rm = T), hdi_fun(x, level))),round)
  colnames(mae_df) <- obs_nodes
  colnames(rsq_df) <- obs_nodes
  rownames(rsq_df) <- rownames(mae_df) <- c("median", "ll", "ul")

  list(MAE=mae_df, Rsq=rsq_df)}
