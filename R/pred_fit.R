#' Predictive fit function to assess performance
#'
#' @param object  A foundational or updated PPMN
#' @param new_data New dataset with all root nodes and least one child node present
#' @param nsim Number of simulations
#' @param level Credibility level (default level = 0.9)
#' @param rotate_x The angle to rotate the variable names on x-axis
#' @param hjust A number determining the height of the variable names on the x-axis
#' @param seed Seed value 123
#'
#' @export
pred_fit <- function(object, new_data, nsim = 1000, level=0.9,
                     rotate_x=0, hjust=0, seed=123){

  if(nsim<1){stop("number of simulations cannot be smaller than 1")}
  if(level>0.99999){stop("level cannot be larger than .99999")}
  if(level<0.00001){stop("level cannot be smaller than .00001")}

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

  #Make predictions
  set.seed(seed)
  preds1 <- predict_ppmn(object, new_data, nsim = nsim)

  #select predicted nodes
  pred_nodes_all <- setdiff(names(preds1$Variance), preds1$Roots)

  #select observed nodes/variables
  obs_nodes <- intersect(colnames(x_obs), pred_nodes_all)

  if (length(obs_nodes) == 0) {stop("No observed (non-root) nodes in new_data match model predicted nodes")}

  #pre-extract observed vectors for speed
  obs_mat <- as.matrix(x_obs[, obs_nodes, drop = FALSE])

  maxpost              <- function(x){d <- density(x); d$x[which.max(d$y)]}

  cors <- lapply(seq_len(nsim), function(i){
    x_hat <- do.call(cbind, lapply(preds1$Variance[obs_nodes], `[[`, i))
    unlist(lapply(seq_along(obs_nodes), function(j){
      ok <- is.finite(obs_mat[, j]) & is.finite(x_hat[, j])
      if(sum(ok) < 3 || sd(obs_mat[ok, j]) == 0 || sd(x_hat[ok, j]) == 0) return(NA_real_)
      cor(x_hat[ok, j], obs_mat[ok, j])}))})
  cors <- do.call(rbind, cors)

  cors_df <- round(data.frame(row.names = obs_nodes,
                              median = apply(cors, 2, median, na.rm = TRUE),
                              map = apply(cors, 2, maxpost),
                              mu = apply(cors, 2, mean, na.rm = TRUE),
                              ll = apply(cors, 2, quantile, probs = (1 - level)/2, na.rm = TRUE),
                              ul = apply(cors, 2, quantile, probs = 1 - (1 - level)/2, na.rm = TRUE)),2)

  mae <- lapply(seq_len(nsim), function(i){
    x_hat <- do.call(cbind, lapply(preds1$Variance[obs_nodes], `[[`, i))
    unlist(lapply(seq_len(length(obs_nodes)), function(j){mean(x_hat[,j]-obs_mat[,j], na.rm=T)}))})
  mae <- do.call(rbind, mae)
  mae <- sweep(mae, 2,  apply(mae, 2, sd), "/")

  mae_df <- round(data.frame(row.names = obs_nodes,
                             median = apply(mae, 2, median, na.rm = TRUE),
                             map = apply(mae, 2, maxpost),
                             mu = apply(mae, 2, mean, na.rm = TRUE),
                             ll = apply(mae, 2, quantile, probs = (1 - level)/2, na.rm = TRUE),
                             ul = apply(mae, 2, quantile, probs = 1 - (1 - level)/2, na.rm = TRUE)),2)

  pbal <- lapply(seq_len(nsim), function(i){
    x_hat <- do.call(cbind, lapply(preds1$Variance[obs_nodes], `[[`, i))
    unlist(lapply(seq_along(obs_nodes), function(j){
      mean(obs_mat[, j] > x_hat[, j], na.rm = TRUE)
    }))})
  pbal <- do.call(rbind, pbal)

  pbal_df <- round(data.frame(row.names = obs_nodes,
                              median = apply(pbal, 2, median, na.rm = TRUE),
                              map = apply(pbal, 2, maxpost),
                              mu = apply(pbal, 2, mean, na.rm = TRUE),
                              ll = apply(pbal, 2, quantile, probs = (1 - level)/2, na.rm = TRUE),
                              ul = apply(pbal, 2, quantile, probs = 1 - (1 - level)/2, na.rm = TRUE)),2)

  mae_plot <- ggplot(mae_df, aes(x=rownames(mae_df), y=median))+
    geom_point(pch=19, size=3)+xlab("Variables")+
    geom_errorbar(aes(ymin = ll, ymax = ul),
                  width = 0, linewidth = 0.5)+
    theme_classic()+ylab("Standardized Residual Error (z)")+
    geom_hline(yintercept = 0, lty=2, lwd=0.6, col="tomato3")

  bal_plot <- ggplot(pbal_df, aes(x=rownames(pbal_df), y=median))+
    geom_point(pch=19, size=3)+
    geom_errorbar(aes(ymin = ll, ymax = ul),
                  width = 0, linewidth = 0.5)+
    geom_hline(yintercept = 0.5, lwd=0.6, col="tomato3", lty=2)+
    theme_classic()+ylab("Balance \nP(observed>predicted)")+
    theme(axis.text.x = element_text(angle=rotate_x, hjust = hjust),
          axis.title.x = element_blank())

  r_plot <- ggplot(cors_df, aes(x=rownames(cors_df), y=median))+
    geom_point(pch=19, size=3)+
    geom_errorbar(aes(ymin = ll, ymax = ul),
                  width = 0, linewidth = 0.5)+
    geom_hline(yintercept = 0, lwd=0.6, col="tomato3", lty=2)+
    theme_classic()+ylab("Correlation \n R(observed, predicted)")+
    theme(axis.text.x = element_text(angle=rotate_x, hjust = hjust),
          axis.title.x = element_blank())

  list(Stdz=mae_df, Balance=pbal_df, correlation=cors_df, plots=list(mae=mae_plot, balance=bal_plot, cor=r_plot))}
