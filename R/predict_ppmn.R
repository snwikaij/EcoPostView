#' Predict function for a PPMN object
#'
#' @param object A foundational or updated PPMN
#' @param new_data New data for prediction
#' @param nsim Number of simulations
#' @param print_drop Display the simulation progress (default print_drop = FALSE)
#'
#' @export
predict_ppmn <- function(object, new_data, nsim = 1000,
                         print_drop = FALSE) {

  input_values <- as.list(new_data)
  input_values <- input_values[object$Roots]

  node_names  <- igraph::V(object$Structure$visual_dag$graph)$name
  given_nodes <- intersect(names(new_data), node_names)

  #mu lower and upper truncation
  mu_global    <- unlist(lapply(object$Parameters, function(edge) {sapply(edge, function(p) p["mu"])}))
  lower_global <- unlist(lapply(object$Parameters, function(edge) {sapply(edge, function(p) p["a"])}))
  upper_global <- unlist(lapply(object$Parameters, function(edge) {sapply(edge, function(p) p["b"])}))

  #set names to extracted mu a and b
  names(mu_global)    <- rownames(object$Sigma)
  names(lower_global) <- rownames(object$Sigma)
  names(upper_global) <- rownames(object$Sigma)

  ################
  #MC simulations#
  ################
  if (all(is.infinite(lower_global)) && all(is.infinite(upper_global))) {
    theta_draws <- MASS::mvrnorm(n=nsim, mu=mu_global,Sigma=object$Sigma)
  }else{
    theta_draws <- tmvtnorm::rtmvnorm(n=nsim, mean=mu_global, sigma=object$Sigma,lower=lower_global, upper=upper_global)}

  #store the draws in matrix and give names
  theta_draws <- matrix(theta_draws,
                        nrow = nsim,
                        dimnames = list(seq_len(nsim), names(mu_global)))

  #select the dag structure
  g          <- object$Structure$visual_dag$graph
  edge_table <- as.data.frame(object$Structure$visual_dag$edges)
  funcs      <- object$Structure$predict_dag$functions
  trans      <- object$Structure$predict_dag$transformations
  topo_order <- igraph::topo_sort(g, mode = "out")

  #store predictions and parameters of mc simulations
  predictions <- vector("list", nsim)
  parameters  <- vector("list", nsim)

  for (s in seq_len(nsim)) {
    predictions[[s]] <- input_values
    parameters[[s]]  <- list()}

  #edge equation formula selection
  get_formula_for_node <- function(node) {
    parents <- edge_table$indep[edge_table$dep == node]
    sig     <- paste(sort(parents), collapse = "+")
    f       <- names(funcs)[sapply(names(funcs), function(x){
      p <- strsplit(x, "~")[[1]]
      p[1] == node && paste(sort(strsplit(p[2], "\\+")[[1]]), collapse = "+") == sig})]
    f[1]}

  ##########
  #DAG loop#
  ##########
  Expected <- list()
  Variance <- list()
  for (i in topo_order) {

    node         <- igraph::V(g)$name[i]
    parents      <- edge_table$indep[edge_table$dep == node]
    if (length(parents) == 0){next}

    edge_name <- get_formula_for_node(node)
    edge_pars <- object$Parameters[[edge_name]]

    param_idx <- paste0(edge_name, "_", names(edge_pars))

    for (s in seq_len(nsim)) {

      parent_vals <- lapply(parents, function(p) predictions[[s]][[p]])
      names(parent_vals) <- parents

      x    <- do.call(cbind,lapply(parents, function(p) predictions[[s]][[p]]))
      beta <- as.numeric(theta_draws[s, param_idx, drop = FALSE])

      parameters[[s]][[node]]  <- beta
      yhat                     <- predict_edge(funcs[[edge_name]], beta, trans[[edge_name]], x)
      yhat[!is.finite(yhat)]   <- NA_real_
      predictions[[s]][[node]] <- yhat}

    beta_mu          <- mu_global[param_idx]
    Expected[[node]] <- as.numeric(
      predict_edge(funcs[[edge_name]], beta_mu, trans[[edge_name]], x))
  }

  Expected <- c(input_values, Expected)

  #summarize variance and parameters
  Params_used <- list()
  for (node in names(predictions[[1]])) {
    Variance[[node]]    <- lapply(predictions, `[[`, node)
    Params_used[[node]] <- lapply(parameters, `[[`, node)}

  if(is.null(dim(new_data))){Data <- new_data[given_nodes]
  }else{Data <-new_data[, given_nodes, drop = FALSE]}

  #output
  list(Expected = Expected, Variance = Variance, Data = Data,
       Structure = object$Structure, Parameters = Params_used, nsim = nsim, Roots = object$Roots)}
