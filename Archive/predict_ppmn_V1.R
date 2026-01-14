
predict_ppmn <- function(object, new_data, nsim=1000, interval_type="CI", print_drop=F){

  #Expected value and variance of y E(y) Var(E(y))
  expected_y     <- function(func, params, trans, x){

    betas <- sapply(params, function(p) p["mu"])

    if(is.na(trans)) {
      x_t <- x
    }else if(trans == "log"){
      x_t <- log(x)
    }else if (trans == "sqrt"){
      x_t <- sqrt(x)
    }else{
      return(NA)}

    if (is.na(func)){func <- "identity"}

    if(func %in%  c("identity", "log", "logit")){
      if(length(betas) == 1){stop("At least two parameters needed for a linear model: b0 and b1.")}
      if (is.numeric(x_t) && is.null(dim(x_t))) {
        eta <- betas[1] + betas[-1] * x_t
      } else if (is.matrix(x_t)) {
        eta <- betas[1] + drop(x_t %*% betas[-1])
      } else {
        stop("x must be a numeric vector or matrix.")
      }}

    if(func == "identity") {
      eta
    }else if (func == "log") {
      exp(eta)
    }else if (func == "logit") {
      plogis(eta)
    }else if (func == "sigmoidal") {
      if (length(betas) != 3){stop("Three parameters needed for a sigmoidal model: b0, b1 and b2.")}
      as.numeric(betas[1] / (1 + exp((betas[2] - x_t) / betas[3])))
    }else if (func == "gaussian") {
      if (length(betas) != 3){stop("Three parameters needed for a gaussian model: b0, b1 and b2.")}
      as.numeric(betas[1] * exp(-0.5 * ((x - betas[2])/betas[3])^2))
    } else if (func == "class") {
      if (length(betas) != 5){stop("Five parameters needed for the Bayesian classifier: b0, b1, b2, b3 and b4.")}
      log_num <- dnorm(x_t, betas[1], betas[2], log = TRUE)
      log_den <- dnorm(x_t, betas[3], betas[4], log = TRUE)
      log_likelihood <- log_num - log_den
      as.numeric((exp(log_likelihood) * betas[5]) / ((exp(log_likelihood) * betas[5]) + (1 - betas[5])))
    } else {
      NA
    }
  }
  variance_y     <- function(func, params, trans, input_value, interval_type){

    x     <- input_value

    if(interval_type == "PI"){
      betas <- sapply(params, function(p) truncnorm::rtruncnorm(1, mean = p["mu"], sd = p["sd"], a = p["a"], b = p["b"]))
    }else if(interval_type == "CI"){
      betas <- sapply(params, function(p) truncnorm::rtruncnorm(1, mean = p["mu"], sd = p["se"], a = p["a"], b = p["b"]))
    }

    if(is.na(trans)) {
      x_t <- x
    }else if (trans == "log") {
      x_t <- log(x)
    }else if (trans == "sqrt") {
      x_t <- sqrt(x)
    }else {
      return(NA) }

    if (is.na(func)){func <- "identity"}

    if (func %in% c("identity", "log", "logit")) {
      if (length(betas) == 1){stop("At least two parameters needed for a linear model: b0 and b1.")}

      if (is.numeric(x_t) && is.null(dim(x_t))) {
        eta <- betas[1] + betas[-1] * x_t
      } else if (is.matrix(x_t)) {
        eta <- betas[1] + drop(x_t %*% betas[-1])
      } else {
        stop("x must be a numeric vector or matrix.")
      }

    }

    if (func == "identity") {
      list(eta, betas)
    } else if (func == "log") {
      list(exp(eta), betas)
    } else if (func == "logit") {
      list(plogis(eta), betas)
    } else if (func == "sigmoidal") {
      if (length(betas) != 3){stop("Three parameters needed for a sigmoidal model: b0, b1 and b2.")}
      eta <- as.numeric(betas[1] / (1 + exp((betas[2] - x_t) / betas[3])))
      list(eta, betas)
    }else if (func == "gaussian") {
      if (length(betas) != 3){stop("Three parameters needed for a gaussian model: b0, b1 and b2.")}
      eta <- as.numeric(betas[1] * exp(-0.5 * ((x - betas[2])/betas[3])^2))
      list(eta, betas)
    } else if (func == "class") {
      if (length(betas) != 5){stop("Four parameters needed for the Bayesian classifier: b0, b1, b2, b3 and b4.")}
      log_num <- dnorm(x_t, betas[1], betas[2], log = TRUE)
      log_den <- dnorm(x_t, betas[3], betas[4], log = TRUE)
      log_likelihood <- log_num - log_den
      eta <- as.numeric((exp(log_likelihood) * betas[5]) / ((exp(log_likelihood) * betas[5]) + (1 - betas[5])))
      list(eta, betas)
    } else {
      NA
    }
  }

  #transform to list
  input_values   <- as.list(new_data)

  #Check input_values only against true roots
  missing_inputs <- setdiff(object$Roots, names(input_values))
  if(length(missing_inputs) > 0){
    stop(paste0("The variable(s): ", paste(missing_inputs, collapse=", ")," is/are not given in the input."))}

  #Select nodes drop not used
  node_names  <- V(object$Structure$visual_dag$graph)$name
  given_nodes <- node_names[node_names %in% names(input_values)]
  not_given   <- names(input_values)[!names(input_values) %in% node_names]
  if(length(not_given)==0){not_given <- "none"}

  input_values   <- input_values[object$Roots]

  #Print dropped input
  if(print_drop==T){
    print(paste("the nodes", paste(given_nodes, collapse = ", "), "are given and",
                paste(not_given, collapse = ", "), "are dropped"))}

  #Function to compute predictions in DAG order
  compute_dag_predictions <- function(structure, parameters, input_values, nsim) {

    funcs      <- structure$predict_dag$functions
    trans      <- structure$predict_dag$transformations
    edge_table <- as.data.frame(structure$visual_dag$edges)
    g          <- structure$visual_dag$graph
    topo_order <- topo_sort(g, mode = "out")

    predictions_expected <- list()
    predictions_variance <- list()
    parameter_list       <- list()

    root_nodes <- names(input_values)
    for (root in root_nodes) {
      predictions_expected[[root]] <- input_values[[root]]
      predictions_variance[[root]] <- input_values[[root]]
    }

    get_formula_for_node <- function(node, funcs, edge_table) {
      parents <- edge_table$indep[edge_table$dep == node]

      # no parents → root node
      if (length(parents) == 0) return(NULL)

      parent_signature <- paste(sort(parents), collapse = "+")

      match_idx <- sapply(names(funcs), function(f) {
        parts <- strsplit(f, "~")[[1]]
        lhs   <- parts[1]
        rhs   <- paste(sort(strsplit(parts[2], "\\+")[[1]]), collapse = "+")
        lhs == node && rhs == parent_signature
      })

      if (!any(match_idx)) {
        stop(paste0(
          "No matching formula for node ", node,
          " with parents: ", paste(parents, collapse = ", ")
        ))
      }

      names(funcs)[which(match_idx)]
    }

    #create a list with predictions if there is not list
    normalize_parent <- function(x, nsim = nsim) {
      if (is.list(x)) {
        return(x)
      } else if (is.numeric(x) && length(x) == length(input_values[[1]])) {
        return(replicate(nsim, x, simplify = FALSE))
      } else {
        stop("Unknown parent type")
      }
    }

    # MAIN LOOP THROUGH DAG IN TOPOLOGICAL ORDER
    for (i in topo_order) {

      node <- igraph::V(g)$name[i]

      #skip root nodes
      parents <- edge_table$indep[edge_table$dep == node]
      if (length(parents) == 0) {
        next
      }

      #fixed stuff
      formula_name   <- get_formula_for_node(node, funcs, edge_table)

      edge_func      <- funcs[[formula_name]]
      edge_trans     <- trans[[formula_name]]
      edge_params    <- parameters[[formula_name]]

      #crete input input for the expected values
      input_expected               <- do.call(cbind, predictions_expected[parents])
      predictions_expected[[node]] <- expected_y(edge_func, edge_params, edge_trans, input_expected)

      #create input for variance
      parent_list    <- lapply(predictions_variance[parents], normalize_parent, nsim = nsim)
      df_list        <- lapply(seq_len(nsim), function(i) { do.call(cbind, lapply(parent_list, function(parent) parent[[i]])) })

      #single parent
      if (length(parents) == 1) {

        output_vars <- lapply(df_list, function(mat) {
          variance_y(func = edge_func, params = edge_params,
                     trans = edge_trans, input_value = mat, interval_type = interval_type)
        })

        predictions_variance[[node]] <- lapply(output_vars, `[[`, 1)
        parameter_list[[node]]       <- lapply(output_vars, `[[`, 2)

        #mutliparent
      } else if (length(parents) > 1) {

        output_vars <- lapply(df_list, function(mat) {
          variance_y(func = edge_func, params = edge_params,
                     trans = edge_trans, input_value = mat, interval_type = interval_type)
        })

        nsim_local      <- length(output_vars)
        collapsed       <- vector("list", nsim_local)
        param_collapsed <- vector("list", nsim_local)
        var_collapsed   <- vector("list", nsim_local)

        for (k in seq_len(nsim_local)) {

          mats <- output_vars[[k]]

          var_collapsed[[k]] <- mats[[1]]
          param_collapsed[[k]] <- mats[[2]]
        }

        predictions_variance[[node]] <- var_collapsed
        parameter_list[[node]]       <- param_collapsed
      }
    }

    return(list(predictions = list(expected = predictions_expected, variance = predictions_variance),
                parameters = parameter_list,
                root = input_values))
  }

  pred_res <- compute_dag_predictions(object$Structure, object$Parameters, input_values, nsim = nsim)

  return(list(Expected=pred_res$predictions$expected,
              Variance=pred_res$predictions$variance,
              Data=new_data[,given_nodes],
              Structure=object$Structure,
              Parameters=pred_res$parameters,
              nsim=nsim,
              Roots=object$Roots))}
