#' Build the foundational Posterior Predictive Meta-analytic Network (PPMN)
#'
#' @param formula An list of expressions that defines the network structure
#' @param data A list with estimates and standard errors to be synthesized
#' @param txt_size Size of the text with a visualization of the graph
#' @param node_size Diameter of the vertices in the DAG
#' @param node_width Width of the variable vertices in the Bipartite graph
#' @param node_height Height of the variable vertices in the Bipartite graph
#' @param fun_width Width of the vertices representing the edge-functions in the Bipartite graph
#' @param fun_height Height of the vertices representing the edge-functions in the Bipartite graph
#' @param fun_lwd Width of the edge-function vetices
#' @param node_lwd Width of the outlines variable vertices
#' @param arrow_size Size of the arrow lines
#' @param arrow_offset Distance of the arrows from the vertices
#' @param layout_type Type of layout
#' @param circular Circular layout or not (default circular = FALSE)
#'
#' @importFrom igraph V
#'
#' @importFrom ggraph ggraph
#' @importFrom ggraph circle
#' @importFrom ggraph create_layout
#' @importFrom ggraph geom_edge_link
#' @importFrom ggraph geom_node_point
#' @importFrom ggraph geom_node_text
#' @importFrom ggraph geom_node_tile
#' @importFrom ggforce geom_ellipse
#'
#' @importFrom ggplot2 theme_classic
#' @importFrom ggplot2 coord_cartesian
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 aes
#'
#' @importFrom grid arrow
#' @importFrom grid unit
#'
#' @export
build_ppmn   <- function(formula, data, txt_size=3, node_size=15,
                         node_width=0.21, node_height=0.5,
                         fun_width=0.5, fun_height=0.5,
                         fun_lwd=0.4, node_lwd=0.4,
                         arrow_size=2, arrow_offset=5.5,
                         layout_type="auto", circular=F){

  if(missing(formula))stop("provide a formula for the model")
  if(missing(data))stop("provide a data frame with estimates")

  #Check the formulas for "~".
  e_name <- sapply(formula, function(x) x["edge"])
  if(!all(grepl("~", e_name))){stop("no ~ detected in at least one formula.")}

  #Check if the correct edge functions are given
  e_fun <- sapply(formula, function(x) x["fun"])
  if (!all(e_fun %in% c('identity', 'log', 'logit', 'sigmoidal', 'gaussian', 'gompertz', 'class'))){
    stop("Not the correct edge function given in data set, it should be identity, log, logit, sigmoidal, gaussian or class")}

  #match edge formula with the provided data
  matches        <- lapply(formula, function(x){

    edge_formula <- gsub(" ", "", x["edge"])
    edge_data    <- gsub(" ", "", data$edge)
    which(data$`function` == x["fun"] & edge_data == edge_formula)})

  #Check if if there is a mismatch between equations in the formula and data
  if(any(unlist(lapply(matches, function(x) length(x) == 0)))){
    stop(paste0("Edge name(s): ",sapply(formula[unlist(lapply(matches, function(x) length(x) == 0))], function(f) f["edge"]), ", or function(s): ", sapply(formula[unlist(lapply(matches, function(x) length(x) == 0))], function(f) f["fun"]), ", are different between the formula and data provided."))}

  #Estimate model parameters
  param_est      <- list()
  for(i in 1:length(matches)){

    #For each edge match sub divide the data and split per parameter b0, b1 and b2
    sub_df       <- data[matches[[i]],]
    split_df     <- split(sub_df, sub_df$parameter)

    #Extract edge function and the name of the edge
    edge_fun <- formula[[i]]["fun"]
    edge     <- formula[[i]]["edge"]

    #Select the priors based on the names
    priors        <- formula[[i]][grepl("prior", names(formula[[i]]))]

    #Select the truncations based on the names
    truncation    <- formula[[i]][grepl("trunc", names(formula[[i]]))]

    #Number of independent vars for linear model
    length_iv     <- length(unlist(strsplit(unlist(strsplit(edge, "\\~"))[2], "\\+")))

    #Generate minimal a list of possible prior names and truncations
    if(any(edge_fun  %in%  c('identity', 'log', 'logit'))){
      pos_par       <- paste0(c("mu_b", "se_b"), rep(seq(0, length_iv, 1), each=2))
      pos_trunc     <- paste0(c("a_b", "b_b"), rep(seq(0, length_iv, 1), each=2))
    }else if(edge_fun  %in% c('sigmoidal', 'gompertz')){
      pos_par       <- c("mu_b0", "se_b0", "mu_b1", "se_b1", paste0(c("mu_b", "se_b"), rep(seq(2, length_iv+1, 1), each=2)))
      pos_trunc     <- c("a_b0", "b_b0", "a_b1", "b_b1", paste0(c("a_b", "b_b"), rep(seq(2, length_iv+1, 1), each=2)))
    }else if('class' %in% edge_fun){
      pos_par       <- c("mu_b0", "se_b0", "mu_b1", "se_b1", "mu_b2", "se_b2", "mu_b3", "se_b3", "mu_b4", "se_b4")
      pos_trunc     <- c("a_b0", "b_b0", "a_b1", "b_b1", "a_b2", "b_b2", "a_b3", "b_b3", "a_b4", "b_b4")
    }else{
      pos_par       <- c("mu_b0", "se_b0", "mu_b1", "se_b1", "mu_b2", "se_b2")
      pos_trunc     <- c("a_b0", "b_b0", "a_b1", "b_b1", "a_b2", "b_b2")}

    #Based on the possible prior names select the given prior. If no prior is given set prior to mu=0 and se=1000.
    select_prior       <- lapply(pos_par, function(pattern) {
      selection        <- as.numeric(priors[grepl(pattern, names(priors))])

      if (length(selection) == 0) {
        if (grepl("mu", pattern)) {return(0)
        } else if (grepl("se", pattern)) { return(1000)}
      } else {
        return(selection)
      }})

    #Select possible truncations
    select_truncation  <- lapply(pos_trunc, function(pattern) {
      selection   <- as.numeric(truncation[grepl(pattern, names(truncation))])

      if(length(selection) == 0){
        if(grepl("a", pattern)){return(-Inf)
        }else if(grepl("b", pattern)){return(Inf)}
      }else{
        return(selection)
      }})

    #Generate a list to store the estimations for the parameters of the edge function
    total   <- list()

    #For the number of parameter each has mu and se therefore generate a vector with
    #an interval sequence of 2 this is the same for the truncatin a and b
    seq_par <- seq(1, length(select_prior), by = 2)
    for (j in 1:length(seq_par)) {

      prior_mu    <- select_prior[[seq_par[j]]]
      prior_mu_se <- select_prior[[seq_par[j] + 1]]

      a_trunc     <- select_truncation[[seq_par[j]]]
      b_trunc     <- select_truncation[[seq_par[j] + 1]]

      if(j>length(split_df)){stop(sprintf("Not correct number of parameters given in the data for the edge %s with the function %s",
                                          formula[[i]][2], formula[[i]][1]))}

      total[[c("b0","b1","b2","b3","b4")[j]]]  <- abmeta(split_df[[j]]$estimate, split_df[[j]]$error, prior_mu, prior_mu_se, a=a_trunc, b=b_trunc)
    }

    param_est[[edge]] <- Filter(function(x) length(x) > 0, total)
  }

  #DAG building function
  build_dag <- function(formula){

    edges        <- sapply(formula, function(x) x["edge"])
    funcs        <- sapply(formula, function(x) x["fun"])
    names(funcs) <- edges
    trans        <- sapply(formula, function(x) x["trans"])
    names(trans) <- edges

    # graph for predictions
    start_nodes <- sub("(.*)~", "", edges)
    end_nodes   <- sub("~(.*)", "", edges)
    g_pred      <- igraph::graph_from_edgelist(cbind(start_nodes, end_nodes), directed = TRUE)

    # graph for visualisation: ordinary DAG
    visual_edges <- do.call(rbind, lapply(edges, function(e){
      parts      <- strsplit(e, "~")[[1]]
      dep        <- parts[1]
      indep      <- strsplit(parts[2], "\\+")[[1]]
      cbind(indep = indep, dep = dep)}))

    g_vis <- igraph::graph_from_edgelist(visual_edges, directed = TRUE)

    # graph for visualisation: bipartite node-function graph
    bip_edges <- do.call(rbind, lapply(seq_along(edges), function(i){
      parts   <- strsplit(edges[i], "~")[[1]]
      dep     <- parts[1]
      indep   <- strsplit(parts[2], "\\+")[[1]]

      fun_node<- paste0(funcs[i], "_", edges[i])

      rbind(cbind(from = indep,    to = fun_node),
        cbind(from = fun_node, to = dep))}))

    g_bip <- igraph::graph_from_edgelist(bip_edges, directed = TRUE)

    igraph::V(g_bip)$type  <- ifelse(grepl("~", igraph::V(g_bip)$name), "function", "node")
    igraph::V(g_bip)$label <- ifelse(igraph::V(g_bip)$type == "function", sub("_.*", "", igraph::V(g_bip)$name),
    igraph::V(g_bip)$name)

    # sorting graph
    g_exec     <- igraph::graph_from_edgelist(visual_edges, directed = TRUE)
    exec_order <- names(igraph::topo_sort(g_exec, mode = "out"))

    return(list(
      predict_dag = list(
        graph = g_pred,
        functions = funcs,
        edges = edges,
        transformations = trans),

      visual_dag = list(
        graph = g_vis,
        edges = visual_edges),

      bipartite_dag = list(
        graph = g_bip,
        edges = bip_edges),

      exec_dag = list(
        graph = g_exec,
        order = exec_order)))}

  #Create DAG based on previous function
  dag_info        <- build_dag(formula)

  ######################################
  #Generate figure from dag information#
  ######################################

  #Split in names of edges, dependent and independent vars
  split_edges       <- strsplit(e_name, "\\~")
  dependent_vars    <- sapply(split_edges, `[`, 1)
  independent_vars  <- unlist(strsplit(sapply(split_edges, `[`, 2), "\\+"))
  independent_vars  <- unique(independent_vars)

  #True roots are independent vars never appearing as dependent
  true_roots <- setdiff(unique(independent_vars), unique(dependent_vars))

  #Extract visual DAG
  g_vis      <- dag_info$visual_dag$graph

  #Build a full parameter table
  par_tab    <- lapply(seq_along(param_est), function(i) {

    sub      <- param_est[[i]]

    indep    <- strsplit(strsplit(names(param_est)[i], "\\~")[[1]][2], "\\+")[[1]]
    if(e_fun[[i]] %in% c("sigmoidal", 'gompertz', "gaussian")) {
      indep <- c(NA, NA, indep)
    } else if (e_fun[[i]] %in% c("class")){
      indep <- c(NA, indep, NA, indep, NA)
    } else {
      indep <- c(NA, indep)
    }

    cbind(
      edge = names(param_est)[i],
      `edge function` = e_fun[i],
      dependent = dependent_vars[i],
      independent = indep,
      parameters = names(param_est[[i]]),
      do.call(rbind, sub)
    )
  })

  #Combine small table into one
  full_table           <- do.call(rbind.data.frame, par_tab)
  rownames(full_table) <-  NULL

  #extract edges and create code
  edge_cols      <- igraph::as_data_frame(g_vis, what = "edges")
  edge_cols$code <- paste0(edge_cols$from, "~", edge_cols$to)

  #extract parameter and create code
  col_table      <- full_table[!is.na(full_table$independent),c(2:6)]
  col_table$code <- paste0(col_table$dependent, "~", col_table$independent)

  #split the col
  split_list <- split(col_table, col_table$code)
  col_table  <- col_table[!duplicated(col_table$code),]

  class <- unlist(lapply(split_list, nrow))

  for(i in seq_along(class)) {
    if(class[i] > 1 && length(split_list[[i]]$mu) >= 2){
      col_table[col_table$code %in% names(class)[i],"mu"] <-
        as.numeric(split_list[[i]]$mu[1]) - as.numeric(split_list[[i]]$mu[2])
    }
  }

  # Plot DAG
  dag_fig <- ggraph(g_vis, layout = layout_type, circular = circular) +
    ggraph::geom_edge_link(edge_width=0.4,
                   color = "black",
                   arrow = arrow(length = unit(arrow_size, 'mm')),
                   end_cap = circle(arrow_offset, 'mm')) +
    ggraph::geom_node_point(shape = 21, fill = "white", color = "black",
                    size = node_size, stroke = .8) +
    ggraph::geom_node_text(aes(label = name), color = "black", size = txt_size) +
    ggplot2::theme_classic() +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = unit(c(10,10,10,10), "mm"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank())

  #Extract visual Bipartite
  g_bip       <- dag_info$bipartite_dag$graph
  lay         <- ggraph::create_layout(g_bip, layout = layout_type, circular = circular)

  nodes_var   <- lay[lay$type == "node", ]
  nodes_fun   <- lay[lay$type == "function", ]

  #Plot bipartite
  bip_fig <- ggraph::ggraph(g_bip, layout = layout_type, circular = circular) +
    ggraph::geom_edge_link(edge_width=0.4,
                   arrow = arrow(length = unit(arrow_size, "mm")),
                   end_cap = circle(arrow_offset, "mm")) +
    ggforce::geom_ellipse(linewidth=node_lwd,
                 data = nodes_var,
                 aes(x0 = x,
                     y0 = y,
                     a = node_width,
                     b = node_height,
                     angle = 0),
                 fill = "white", colour = "black")+
    ggraph::geom_node_tile(linewidth=fun_lwd,
                   data = nodes_fun,
                   aes(x = x, y = y),
                   width = fun_width,
                   height = fun_height,
                   fill = "white",
                   colour = "black")+
    ggraph::geom_node_text(aes(label = label), color = "black", size = txt_size) +
    ggplot2::theme_classic() +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme(plot.margin = unit(c(10,10,10,10), "mm"),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank())

  theta_names <- unlist(
    Map(function(edge, pars) paste0(edge, "_", names(pars)),
        names(param_est), param_est))

  theta_var <- unlist(
    lapply(param_est, function(e) sapply(e, `[[`, "se")^2))

  Sigma           <- diag(theta_var)
  dimnames(Sigma) <- list(theta_names, theta_names)

  return(list(Plot_as_DAG=dag_fig,
              Plot_as_BIP=bip_fig,
              `Parameter table`=full_table,
              Formula=formula,
              Parameters=param_est,
              Sigma=Sigma,
              Structure=dag_info,
              Data=data,
              Roots=true_roots,
              Update=0))}
