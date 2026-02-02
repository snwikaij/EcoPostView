root_est    <- function(object, child=NULL, prior=NULL, print_runtime = T){

  if(is.null(child)){stop("a name for the child variable is not provided")}

  if(!child %in% V(object$Structure$visual_dag$graph)$name){stop("child name not found")}

  if(is.null(prior)){stop("prior data frame is not provided")}

  if(!any(object$Roots %in% colnames(prior))){stop("not all root nodes are provided in the prior data frame")}

  predict_random <- function(model, new_data, child) {
    out <- predict_ppmn(model, nsim = 1, new_data = new_data)
    cbind(yhat = unlist(out$Variance[[child]]), new_data)}

  n_cores <- parallel::detectCores()-1
  cl      <- makeCluster(n_cores)
  registerDoParallel(cl)

  start_time <- Sys.time()

  abc_result <- foreach(
    i = 1:nrow(prior),
    .combine = rbind,
    .packages = c("igraph", "truncnorm", "tmvtnorm", "MASS"),
    .export = c("predict_ppmn")
  ) %dopar% {
    predict_random(object, new_data = prior[i, , drop = FALSE], child = child)}

  end_time <- Sys.time()
  stopCluster(cl)

  tot_time <- end_time-start_time

  if(print_runtime == T){
    print(paste("runtime:", tot_time))}

  return(list(simulations=abc_result, root=object$Roots))}
