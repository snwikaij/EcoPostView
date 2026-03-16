seq_ppmn <- function(object, new_data, batch_nr=NULL, nsim = 1000, level = 0.9,
                     update_method = "GBU", target_ess_frac = 0.3,
                     lambda_grid = exp(seq(log(1e-8), log(1e8), length.out = 500)),
                     mad_constant=1.4826, print_progress=F){

  if(is.null(batch_nr)){stop("the batch numbers need to be given")}
  if(!length(unique(batch_nr))>1){stop("at least two unique batches are needed")}

  batch <- split.data.frame(new_data, batch_nr)
  KL    <- vector("list", length(batch))

  for(i in seq_len(length(batch))){

    object <- update_ppmn(object=object, new_data=batch[[i]],
                          update_method = update_method,
                          nsim=nsim, level=level,
                          target_ess_frac = target_ess_frac,
                          lambda_grid = lambda_grid,
                          mad_constant = mad_constant)

    KL[[i]] <- unlist(object$Info$KL_local)}

  return(list(object=object, KL=KL))}
