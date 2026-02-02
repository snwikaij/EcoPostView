seq_ppmn <- function(object, new_data, sequence=NULL, nsim = 1000, level = 0.9,
                     target_ess_frac = 0.2, lambda_grid = exp(seq(log(1e-8), log(1e8), length.out = 500)),
                     mad_constant=1.4826, print_progress=F){

  if(is.null(sequence)){stop("the sequence must be given")}

  batch <- split.data.frame(new_data, sequence)

  for(i in seq_len(length(batch))){

    object <- update_ppmn(object=object, new_data=batch[[i]],
                          nsim=nsim, level=level,
                          target_ess_frac = target_ess_frac,
                          lambda_grid = lambda_grid,
                          mad_constant = mad_constant)}

  return(object)}
