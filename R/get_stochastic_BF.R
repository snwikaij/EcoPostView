#' Stochastic Bayes Factors
#'
#' @param Object Model object from the meta function
#' @param prior_names names to place above the priors in the table
#'
#' @description This simple function only works when one uses at least two priors
#' in the meta function. It returns a table with the number of simulation from each prior
#' or the fraction of the simulates. Dividing one by the other would then generate the
#' "Stochastic" Bayes factor. It is not a 'true' Bayes Factor, because normaly the prior odds
#' Are treated as fixed. Since there is no focus on hypothesis testing via Bayes Factors this is
#' not a major concern, but it gives a good posterior description of what prior model is
#' most plausible.
#'
#' @importFrom tidyr spread
#'
#' @export
get_stochastic_BF <- function(object, prior_names=NULL){
  bf <- mod$Chains_podd[!(mod$Chains_podd$predictor) == "NA", ]

  podds_table <- as.data.frame(table(bf$estimate, paste0(bf$predictor, "_", bf$link, "_", bf$group)))

  podds_table <-  tidyr::spread(podds_table, "Var1", "Freq")
  if(!is.null(prior_names)){
    if(length(prior_names) == c(ncol(podds_table)-1)){
      colnames(podds_table)[1:ncol(podds_table)] <- c("Response", prior_names)
    }else{warning("Prior names not the same as number of priors")}}

  pfrac_table      <- podds_table
  pfrac_table[,-1] <- round(pfrac_table[,-1]/rowSums(pfrac_table[,-1]),3)

  return(list(simulates=podds_table, fraction=pfrac_table))}
