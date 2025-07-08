#' Get Bayes Factor
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
get_BF <- function(object, prior_names=NULL){
  bf <- object$Chains_podd[!(object$Chains_podd$predictor) == "NA", ]

  podds_table <- as.data.frame(table(bf$estimate, paste0(bf$predictor, "_", bf$link, "_", bf$group)))

  podds_table <-  tidyr::spread(podds_table, "Var1", "Freq")
  if(!is.null(prior_names)){
    if(length(prior_names) == c(ncol(podds_table)-1)){
      colnames(podds_table)[1:ncol(podds_table)] <- c("Response", prior_names)
    }else{warning("Prior names not the same as number of priors")}}

  pfrac_table      <- podds_table
  pfrac_table[,-1] <- round(pfrac_table[,-1]/rowSums(pfrac_table[,-1]),3)

  num_cols <- pfrac_table[,-1]
  pairs    <- combn(names(num_cols), 2, simplify = FALSE)
  bf_df <- pfrac_table["Var2"]

  for (pair in pairs) {
    col1 <- pair[1]
    col2 <- pair[2]

    bf_name1 <- paste0("BF", col1, "/", col2)
    bf_name2 <- paste0("BF", col2, "/", col1)

    bf_df[[bf_name1]] <- num_cols[[col1]] / num_cols[[col2]]
    bf_df[[bf_name2]] <- num_cols[[col2]] / num_cols[[col1]]
  }


  return(list(bayesfactor=bf_df, simulates=podds_table, fraction=pfrac_table))}
