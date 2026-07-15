#' Displacement test for PPMN
#'
#' @param object A foundational or updated PPMN
#' @param new_data New dataset with all root nodes and least one child node present
#' @param nsim Number of simulations needed for updating
#' @param level Credibility level (default level = 0.9)
#'
#' @export
disp_ppmn <- function(object, new_data, nsim=1000, level=0.9){

  pred1    <- predict_ppmn(object, new_data, nsim = nsim)
  ppmn_100 <- update_ppmn(object, kl_frac = 1, val_data, nsim=nsim)
  pred_100 <- predict_ppmn(ppmn_100, new_data, nsim = nsim)

  #select predicted nodes
  pred_nodes_all <- setdiff(names(pred1$Variance), pred1$Roots)

  #select observed nodes/variables
  obs_nodes <- intersect(colnames(new_data), pred_nodes_all)

  #stack all prediction for standard mod
  pred1_stack <- do.call(rbind,
                         lapply(seq_len(nsim), function(i) {
                           do.call(cbind, lapply(pred1$Variance[obs_nodes], `[[`, i))}))

  #stack all predictions for updated post
  pred100_stack <- do.call(rbind,
                           lapply(seq_len(nsim), function(i) {
                             do.call(cbind, lapply(pred_100$Variance[obs_nodes], `[[`, i))}))

  #substract difference
  diffpred <- pred1_stack-pred100_stack
  #standardize per standard mod
  sdpred1  <- apply(pred1_stack, 2, sd)
  #divide diff so by sd standard model
  stdzdiff <- sweep(diffpred, 2, sdpred1, "/")

  #prb per link
  m_link <- colSums(stdzdiff>0)/nrow(stdzdiff)

  #map
  maxpost            <- function(x){d <- density(x); d$x[which.max(d$y)]}

  #density curves
  forestplot <- do.call(rbind, lapply(seq_len(ncol(stdzdiff)), function(i){
    d <- stdzdiff[,i]
    q <- quantile(d, c((1-level)/2, 1-(1-level)/2))
    data.frame(Variable=obs_nodes[i], GD=m_link[i], median=median(d), map=maxpost(d), mu=mean(d), ll=q[1], ul=q[2])}))

  #create final results
  plot <- ggplot(forestplot, aes(paste0(rownames(forestplot), "\nGD=",round(GD,2)), median)) +
    geom_point()+
    geom_errorbar(aes(ymin = ll, ymax = ul),
                  width = 0, linewidth = 0.5)+
    geom_hline(yintercept = 0, linetype = 2,
               linewidth = .6, colour = "tomato3") +
    labs(x = NULL, y = "Displacement") +
    theme_classic()

  list(plot=plot, results=forestplot)}
