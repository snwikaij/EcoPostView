
pred_check  <- function(object, variable=NULL, interval=0.9){

  if(is.null(variable)){stop("variable to display predictive performance needs to be selected")}

  pred_obs_df <- na.omit(data.frame(
    obs=object$Data[[variable]],
    pred=object$Expected[[variable]]))

  predvar <- lapply(1:object$nsim, function(i) object[["Variance"]][[variable]][[i]] - object$Data[[variable]])
  predvar <- do.call(cbind, predvar)

  resid_df <- na.omit(data.frame(predicted=object[["Expected"]][[variable]], residuals=(object$Data[[variable]] - object[["Expected"]][[variable]])))

  cint    <- t(apply(predvar, 1, quantile, probs = c(0.5-(interval/2), 0.5+(interval/2)), na.rm = TRUE))

  diff_summary <- data.frame(
    obs_id = seq_len(nrow(cint)),
    ll = cint[, 1],
    ul = cint[, 2])
  diff_summary$cover = sign(diff_summary$ll) != sign(diff_summary$ul)

  diff_summary   <- na.omit(diff_summary)
  diff_summary$x <- seq(min(pred_obs_df[,2]), max(pred_obs_df[,2]), length.out=nrow(diff_summary))

  rmse_per_type     <- data.frame(mse=mean(abs((pred_obs_df$obs-pred_obs_df$pred)), na.rm=T),
                                  sd=sd(abs(pred_obs_df$obs-pred_obs_df$pred), na.rm=T),
                                  n=nrow(pred_obs_df))

  label <- round(c("MAE"=rmse_per_type$mse,
                   "SD"= rmse_per_type$sd,
                   #"interval coverage"=sum(sign(diff_summary$ll) != sign(diff_summary$ul)) / nrow(diff_summary),
                   "n"=rmse_per_type$n),3)

  plot1 <- ggplot(pred_obs_df, aes(obs, pred))+
    theme_classic()+
    theme(legend.position = "none")

  plot1_vals <- ggplot_build(plot1)
  xlims      <- plot1_vals$layout$panel_params[[1]]$x.range
  ylims      <- plot1_vals$layout$panel_params[[1]]$y.range

  plot1 <- plot1+annotate("text", x=xlims[1]+diff(xlims)*0.7,
                          y=ylims[1]+diff(ylims)*0.25,
                          label=paste0(names(label),"=", label, collapse = "\n"))+
    geom_point(size=3, pch=19, alpha=0.1)+xlab("Observed")+ylab("Predicted")+
    geom_line(data=diff_summary, aes(x, x), lwd=1.2, lty=2, col="tomato3", inherit.aes = F)+
    xlim(xlims)+ylim(ylims)

  plot2 <- ggplot(resid_df, aes(predicted, residuals))+
    geom_point(size=3, pch=19, alpha=0.1)+xlab("Predicted")+ylab("Residuals")+
    geom_hline(yintercept = 0, lty=2, lwd=1.2, col="tomato3")+
    theme_classic()+
    theme(legend.position = "none")

  plot2_vals  <- ggplot_build(plot2)
  xlims2      <- plot2_vals$layout$panel_params[[1]]$x.range
  ylims2      <- plot2_vals$layout$panel_params[[1]]$y.range

  plot2 <- plot2+annotate("text", x=xlims2[1]+diff(xlims2)*0.7,
                          y=ylims2[1]+diff(ylims2)*0.75,
                          label=paste0(names(label),"=", label, collapse = "\n"))

  return(list(`obs vs. pred`=plot1, `pred vs. resid`=plot2))}
