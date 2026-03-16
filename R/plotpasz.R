#' Plot p-values as z-values
#'
#' @param p a vector of numeric p-values
#' @param operator an operator indication < or >
#' @param line_position postion of the vertical line
#' @param bin_limit maximum value on the x-axis for the bins
#' @param bin_width width of the of the bins
#'
#' @returns
#' @export
#'
#' @examples
plotpasz <- function(p, operator=NULL,
                     line_position=1.96,
                     bin_limit=6, bin_width=0.5){

  catp <- function(x, p) {

    if(p == 0){p <- 0.0001}
    if(p == 1){p <- 0.9999}

    if (is.na(x)) {
      return(qnorm(1 - p/2))
    } else if (x == "<") { p <- runif(1, 0, p)
    } else if (x == ">"){ p <- runif(1, p, 1)}

    return(qnorm(1 - p/2))}

  if(!is.null(operator)){

    df <- data.frame(p=p, operator=operator)

    z_list <- list()

    for(j in 1:1000){

      z_vec           <- setNames(cbind.data.frame(df$operator, mapply(catp, df$operator, df$p), NA), c("operator", "z", "bin"))
      rownames(z_vec) <- NULL
      z_vec[,3]       <- cut(as.numeric(z_vec[,2]), seq(0, bin_limit, bin_width))
      z_vec           <- as.data.frame(table(z_vec$bin, z_vec$operator))
      z_list[[j]] <- z_vec}

    long_df           <- do.call(rbind.data.frame, z_list)
    colnames(long_df) <- c("bin", "operator", "count")
    df_bins           <- aggregate(data=long_df, count~bin+operator, mean)

    pos_line <- (line_position-(bin_width/2))/bin_width+1

    plhist <- ggplot(df_bins, aes(x = bin, y = count, fill = operator)) +
      geom_bar(stat = "identity", col="black", width = 1) + scale_fill_grey(start = 0.85, end = 0.15)+
      theme_classic()+xlab("z-value")+geom_vline(xintercept=pos_line, lty=2, col="tomato3", lwd=0.8)+
      theme(axis.text.x = element_text(hjust=1, angle=45))

  }else{

    df <- data.frame(p=p)

    df$z <- qnorm(1 - df$p/2)

    plhist <-   ggplot(df, aes(x = z)) +
      geom_histogram(col = "black", fill = "grey70", binwidth = bin_width, boundary = 0, closed = "left") +
      theme_classic() +
      xlab("z-value") + ylab("count")+
      geom_vline(xintercept = line_position, lty = 2, col = "tomato3", lwd = 0.8) +
      theme(axis.text.x = element_text(hjust = 1, angle = 45))
    }

plhist}
