#' Calculate normalization factors for BRIC-seq datasets.
#'
#' \code{BridgeRNormalizationFactors} returns the dataframe of
#' the normalization factors for BRIC-seq datasets.
#'
#' @param inputFile The vector of tab-delimited matrix file.
#'
#' @param group The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param inforColumn The number of information columns.
#'
#' @param save Whether to save the output matrix file.
#'
#' @param YMin Y-axis min.
#'
#' @param YMax Y-axis max.
#'
#' @param makeFig Whether to save the figure of normalization factor.
#'
#' @param figOutputPrefix The prefix for the name of figure output.
#'

BridgeRNormalizationFactors <- function(inputFile,
                                        group = c("Control","Knockdown"),
                                        hour = c(0, 1, 2, 4, 8, 12),
                                        inforColumn = 4,
                                        save = T,
                                        YMin = -2,
                                        YMax = 2,
                                        makeFig = FALSE,
                                        cutoffQuantile = 0.975,
                                        figOutputPrefix = "BridgeR_3_fig_",
                                        factorOutputPrefix = "BridgeR_3_"){
  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.logical(save))
  stopifnot(is.numeric(YMin))
  stopifnot(is.numeric(YMax))
  stopifnot(is.logical(makeFig))
  stopifnot(is.numeric(cutoffQuantile))
  stopifnot(all(cutoffQuantile >= 0))
  stopifnot(all(cutoffQuantile <= 1))
  stopifnot(is.character(figOutputPrefix))

  group_number <- length(group)
  time_points <- length(hour)

  quantile_table <- NULL
  for (group_index in 1:group_number) {
    # information column
    infor_st_ed <- generate_infor_st_ed(group_index,
                                        time_points,
                                        inforColumn)
    infor_st <- infor_st_ed[1]
    infor_ed <- infor_st_ed[2]
    exp_st <- infor_ed + 1
    exp_ed <- infor_ed + time_points

    # hour label
    hour_label <- generate_hour_label(group,
                                      hour,
                                      group_index)

    # Calc quantile
    exp_data <- inputFile[, exp_st:exp_ed, with = F]
    quantile_calc <- function(vec){
      quantile_data <- quantile(as.numeric(vec),
                                prob = cutoffQuantile, na.rm = T)
      return(quantile_data)
    }
    quantile_data <- apply(exp_data, 2, quantile_calc)
    quantile_table <- rbind(quantile_table, quantile_data)

    # Plotting
    if (makeFig == TRUE) {
      # BRIC-seq data
      fig_data <- generate_fig_log10_matrix(exp_data,
                                            hour)
      class_flg <- unlist(lapply(1:length(exp_data[[1]]), function(x) rep(x, time_points)))
      color_flg <- rep("black", length(exp_data[[1]]))
      fig_data <- cbind(fig_data, class_flg, color_flg)

      # Normalization factor data
      quantile_fig_exp <- generate_fig_log10_matrix(quantile_data,
                                                    hour)
      class_flg <- rep(0, time_points)
      color_flg <- rep("red", time_points)
      quantile_fig <- cbind(quantile_fig_exp,
                            class_flg,
                            color_flg)

      # merge both datasets
      fig_data <- rbind(quantile_fig, fig_data)

      # Fig information
      fig_name <- paste(figOutputPrefix, "_lineGraph_Rel_RPKM_",
                        group[group_index],".png",sep="")
      png(filename=figOutputPrefix, width = 1200, height = 1200)

      # plotting
      p <- BridgeRCheckLineGraph(fig_data)
      plot(p)
      dev.off()    # close fig
      plot.new()
    }

  }

  # result
  rownames(quantile_table) <- group
  if (save == TRUE) {
    write.table(quantile_table, quote = F, sep = "\t",
                file = paste(factorOutputPrefix,
                             "_normalization_factor.txt", sep=""))
  }

  return(quantile_table)
}

BridgeRCheckLineGraph <- function(fig_data){
  p <- ggplot(data = data.table(fig_data2[1:6000,]),
              aes(x = as.numeric(as.vector(label)),
                  y = exp,
                  class = factor(class_flg),
                  colour = color_flg)
              )
  p <- p + geom_line()
                 # size=0.02,
                 # alpha=0.01)
  p <- p + xlim(0,max(as.numeric(as.vector(fig_data$label)))) + ylim(YMin, YMax)
  p <- p + ggtitle("All genes distribution")
  p <- p + xlab("Time course") + ylab("Relative RPKM (Time0 = 1)")
  return(p)
}

# Testing
exp_data <- test_table[, 5:10, with = F]
hour <- c(0,1,2,4,8,12)
time_points <- 6
YMin <- -2
YMax <- 2
fig_data <- generate_fig_log10_matrix(exp_data,
                                      hour)
quantile_table <- BridgeRNormalizationFactors(test_table,
                                              makeFig = T)
