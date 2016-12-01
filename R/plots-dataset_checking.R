#' BRIC-seq Dataset checker
#'
#' \code{BridgeRDatasetChecker} returns several dataset information
#'
#' @param inputFile Input matrix object.
#'
#' @param group The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param inforColumn The number of information columns.
#'
#' @param percentile Percentile numbers.
#'
#' @param save Whether to save the output fig file.
#'
#' @param outputPrefix The prefix for the name of the output.
#'
#' @export
#'
#' @import data.table ggplot2
#' @importFrom grDevices dev.off png
#' @importFrom graphics layout plot plot.new
#' @importFrom stats quantile

BridgeRDatasetChecker <- function(inputFile,
                                  group = c("Control","Knockdown"),
                                  hour = c(0, 1, 2, 4, 8, 12),
                                  inforColumn = 4,
                                  percentile = c(0.99,
                                                 0.95,
                                                 0.90,
                                                 0.80,
                                                 0.70,
                                                 0.60,
                                                 0.50,
                                                 0.40,
                                                 0.30,
                                                 0.20,
                                                 0.10,
                                                 0.05),
                                  save = T,
                                  outputPrefix = "BridgeR_2_raw"){

  # prepare datasets
  sample_size <- length(group)
  time_points <- length(hour)
  time_label_number <- time_points - 1

  checkData <- BridgeRCheckDataPrep(inputFile,
                                    group,
                                    hour,
                                    inforColumn,
                                    percentile)

  merge_fig_data <- data.table(checkData[[1]])
  merge_fig_percentile_data <- data.table(checkData[[2]])
  setkey(merge_fig_data, "label")
  setkey(merge_fig_percentile_data, "name")
  merge_time_label <- checkData[[3]]

  for (sample_index in 1:sample_size) {
    # prepare datasets
    time_label_st <- time_label_number * (sample_index - 1) + 1
    time_label_ed <- time_label_number * sample_index
    time_label <- merge_time_label[time_label_st:time_label_ed]
    exp_percentile_data <- merge_fig_percentile_data[time_label]
    fig_data <- merge_fig_data[time_label]

    # scattered plot fig - percentile
    # Fig information
    fig_name <- paste(outputPrefix, "_Scattered_percentile_",
                      group[sample_index], ".png", sep="")
    fig_width <- 150 * (time_points - 1)
    if (save == TRUE) {
      png(filename = fig_name, width = fig_width, height = 1200)
    }


    # plotting
    p1 <- BridgeRCheckScatter(exp_percentile_data)
    if (save == TRUE) {
      suppressWarnings(plot(p1))
      dev.off()    # close fig
      plot.new()
    }

    # boxplot plot fig - percentile
    # Fig information
    fig_name <- paste(outputPrefix, "_Boxplot_Rel_RPKM_",
                      group[sample_index],".png",sep="")
    # fig_width <- 150 * (time_points - 1)
    if (save == TRUE) {
      png(filename = fig_name, width = fig_width, height = 1200)
    }

    # plotting
    p2 <- BridgeRCheckboxplot(fig_data)
    if (save == TRUE) {
      suppressWarnings(plot(p2))
      dev.off()    # close fig
      plot.new()
    }

    # density plot fig - relative RNA remaining compared with 0hr
    # Fig information
    fig_name <- paste(outputPrefix, "_density_Rel_RPKM_",
                      group[sample_index],".png",sep="")
    if (save == TRUE) {
      png(filename = fig_name, width = 1300, height = 1000)
    }

    # plotting
    p3 <- BridgeRCheckdensity(fig_data)
    if (save == TRUE) {
      suppressWarnings(plot(p3))
      dev.off()    # close fig
      plot.new()
    }
  }

  # boxplot plot fig for all samples
  # prepare fig_name
  fig_name_func <- function(sample_size,
                            outputPrefix,
                            group,
                            figname){
    for(sample_index in 1:sample_size){
      if(sample_index == 1){
        fig_name <- paste(outputPrefix, "_", figname, "_Rel_RPKM_",
                          group[sample_index], sep="")
      }else{
        fig_name <- paste(fig_name, "_", group[sample_index], sep="")
      }
    }
    fig_name <- paste(fig_name, ".png", sep="")
    return(fig_name)
  }
  fig_name <- fig_name_func(sample_size = sample_size,
                            outputPrefix = outputPrefix,
                            group = group,
                            figname = "Boxplot")

  # Fig information
  fig_width <- 150 * (time_points - 1) * sample_size
  if (save == TRUE) {
    png(filename=fig_name, width = fig_width, height = 1200)
  }

  # plotting
  label_sort <- sort(unique(as.character(merge_fig_data$label)))
  merge_fig_data$label <- factor(merge_fig_data$label,
                                 levels = label_sort)
  p4 <- BridgeRCheckboxplot(merge_fig_data)
  if (save == TRUE) {
    suppressWarnings(plot(p4))
    dev.off()    # close fig
    plot.new()
  }

  # scattered plot fig for all samples
  # prepare fig_name
  fig_name <- fig_name_func(sample_size = sample_size,
                            outputPrefix = outputPrefix,
                            group = group,
                            figname = "scattered")

  # Fig information
  fig_width <- 110 * (time_points - 1) * sample_size
  if (save == TRUE) {
    png(filename=fig_name, width = fig_width, height = 1200)
  }

  # plotting
  label_sort <- sort(unique(as.character(merge_fig_percentile_data$name)))
  merge_fig_percentile_data$name <- factor(merge_fig_percentile_data$name,
                                 levels = label_sort)

  p5 <- BridgeRCheckScatter(merge_fig_percentile_data)
  if (save == TRUE) {
    suppressWarnings(plot(p5))
    dev.off()    # close fig
    plot.new()
  }

  return(list(p1, p2, p3, p4, p5))
}


BridgeRCheckDataPrep <- function(inputFile,
                                 group = c("Control","Knockdown"),
                                 hour = c(0, 1, 2, 4, 8, 12),
                                 inforColumn = 4,
                                 percentile = c(0.99,
                                                0.95,
                                                0.90,
                                                0.80,
                                                0.70,
                                                0.60,
                                                0.50,
                                                0.40,
                                                0.30,
                                                0.20,
                                                0.10,
                                                0.05)){

  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.numeric(percentile) && is.vector(percentile))
  stopifnot(all(percentile >= 0))
  stopifnot(all(percentile <= 1))

  # Calc percentile
  test_q <- function(x,y){
    q_func <- function(vec){
      q_value <- as.vector(quantile(x, prob=vec, na.rm=T))
      return(q_value)
    }
    label_func <- function(vec){
      label <- vec * 100
      if (label < 10){
        return(paste("0", label, "%", sep=""))
      } else {
        return(paste(vec*100, "%", sep=""))
      }
    }
    q_vec <- sapply(percentile, q_func)
    factor_label <- sapply(percentile, label_func)
    label <- rep(y,length(percentile))
    q_table <- data.frame(name=label,q=q_vec,factor=factor_label)
    return(q_table)
  }

  # Main
  time_points <- length(hour)
  sample_size <- length(group)
  test_data <- NULL    # Input data for fig

  # Create boxplot for each sample
  merge_fig_data <- NULL
  merge_fig_percentile_data <- NULL
  merge_time_label <- NULL

  for (sample_index in 1:sample_size) {
    # Infor data
    infor_st_ed <- generate_infor_st_ed(sample_index,
                                        time_points,
                                        inforColumn)
    infor_st <- infor_st_ed[1]
    infor_ed <- infor_st_ed[2]

    # choose 0hr index
    exp_0h <- infor_ed + 1

    # filtered data
    test_data <- inputFile[inputFile[[exp_0h]] == 1,]

    # information & exp_data column
    exp_st <- infor_ed + 1
    exp_ed <- infor_ed + time_points

    # hour label
    hour_label <- generate_hour_label(group,
                                      hour,
                                      sample_index)
    merge_time_label <- append(merge_time_label, hour_label)

    # Prepare exp_data
    exp_st <- exp_st + 1    # Skip 0hr
    exp_data <- test_data[, exp_st:exp_ed, with = F]    # Except 0hr

    exp_percentile_data <- NULL
    time_points_for_fig <- time_points - 1    # Except 0hr
    for (time_index in 1:time_points_for_fig) {
      q_data <- test_q(log10(as.numeric(exp_data[[time_index]])),
                       hour_label[[time_index]])
      if (time_index == 1) {
        exp_percentile_data <- q_data
      } else {
        exp_percentile_data <- rbind(exp_percentile_data, q_data)
      }
    }

    fig_data <- generate_fig_log10_matrix(exp_data,
                                          hour_label)

    if (sample_index == 1) {
      merge_fig_data <- fig_data
      merge_fig_percentile_data <- exp_percentile_data
    }else{
      merge_fig_data <- rbind(merge_fig_data, fig_data)
      merge_fig_percentile_data <- rbind(merge_fig_percentile_data,
                                         exp_percentile_data)
    }

    # Matrix data for plotting
    # - fig_data
    # - exp_percentile_data
    # - merge_fig_data
    # - merge_fig_percentile_data

  }
  return(list(merge_fig_data,
              merge_fig_percentile_data,
              merge_time_label))
}


BridgeRCheckScatter <- function(exp_percentile_data){
  # Fig plotting
  p <- ggplot(data=exp_percentile_data,
              aes_string(x="name", y="q", colour="factor"))
  p <- p + geom_point(size = 5,
                      shape = 19)
  p <- p + xlab("") + ylab("Relative RPKM (Time0 = 1)")
  p <- p + ylim(-1.5,1.5)
  return(p)
}

BridgeRCheckboxplot <- function(fig_data){
  # Fig plotting
  p <- ggplot(data=fig_data,
              aes_string(x="label",y="exp"))
  p <- p + geom_boxplot()
  p <- p + xlab("") + ylab("Relative RPKM (Time0 = 1)")
  p <- p + ylim(-2,2)
  return(p)
}

BridgeRCheckdensity <- function(fig_data){
  # Fig plotting
  p <- ggplot(data=fig_data,
              aes_string(x="exp",colour="label"))
  p <- p + geom_density(size=1.2)
  p <- p + xlab("") + ylab("Relative RPKM (Time0 = 1)")
  p <- p + xlim(-2,2) + ylim(0,7)
  return(p)
}
