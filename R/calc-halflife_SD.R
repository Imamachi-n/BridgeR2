#' Calculate RNA half-life for each gene.
#'
#' \code{BridgeRHalfLifeCalcR2Select} returns the dataframe of
#' RNA half-life for each gene.
#'
#' @param inputFile The dataframe of halflife table.
#'
#' @param rawFile The dataframe of RPKM table.
#'
#' @param group The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param inforColumn The number of information columns.
#'
#' @param save Whether to save the output matrix file.
#'
#' @param figSave Whether to save the output fig file.
#'
#' @param outputPrefix The prefix for the name of the output.
#'
#' @export
#'
#' @import data.table ggplot2
#' @importFrom stats sd

CalcHalflifeDeviation <- function(inputFile,
                                  rawFile,
                                  group = c("CTRL_PUM1",
                                            "CTRL_PUM2",
                                            "CTRL_DKD"),
                                  hour = c(0, 1, 2, 4, 8, 12),
                                  save = T,
                                  figSave = F,
                                  inforColumn = 4,
                                  outputPrefix = "BridgeR_7"){

  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.logical(save))
  stopifnot(is.character(outputPrefix))

  # data infor
  time_points <- length(hour)
  group_number <- length(group)
  infor_data_table <- inputFile[, 1:inforColumn, with = F]
  infor_header <- names(infor_data_table)

  # RPKM infor
  RPKM_index <- sapply(1:group_number,
                       function(i) (i - 1) * (time_points + inforColumn) + (inforColumn + 1))
  RPKM_data_table <- rawFile[, RPKM_index, with = F]

  # RNA half-life infor
  halflife_index <- sapply(1:group_number,
                           function(i) i * (time_points + inforColumn + 3))
  half_data_table <- inputFile[, halflife_index, with = F]

  # RNA half-life variance
  calc_halflife_variance <- function(i){
    halflife_list <- suppressWarnings(as.numeric(i))
    halflife_list <- halflife_list[!is.na(halflife_list)]

    # calc RNA half-life variance
    variance_list <- NULL
    if (length(halflife_list) == 0 || length(halflife_list) == 1) {
      variance_list <- c(NA, NA, NA)
    } else {
      variance_list <- halflife_list - mean(halflife_list)
    }
    return(variance_list)
  }

  halflife_variance <- unlist(t(apply(half_data_table, 1, calc_halflife_variance)))
  halflife_variance <- halflife_variance[!is.na(halflife_variance)]
  halflife_variance <- as.data.frame(halflife_variance)

  # density plot fig - RNA half-life variancer
  if (figSave == TRUE) {
    # Fig information
    fig_name <- paste(outputPrefix, "_Histgram_Halflife_variance.png", sep="")
    png(filename = fig_name, width = 600, height = 600)

    # Fig_data: RPKM(mean) vs Half-life(SD)
    p <- drawHalfVariance(halflife_variance)
    suppressWarnings(plot(p))
    dev.off()    # close fig
    plot.new()
  }


  # Calc RNA half-life/RPKM mean & SD
  test_func <- function(vector){
    c(means = mean(as.numeric(vector), na.rm = T),
      sds = sd(as.numeric(vector), na.rm = T) / mean(as.numeric(vector), na.rm = T))
  }

  # RNA half-life mean & SD
  halflife_sd <- as.data.frame(t(apply(half_data_table, 1, test_func)))

  # RPKM mean & SD
  RPKM_sd <- as.data.frame(t(apply(RPKM_data_table, 1, test_func)))

  # all data
  fig_data <- cbind(halflife_sd,RPKM_sd)
  header_label <- c("Half_mean","Half_SD","RPKM_mean","RPKM_SD")
  colnames(fig_data) <- header_label

  # all data collection
  result_data <- cbind(infor_data_table,
                       half_data_table,
                       RPKM_data_table,
                       fig_data)
  half_names <- paste("HalfLife_", group, sep="")
  RPKM_names <- paste("RPKM_", group, sep="")
  colnames(result_data) <- c(infor_header,
                             half_names,
                             RPKM_names,
                             header_label)

  # Fig_data: RPKM(mean) vs Half-life(SD)
  if (figSave == TRUE) {
    fig_name <- paste(outputPrefix, "_RPKM_mean_vs_HalfLife_SD.png", sep="")
    png(filename = fig_name, width = 600, height = 600)
    p <- draw_rpkm_vs_halflife_sd(fig_data)
    suppressWarnings(plot(p))
    dev.off()    # close fig
    plot.new()
  }

  # Fig_data: Half-life(mean) vs Half-life(SD)
  if (figSave == TRUE) {
    fig_name <- paste(outputPrefix, "_HalfLife_mean_vs_HalfLife_SD.png", sep="")
    png(filename = fig_name, width = 600, height = 600)
    p <- draw_halflife_mean_vs_sd(fig_data)
    suppressWarnings(plot(p))
    dev.off()    # close fig
    plot.new()
  }

  if (save == TRUE) {
    write.table(result_data, quote = F, sep = "\t", row.names = F,
                file = paste(outputPrefix, "_halflife_mean_SD.txt", sep=""))
  }

  return(result_data)
}

#' RNA half-life Grubbs test.
#'
#' \code{BridgeRGrubbsTest} returns the dataframe of
#' RNA half-life for each gene.
#'
#' @param controlFile The dataframe of halflife table.
#'
#' @param compFile The dataframe of RPKM table.
#'
#' @param controlGroup The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param inforColumn The number of information columns.
#'
#' @param compIndex The number of information columns.
#'
#' @param save Whether to save the output matrix file.
#'
#' @param outputPrefix The prefix for the name of the output.
#'
#' @export
#'
#' @import data.table ggplot2
#'
#' @importFrom outliers grubbs.test

BridgeRGrubbsTest <- function(controlFile,
                              compFile,
                              hour = c(0,1,2,4,8,12),
                              controlGroup = c("CTRL_PUM1",
                                               "CTRL_PUM2",
                                               "CTRL_DKD"),
                              inforColumn = 4,
                              compIndex = 2,
                              save = T,
                              outputPrefix = "BridgeR_8"){

  # infor data
  time_points <- length(hour)
  infor_data_table <- controlFile[, 1:(inforColumn), with = F]
  infor_header <- names(infor_data_table)
  halflife_header <- c(controlGroup, "Treated")
  grubbs_header <- c("p-value", "Status")

  # Control RNA half-life infor
  half_st <- 1 + inforColumn
  half_ed <- inforColumn + length(controlGroup)
  half_ctrl_data_table <- controlFile[, half_st:half_ed, with = F]

  # treated RNA half-life infor
  halflife_index <- compIndex * (time_points + inforColumn + 3)
  half_treated_data_table <- compFile[, halflife_index, with = F]

  # all half data
  half_all_data_table <- cbind(half_ctrl_data_table,
                               half_treated_data_table)

  # Grabbs test
  grabbs_halflife_test <- function(vector){
    # halflife list prep
    half_list <- as.numeric(vector)
    treated_halflife <- half_list[length(half_list)]
    if (is.na(treated_halflife) || is.nan(treated_halflife)) {
      return(c("NA", "Notest"))
    }
    half_list <- half_list[!is.na(half_list)]
    if (length(half_list) == 0) {
      return(c("NA", "Notest"))
    }

    # Grubbs test

    grubbs_result <- suppressWarnings(grubbs.test(half_list))
    grubbs_alternative <- grubbs_result$alternative
    grubbs_alternative <- gsub("highest value ","",grubbs_alternative)
    grubbs_alternative <- gsub("lowest value ","",grubbs_alternative)
    grubbs_alternative <- as.numeric(gsub(" is an outlier","",grubbs_alternative))

    if(grubbs_alternative == treated_halflife){
      grubbs_pvalue <- grubbs_result$p.value
      return(c(grubbs_pvalue, "Grubbs"))
    }else{
      return(c("NA", "Notest"))
    }
  }

  # testing
  grubbs_table <- as.data.frame(t(apply(half_all_data_table, 1, grabbs_halflife_test)))
  result_table <- cbind(infor_data_table, half_all_data_table, grubbs_table)
  colnames(result_table) <- c(infor_header, halflife_header, grubbs_header)

  if (save == TRUE) {
    write.table(result_table, quote = F, sep = "\t", row.names = F,
                file = paste(outputPrefix, "_halflife_Grubbs_test.txt", sep=""))
  }

  return(result_table)
}

# ggplot2 wrapper
drawHalfVariance <- function(table){
  p <- ggplot(data = table,
              aes(x = halflife_variance))
  p <- p + geom_histogram(fill = "steelblue", color = "white", binwidth = 0.1)
  p <- p + xlim(-5,5)
  p <- p + xlab("RNA half-life variance") + ylab("Count")
  return(p)
}

draw_rpkm_vs_halflife_sd <- function(table){
  p <- ggplot(data = table, aes(x = as.numeric(RPKM_mean),
                                y = as.numeric(Half_SD)))
  p <- p + geom_point(alpha = 0.2)
  p <- p + geom_smooth(method = "loess")
  p <- p + xlab("RPKM mean") + ylab("RNA Half-life SD / RNA Half-life mean")
  p <- p + xlim(0,100) + ylim(0,2)
  return(p)
}

draw_halflife_mean_vs_sd <- function(table){
  p <- ggplot(data = table,
              aes(x = as.numeric(Half_mean),
                  y = as.numeric(Half_SD)))
  p <- p + geom_point(alpha = 0.2)
  p <- p + geom_smooth(method = "loess")
  p <- p + xlab("RNA Half-life mean") + ylab("Half-life SD / Half-life mean")
  p <- p + xlim(0,24) + ylim(0,2)
  return(p)
}


# testing
# library(data.table)
# library(ggplot2)
# library(outliers)
# half_sd_table <- CalcHalflifeDeviation(halflife_table, raw_table,
#                                        outputPrefix = "data/BridgeR_7")
#
# controlFile <- half_sd_table
# compFile <- fread("C:/Users/Naoto/OneDrive/Shiny_app/For_Git/BridgeR2/BridgeR_6_halflife_pvalue_evaluation.txt", header = T)
# grubbs_table <- BridgeRGrubbsTest(half_sd_table,
#                                   compFile = compFile)
