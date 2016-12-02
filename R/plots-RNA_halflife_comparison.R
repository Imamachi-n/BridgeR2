#' BRIC-seq result checker
#'
#' \code{BridgeRResultChecker} returns  .
#'
#' @param inputFile The vector of tab-delimited matrix file.
#'
#' @param group The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param inforColumn The number of information columns.
#'
#' @param save Whether to save the output fig file.
#'
#' @param outputPrefix The prefix for the name of the output.
#'
#' @return list object about ggplot2 fig data.
#'
#' @examples
#' library(data.table)
#' normalized_table <- data.table(gr_id = c(8, 9, 14),
#'                                symbol = c("AAAS", "AACS", "AADAT"),
#'                                accession_id = c("NM_015665", "NM_023928", "NM_182662"),
#'                                locus = c("chr12", "chr12", "chr4"),
#'                                CTRL_1_0h = c(1.00, 1.00, 1.00),
#'                                CTRL_1_1h = c(1.00, 0.86, 0.96),
#'                                CTRL_1_2h = c(1.00, 0.96, 0.88),
#'                                CTRL_1_4h = c(1.00, 0.74, 0.85),
#'                                CTRL_1_8h = c(1.00, 0.86, 0.68),
#'                                CTRL_1_12h = c(1.01, 0.65, 0.60),
#'                                gr_id = c(8, 9, 14),
#'                                symbol = c("AAAS", "AACS", "AADAT"),
#'                                accession_id = c("NM_015665", "NM_023928", "NM_182662"),
#'                                locus = c("chr12", "chr12", "chr4"),
#'                                KD_1_0h = c(1.00, 1.00, 1.00),
#'                                KD_1_1h = c(1.01, 0.73, 0.71),
#'                                KD_1_2h = c(1.01, 0.77, 0.69),
#'                                KD_1_4h = c(1.01, 0.72, 0.67),
#'                                KD_1_8h = c(1.01, 0.64, 0.38),
#'                                KD_1_12h = c(1.00, 0.89, 0.63))
#' group <- c("Control", "Knockdown")
#' hour <- c(0, 1, 2, 4, 8, 12)
#' halflife_table <- BridgeRHalfLifeCalcR2Select(normalized_table,
#'                                               group = group,
#'                                               hour = hour,
#'                                               save = FALSE)
#' pvalue_table <- BridgeRPvalueEvaluation(halflife_table,
#'                                         group = group,
#'                                         hour = hour,
#'                                         save = FALSE)
#' result_fig <- BridgeRResultChecker(pvalue_table,
#'                                    save = FALSE)
#'
#' @export
#'
#' @import data.table ggplot2

BridgeRResultChecker <- function(inputFile,
                                 group = c("Control","Knockdown"),
                                 hour = c(0, 1, 2, 4, 8, 12),
                                 inforColumn = 4,
                                 save = T,
                                 outputPrefix = "BridgeR_9"){

  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.logical(save))
  stopifnot(is.character(outputPrefix))

  # file infor prep
  time_points <- length(hour)
  group_number <- length(group)
  input_matrix <- inputFile

  # RNA half-life infor for scattered plot
  halflife_index <- sapply(1:group_number,
                           function(i) i * (time_points + inforColumn + 3))
  half_data_table_pre <- inputFile[, halflife_index, with = F]
  colnames(half_data_table_pre) <- c("half_1", "half_2")
  half_data_table <- log2(half_data_table_pre[, names(half_data_table_pre):=lapply(.SD, as.numeric)])

  # RNA half-life infor for density plot
  half_1_data_table_pre <- inputFile[, halflife_index[1], with = F]
  half_1_data_table <- half_1_data_table_pre[, names(half_1_data_table_pre):=lapply(.SD, as.numeric)]
  half_1_data_table_label <- rep(group[1], length(nrow(half_1_data_table)))
  half_1_data_table <- cbind(halflife = half_1_data_table, label = half_1_data_table_label)
  half_2_data_table_pre <- inputFile[, halflife_index[2], with = F]
  half_2_data_table <- half_2_data_table_pre[, names(half_2_data_table_pre):=lapply(.SD, as.numeric)]
  half_2_data_table_label <- rep(group[2], length(nrow(half_2_data_table)))
  half_2_data_table <- cbind(halflife = half_2_data_table, label = half_2_data_table_label)
  result_density <- rbind(half_1_data_table, half_2_data_table)
  colnames(result_density) <- c("halflife", "label")

  # FC/p-value values
  fold_change_index <- group_number * (time_points + inforColumn + 3) + 1
  p_value_index <- group_number * (time_points + inforColumn + 3) + 2

  FC_p_value <- inputFile[, c(fold_change_index, p_value_index), with = F]

  check_FC_p_value <- function(vector){
    FC <- vector[1]
    p_value <- vector[2]

    flg <- NULL
    if (p_value == "NA" || p_value == "NaN" || is.na(p_value) || is.nan(p_value)) {
      flg <- 0
      return(flg)
    }
    if (FC >= 1 && p_value < 0.05){
      flg <- 1
    } else if (FC <= -1 && p_value < 0.05) {
      flg <- 2
    } else {
      flg <- 0
    }
    return(flg)
  }

  flg_list <- as.factor(apply(FC_p_value, 1, check_FC_p_value))
  result <- cbind(half_data_table, flg_list)

  # RNA half-life comparison
  fig_name <- paste(outputPrefix, "_RNA_HalfLife_comparison_",
                    group[1], "_vs_", group[2], ".png", sep="")
  if (save == TRUE) {
    png(filename = fig_name, width = 600, height = 600)
  }
  p1 <- draw_halflife_comparison(result, group)
  if (save == TRUE) {
    suppressWarnings(plot(p1))
    dev.off()    # close fig
    plot.new()

  }

  # RNA half-life distribution
  fig_name <- paste(outputPrefix, "_RNA_HalfLife_distribution_",
                    group[1], "_vs_", group[2], ".png", sep="")
  if (save == TRUE) {
    png(filename = fig_name, width = 600, height = 600)
  }
  p2 <- draw_halflife_distribution(result_density, group)
  if (save == TRUE) {
    suppressWarnings(plot(p2))
    dev.off()    # close fig
    plot.new()
  }

  return(list(p1, p2))
}


draw_halflife_comparison <- function(table, group){
  p <- ggplot(data = table,
              aes_string(x = "half_1",
                         y = "half_2",
                         colour = "flg_list"))
  p <- p + geom_point(alpha = 0.5)
  p <- p + scale_color_manual(values = c("gray", "red", "blue"))
  p <- p + xlab(paste("log2(RNA half-life [", group[1], "])", sep = ""))
  p <- p + ylab(paste("log2(RNA half-life [", group[2], "])", sep = ""))
  p <- p + xlim(0, 5) + ylim(0, 5)
  p <- p + theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  p <- p + theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))
  p <- p + theme(legend.position = "NONE")
  return(p)
}

draw_halflife_distribution <- function(table, group){
  p <- ggplot(data = table,
              aes_string(x = "halflife",
                         colour = "label"))
  p <- p + geom_freqpoly(binwidth = 1, size = 1.0)
  p <- p + xlab("RNA half-life")
  p <- p + ylab("The number of genes")
  p <- p + scale_x_continuous(breaks = seq(0, 24, 2))
  p <- p + theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  p <- p + theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))
  return(p)
}
