#' Draw.
#'
#' \code{BridgeRResultChecker} returns  .
#'
#'

BridgeRResultChecker <- function(inputFile,
                                 group = c("Control","Knockdown"),
                                 hour = c(0, 1, 2, 4, 8, 12),
                                 inforColumn = 4,
                                 outputPrefix = "BridgeR_9"){

  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.logical(logScale))
  stopifnot(is.character(outputPrefix))

  # file infor prep
  time_points <- length(hour)
  group_number <- length(group)
  input_matrix <- inputFile

  # RNA half-life infor for scattered plot
  halflife_index <- sapply(1:group_number,
                           function(i) i * (time_points + inforColumn + 3))
  half_data_table <- inputFile[, halflife_index, with = F]
  colnames(half_data_table) <- c("half_1", "half_2")

  # RNA half-life infor for density plot
  half_1_data_table <- inputFile[, halflife_index[1], with = F]
  half_1_data_table_label <- rep(group[1], length(nrow(half_1_data_table)))
  half_1_data_table <- cbind(halflife = half_1_data_table, label = half_1_data_table_label)
  half_2_data_table <- inputFile[, halflife_index[2], with = F]
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
  png(filename = fig_name, width = 600, height = 600)
  p <- draw_halflife_comparison(result, group)
  plot(p)
  dev.off()    # close fig
  plot.new()

  # RNA half-life distribution
  fig_name <- paste(outputPrefix, "_RNA_HalfLife_distribution_",
                    group[1], "_vs_", group[2], ".png", sep="")
  png(filename = fig_name, width = 600, height = 600)
  p <- draw_halflife_distribution(result_density, group)
  plot(p)
  dev.off()    # close fig
  plot.new()

}


draw_halflife_comparison <- function(table, group){
  p <- ggplot(data = table,
              aes(x = log2(as.numeric(half_1)),
                  y = log2(as.numeric(half_2)),
                  colour = factor(flg_list)))
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
              aes(x = as.numeric(halflife),
                  colour = factor(label)))
  p <- p + geom_freqpoly(binwidth = 1, size = 1.0)
  p <- p + xlab("RNA half-life")
  p <- p + ylab("The number of genes")
  p <- p + scale_x_continuous(breaks = seq(0, 24, 2))
  p <- p + theme(axis.title.x = element_text(size=12), axis.title.y = element_text(size=12))
  p <- p + theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))
  return(p)
}

# testing
library(data.table)
p_value_table <- fread("C:/Users/Naoto/OneDrive/Shiny_app/For_Git/BridgeR2/tmp/BridgeR_6_halflife_pvalue_evaluation.txt",
                       header = T)
inputFile <- p_value_table
group = c("Control","Knockdown")
hour = c(0, 1, 2, 4, 8, 12)
inforColumn = 4
logScale = F
outputPrefix = "BridgeR_9"
