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
#' @param downsamplingFig the factor for downsampling.
#'
#' @param makeFig Whether to save the figure of normalization factor.
#'
#' @param cutoffQuantile cutoff value of quantile.
#'
#' @param figOutputPrefix The prefix for the name of figure output.
#'
#' @param factorOutputPrefix The prefix for the name of factor output.
#'
#' @export
#'
#' @import data.table ggplot2

BridgeRNormalizationFactors <- function(inputFile,
                                        group = c("Control","Knockdown"),
                                        hour = c(0, 1, 2, 4, 8, 12),
                                        inforColumn = 4,
                                        save = T,
                                        YMin = -2,
                                        YMax = 2,
                                        downsamplingFig = 0.2,
                                        makeFig = FALSE,
                                        cutoffQuantile = 0.975,
                                        figOutputPrefix = "BridgeR_3_fig",
                                        factorOutputPrefix = "BridgeR_3"){
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
      select_gene_num <- round(nrow(exp_data) * downsamplingFig)
      test_data <- exp_data[sample(nrow(exp_data), select_gene_num),]
      fig_data <- generate_fig_log10_matrix(test_data,
                                            hour)
      class_flg <- unlist(lapply(1:nrow(test_data),
                                 function(x) rep(x, time_points)))
      color_flg <- rep("black", nrow(test_data))
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
      fig_name <- paste(figOutputPrefix, "lineGraph_Rel_RPKM_",
                        group[group_index],".png",sep="")
      png(filename=fig_name, width = 1200, height = 1200)

      # plotting
      p <- BridgeRCheckLineGraph(fig_data, YMin, YMax)
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


#' Calculate normalization factors from house-keeping genes.
#'
#' \code{BridgeRNormalizationFactorsHK} returns the dataframe of
#' the normalization factors from house-keeping genes.
#'
#' @param inputFile The vector of tab-delimited matrix file.
#'
#' @param group The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param inforColumn The number of information columns.
#'
#' @param inforHKGenesRow The column number of house-keeping gene information.
#'
#' @param HKGenes The vector of house-keeping genes.
#'
#' @param save Whether to save the output matrix file.
#'
#' @param factorOutputPrefix The prefix for the name of factor output.
#'
#' @export
#'
#' @import data.table

BridgeRNormalizationFactorsHK <- function(inputFile,
                                          group = c("Control","Knockdown"),
                                          hour = c(0, 1, 2, 4, 8, 12),
                                          inforColumn = 4,
                                          inforHKGenesRow = "symbol",
                                          HKGenes = c("GAPDH",
                                                      "PGK1",
                                                      "PPIA",
                                                      "ENO1",
                                                      "ATP5B",
                                                      "ALDOA"),
                                          save = T,
                                          factorOutputPrefix = "BridgeR_3"){

  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.character(inforHKGenesRow))
  stopifnot(is.character(HKGenes) && is.vector(HKGenes))
  stopifnot(is.logical(save))
  stopifnot(is.character(factorOutputPrefix))

  # data infor
  group_number <- length(group)
  time_points <- length(hour)

  # key setting
  setkeyv(inputFile, inforHKGenesRow)

  # extract house-keeping gene infor
  hk_input <- inputFile[HKGenes]

  hk_table <- NULL
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

    # Calc norm factor
    exp_data <- hk_input[, exp_st:exp_ed, with = F]
    hk_data <- t(apply(exp_data, 2,
                       function(x) prod(as.numeric(x)) ^ (1/length(x))))
    hk_table <- rbind(hk_table, hk_data)
  }

  # result
  rownames(hk_table) <- group
  if (save == TRUE) {
    write.table(hk_table, quote = F, sep = "\t",
                file = paste(factorOutputPrefix,
                             "_HK_normalization_factor.txt", sep=""))
  }

  return(hk_table)
}


BridgeRCheckLineGraph <- function(fig_data, YMin, YMax){
  p <- ggplot(data = data.table(fig_data),
              aes_string(x = "label",
                         y = "exp",
                         class = "class_flg",
                         colour = "color_flg")
              )
  p <- p + geom_line(# size=0.02,
                     alpha=0.2)
  p <- p + xlim(0,max(as.numeric(as.vector(fig_data$label)))) + ylim(YMin, YMax)
  p <- p + ggtitle("All genes distribution")
  p <- p + xlab("Time course") + ylab("Relative RPKM (Time0 = 1)")
  return(p)
}

#' Calculate the normalized RPKM for BRIC-seq dataset.
#'
#' \code{BridgeRNormalization} returns the dataframe of
#' the normalized RPKM values.
#'
#' @param inputFile The vector of tab-delimited matrix file.
#'
#' @param normFactorFile The vector of tab-delimited normalization factor file.
#'
#' @param group The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param inforColumn The number of information columns.
#'
#' @param save Whether to save the output matrix file.
#'
#' @param outputPrefix The prefix for the name of the output.
#'
#' @export
#'
#' @import data.table

BridgeRNormalization <- function(inputFile,
                                 normFactorFile,
                                 group = c("Control","Knockdown"),
                                 hour = c(0, 1, 2, 4, 8, 12),
                                 inforColumn = 4,
                                 save = T,
                                 outputPrefix = "BridgeR_4"){
  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.logical(save))
  stopifnot(is.character(outputPrefix))

  # prepare files
  group_number <- length(group)
  time_points <- length(hour)

  input_matrix <- inputFile
  norm_matrix <- normFactorFile

  # Calc normalized RPKM
  calc_norm_exp <- function(data){
    data_vector <- NULL
    # Infor data
    infor_st_ed <- generate_infor_st_ed(group_index,
                                        time_points,
                                        inforColumn)
    infor_st <- infor_st_ed[1]
    infor_ed <- infor_st_ed[2]
    gene_infor <- data[infor_st:infor_ed]

    # Exp data
    exp_st <- infor_ed + 1
    exp_ed <- infor_ed + time_points
    exp <- data[exp_st:exp_ed]
    exp <- as.numeric(exp)

    # Normalized
    normalized_exp <- exp / norm_factor
    data_vector <- append(data_vector, c(gene_infor, normalized_exp))
    return(data_vector)
  }

  output_matrix <- NULL
  for (group_index in 1:group_number) {
    # header prep
    colname_st <- 1 + (group_index - 1) * (inforColumn + time_points)
    colname_ed <- group_index * (inforColumn + time_points)
    header_label <- colnames(input_matrix)[colname_st:colname_ed]

    # norm factor prep
    norm_factor <- norm_matrix[group_index,]

    # output result
    result_matrix <- t(apply((input_matrix), 1, calc_norm_exp))
    colnames(result_matrix) <- header_label
    output_matrix <- cbind(output_matrix, result_matrix)
  }
  output_matrix <- data.table(output_matrix)

  if (save == T) {
    write.table(output_matrix, quote = F, sep = "\t", row.names = F,
                file = paste(outputPrefix, "_normalized_expression_dataset.txt", sep=""))
  }

  return(output_matrix)
}
