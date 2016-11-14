#' Calculate relative RPKM expression.
#'
#' \code{BridgeRDataSetFromMatrix} returns the dataframe of
#' the relative RPKM values compared with 0hr.
#'
#' @param inputFile The vector of tab-delimited matrix file.
#'
#' @param group The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param cutoff Cutoff value of RPKM at 0hr.
#'
#' @param cutoffBelow Cutoff value of RPKM at all time points.
#'
#' @param inforColumn The number of information columns.
#'
#' @param save Whether to save the output matrix file.
#'
#' @param outputPrefix The prefix for the name of the output.
#'

BridgeRDataSetFromMatrix <- function(inputFile,
                                     group = c("Control","Knockdown"),
                                     hour = c(0, 1, 2, 4, 8, 12),
                                     cutoff = 0.1,
                                     cutoffBelow = 0.1,
                                     inforColumn = 4,
                                     save = T,
                                     outputPrefix = "BridgeR_1"){

  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.numeric(cutoff))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.logical(save))
  stopifnot(is.character(outputPrefix))

  # prepare files
  time_points <- length(hour)
  input_file_numbers <- length(inputFile)

  input_matrix <- NULL
  for (file in inputFile) {
    if (is.null(input_matrix)) {
      input_matrix <- suppressWarnings(fread(file, header = T))
    }else{
      input_matrix <- cbind(input_matrix,
                            suppressWarnings(fread(file, header = T)))
    }
  }

  # header label
  header_label <- NULL
  hour_label <- NULL

  for (group_index in 1:length(group)) {
    hour_label <- NULL
    for (x in hour) {
      label <- x
      if (x < 10) {
        label <- paste("0",x,sep="")
      }
      hour_label <- append(hour_label,
                           paste("T", label, "_", group_index, sep=""))
    }
    infor_st_ed <- generate_infor_st_ed(group_index, time_points, inforColumn)
    infor_st <- infor_st_ed[1]
    infor_ed <- infor_st_ed[2]
    infor <- colnames(input_matrix)[infor_st:infor_ed]
    header_label <- append(header_label, c(infor, hour_label))
  }

  # Read data matrix
  gene_number <- nrow(input_matrix) # Total number of genes
  sample_size <- length(group)

  # FUNC: calc_rel_exp
  calc_rel_exp <- function(data){
    data_vector <- NULL

    for (sample_index in 1:sample_size) {
      # Infor data
      infor_st_ed <- generate_infor_st_ed(sample_index,
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
      start_time <- exp[1]

      ## IF: RPKM_cutoff < 0.1
      if (start_time <= cutoff) {
        exp <- rep(0, time_points)
      }

      ## Replace the value under cutoff value with 0
      exp <- replace(exp, which(exp < cutoffBelow), 0)

      rel_exp <- exp/start_time
      rel_exp <- replace(rel_exp, which(is.nan(rel_exp)), 0)
      data_vector <- append(data_vector, c(gene_infor, rel_exp))
    }

    return(data_vector)
  }

  output_matrix <- t(apply((input_matrix), 1, calc_rel_exp))
  colnames(output_matrix) <- header_label
  output_matrix <- data.table(output_matrix)

  if (save == T) {
    write.table(output_matrix, quote = F, sep = "\t", row.names = F,
                file = paste(outputPrefix, "_Relative_expression_dataset.txt", sep=""))
  }

  return(output_matrix)
}

# Test
library(data.table)
inputFile <- c("C:/Users/Naoto/OneDrive/Shiny_app/bridger2/data/PUM1_study_siStealth_compatible_genes_RefSeq_result_mRNA.fpkm_table",
               "C:/Users/Naoto/OneDrive/Shiny_app/bridger2/data/PUM1_study_siPUM1_compatible_genes_RefSeq_result_mRNA.fpkm_table")

group <- c("CTRL","PUM1KD")
hour <- c(0,1,2,4,8,12)
test_table <- BridgeRDataSetFromMatrix(inputFile = inputFile,
                                       group = group,
                                       hour = hour,
                                       cutoff = 0.1,
                                       inforColumn = 4)
