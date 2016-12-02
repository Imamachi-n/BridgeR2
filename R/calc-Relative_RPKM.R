#' Calculate relative RPKM expression from raw data.
#'
#' \code{BridgeRDataSetFromRaw} returns the dataframe of
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
#' @return data.table object about relative RPKM values.
#'
#' @export
#'
#' @import data.table
#' @importFrom utils write.table

BridgeRDataSetFromRaw <- function(inputFile,
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

  input_matrix <- NULL
  for (file in inputFile) {
    if (is.null(input_matrix)) {
      input_matrix <- suppressWarnings(fread(file, header = T))
    }else{
      input_matrix <- cbind(input_matrix,
                            suppressWarnings(fread(file, header = T)))
    }
  }

  output <- BridgeRDataSetInput(inputFile = input_matrix,
                                group = group,
                                hour = hour,
                                cutoff = cutoff,
                                cutoffBelow = cutoffBelow,
                                inforColumn = inforColumn,
                                save = save,
                                outputPrefix = outputPrefix)
  return(output)
}

#' Calculate relative RPKM expression from data.table format.
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
#' @return data.table object about relative RPKM values.
#'
#' @examples
#' library(data.table)
#' rpkm_matrix <- data.table(gr_id = c(8, 9, 14),
#'                           symbol = c("AAAS", "AACS", "AADAT"),
#'                           accession_id = c("NM_015665", "NM_023928", "NM_182662"),
#'                           locus = c("chr12", "chr12", "chr4"),
#'                           CTRL_1_0h = c(41, 5, 5),
#'                           CTRL_1_1h = c(48, 7, 6),
#'                           CTRL_1_2h = c(56, 10, 6),
#'                           CTRL_1_4h = c(87, 12, 10),
#'                           CTRL_1_8h = c(124, 20, 11),
#'                           CTRL_1_12h = c(185, 22, 15),
#'                           gr_id = c(8, 9, 14),
#'                           symbol = c("AAAS", "AACS", "AADAT"),
#'                           accession_id = c("NM_015665", "NM_023928", "NM_182662"),
#'                           locus = c("chr12", "chr12", "chr4"),
#'                           KD_1_0h = c(21, 10, 3),
#'                           KD_1_1h = c(33, 11, 3),
#'                           KD_1_2h = c(42, 15, 4),
#'                           KD_1_4h = c(60, 20, 5),
#'                           KD_1_8h = c(65, 37, 6),
#'                           KD_1_12h = c(70, 42, 6))
#' group <- c("Control", "Knockdown")
#' hour <- c(0, 1, 2, 4, 8, 12)
#' test_table <- BridgeRDataSetFromMatrix(inputFile = rpkm_matrix,
#'                                        group = group,
#'                                        hour = hour,
#'                                        cutoff = 0.1,
#'                                        inforColumn = 4,
#'                                        save = FALSE)
#'
#' @export
#'
#' @import data.table

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

  output <- BridgeRDataSetInput(inputFile = inputFile,
                                group = group,
                                hour = hour,
                                cutoff = cutoff,
                                cutoffBelow = cutoffBelow,
                                inforColumn = inforColumn,
                                save = save,
                                outputPrefix = outputPrefix)
  return(output)
}

BridgeRDataSetInput <- function(inputFile,
                                group = c("Control","Knockdown"),
                                hour = c(0, 1, 2, 4, 8, 12),
                                cutoff = 0.1,
                                cutoffBelow = 0.1,
                                inforColumn = 4,
                                save = T,
                                outputPrefix = "BridgeR_1"){

  # prepare files
  time_points <- length(hour)
  input_file_numbers <- length(inputFile)

  # input file
  input_matrix <- inputFile

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

  # raw data output
  raw_output_matrix <- input_matrix
  colnames(raw_output_matrix) <- header_label
  raw_output_matrix <- data.table(raw_output_matrix)

  # relative data output
  output_matrix <- t(apply((input_matrix), 1, calc_rel_exp))
  colnames(output_matrix) <- header_label
  output_matrix <- data.table(output_matrix)

  if (save == T) {
    write.table(output_matrix, quote = F, sep = "\t", row.names = F,
                file = paste(outputPrefix, "_Relative_expression_dataset.txt", sep=""))
    write.table(raw_output_matrix, quote = F, sep = "\t", row.names = F,
                file = paste(outputPrefix, "_Raw_expression_dataset.txt", sep=""))
  }

  return(list(raw_output_matrix, output_matrix))
}
