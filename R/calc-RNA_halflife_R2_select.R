#' Calculate the normalized RPKM for BRIC-seq dataset.
#'
#' \code{BridgeRNormalization} returns the dataframe of
#' the normalized RPKM values.
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
#' @param outputPrefix The prefix for the name of the output.
