#' Calculate RNA half-life for each gene.
#'
#' \code{BridgeRHalfLifeCalcR2Select} returns the dataframe of
#' RNA half-life for each gene.
#'
#' @param inputFile The vector of tab-delimited matrix file.
#'
#' @param group The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param inforColumn The number of information columns.
#'
#' @param CutoffTimePointNumber The number of minimum time points for calc.
#'
#' @param save Whether to save the output matrix file.
#'
#' @param outputPrefix The prefix for the name of the output.
#'

BridgeRHalfLifeCalc3models <- function(inputFile,
                                       group = c("Control","Knockdown"),
                                       hour = c(0, 1, 2, 4, 8, 12),
                                       inforColumn = 4,
                                       CutoffTimePointNumber = 4,
                                       save = T,
                                       outputPrefix = "BridgeR_5"){

  # check arguments
  stopifnot(is.character(group) && is.vector(group))
  stopifnot(is.numeric(hour) && is.vector(hour))
  stopifnot(is.numeric(CutoffTimePointNumber))
  stopifnot(is.numeric(inforColumn))
  stopifnot(is.logical(save))
  stopifnot(is.character(outputPrefix))
}
