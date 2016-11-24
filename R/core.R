############################################################
#
# BridgeR organization of R files
#
# core.R ................................. Core wrapper code for calculating and comparing RNA half-life
# calc-Relative_RPKM.R ................... calculation of relative RPKM (BRIC-seq dataset)
# calc-Normalization.R ................... Normalization of relative RPKM (BRIC-seq dataset)
# calc-RNA_halflife_3models.R ............ Calculaton of RNA half-life using "3 models" method
# calc-RNA_halflife_R2_selection.R ....... Calculation of RNA half-life using "R2 selection" method
# calc-pvalue_evaluation.R ............... Calculation of p-value for Fold-changes of RNA half-lives between two conditions
# calc-halflife_SD.R ..................... Calculation of RNA half-life SD and Grubbs test
# plots-dataset_checking.R ............... Plotting functions for checking BRIC-seq dataset
# plots-RNA_halflife_comparison.R ........ Plotting functions for RNA half-life comparison
# reporting.R ............................ Reporting of RNA half-life and fitting curve with shiny and plotly
# utils.R ................................ Utility codes
#

#' BridgeR Core function
#'
#' \code{BridgeRCore} returns RNA half-life for each gene.
#'
#' @param inputFile The vector of tab-delimited matrix file.
#'
#' @param inforColumn The number of information columns.
#'
#' @param group The vector of group names.
#'
#' @param hour The vector of time course about BRIC-seq experiment.
#'
#' @param RPKMcutoff Cutoff value of RPKM at 0hr.
#'
#' @param save Whether to save the output matrix file.
#'
#' @param outputPrefix The prefix for the name of the output.
#'
#' @export
#'
#' @import data.table ggplot
#'
#' @importFrom BSDA tsum.test


BridgeRCore <- function(inputFile,
                        inforColumn = 4,
                        group = c("Control","Knockdown"),
                        hour = c(0, 1, 2, 4, 8, 12),
                        RPKMcutoff = 0.1,
                        save = TRUE,
                        outputPrefix = "BridgeR"){
  # check arguments

  # Calculate relative RPKM values compared with 0hr.
  test_table <- BridgeRDataSetFromMatrix(inputFile = inputFile,
                                         group = group,
                                         hour = hour,
                                         cutoff = RPKMcutoff,
                                         inforColumn = inforColumn)

  # Check BRIC-seq datasets
  BridgeRDatasetChecker(inputFile = test_table)

  # Calc Normalization factors
  factor_table <- BridgeRNormalizationFactors(inputFile = test_table)

  # Calc Normalized RPKM values
  normalized_table <- BridgeRNormalization(test_table, factor_table)

  # Calc RNA half-life for each gene
  halflife_table <- BridgeRHalfLifeCalcR2Select(normalized_table)

  return(halflife_table)
}


