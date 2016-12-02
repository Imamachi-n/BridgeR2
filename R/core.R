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
# data.R ................................. Data information


#' BridgeR basic function for calculating RNA half-life from BRIC-seq data
#'
#' \code{BridgeRCore} is a basic function for calculating RNA half-life
#' BRIC-seq data and a wrapper of the other individual bridger2 functions.
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
#' @param cutoffBelow Cutoff value of RPKM at all time points.
#'
#' @param YMin Y-axis min.
#'
#' @param YMax Y-axis max.
#'
#' @param downsamplingFig the factor for downsampling.
#'
#' @param makeFig Whether to save the figure of normalization factor.
#'
#' @param cutoffQuantile cutoff value of quantile.#' @param save Whether to save the output matrix file.
#'
#' @param inforHKGenesRow The column number of house-keeping gene information.
#'
#' @param HKGenes The vector of house-keeping genes.
#'
#' @param CutoffTimePointNumber The number of minimum time points for calc.
#'
#' @param R2_criteria The cutoff of R2 for R2 selection.
#'
#' @param TimePointRemoval1 The candicate_1 of time point removal.
#'
#' @param TimePointRemoval2 The candicate_2 of time point removal.
#'
#' @param ThresholdHalfLife1 The cutoff of TimePointRemoval1.
#'
#' @param ThresholdHalfLife2 The cutoff of TimePointRemoval2.
#'
#' @param save Whether to save the output matrix file.
#'
#' @param outputPrefix The prefix for the name of the output.
#'
#' @param normalization select "default" (percentile method) or "house_keeping_genes"
#'
#' @param method select "default" (R2 selection/1st-order) or "3models".
#'
#' @return data.table object including RNA half-life, R2 and the selected fitting model.
#'
#' @examples
#' halflife_table <- BridgeRCore(RNA_halflife_comparison[1:30,],
#'                               save = FALSE)
#' halflife_table <- BridgeRCore(RNA_halflife_comparison_HK[177:206],
#'                               save = FALSE,
#'                               normalization = "house_keeping_genes",
#'                               method = "3models")
#'
#' @export
#'
#' @import data.table ggplot2
#'
#' @importFrom BSDA tsum.test


BridgeRCore <- function(inputFile,
                        inforColumn = 4,
                        group = c("Control","Knockdown"),
                        hour = c(0, 1, 2, 4, 8, 12),
                        RPKMcutoff = 0.1,
                        cutoffBelow = 0.1,
                        YMin = -2,
                        YMax = 2,
                        downsamplingFig = 0.2,
                        makeFig = FALSE,
                        cutoffQuantile = 0.975,
                        inforHKGenesRow = "symbol",
                        HKGenes = c("GAPDH",
                                    "PGK1",
                                    "PPIA",
                                    "ENO1",
                                    "ATP5B",
                                    "ALDOA"),
                        CutoffTimePointNumber = 4,
                        R2_criteria = 0.9,
                        TimePointRemoval1 = c(1,2),
                        TimePointRemoval2 = c(8,12),
                        ThresholdHalfLife1 = 3,
                        ThresholdHalfLife2 = 12,
                        save = TRUE,
                        outputPrefix = "BridgeR",
                        normalization = "default",
                        method = "default"){
  # check arguments

  # Calculate relative RPKM values compared with 0hr.
  test_table <- BridgeRDataSetFromMatrix(inputFile = inputFile,
                                         group = group,
                                         hour = hour,
                                         cutoff = RPKMcutoff,
                                         cutoffBelow = cutoffBelow,
                                         inforColumn = inforColumn,
                                         save = save,
                                         outputPrefix = paste(outputPrefix,
                                                              "_1",
                                                              sep=""))
  raw_table <- test_table[[1]]
  test_table <- test_table[[2]]

  # Calc Normalization factors
  factor_table <- NULL
  if (normalization == "default") {
    factor_table <- BridgeRNormalizationFactors(inputFile = test_table,
                                                group = group,
                                                hour = hour,
                                                inforColumn = inforColumn,
                                                save = save,
                                                YMin = YMin,
                                                YMax = YMax,
                                                downsamplingFig = downsamplingFig,
                                                makeFig = makeFig,
                                                cutoffQuantile = cutoffQuantile,
                                                figOutputPrefix = paste(outputPrefix,
                                                                        "_3_fig",
                                                                        sep=""),
                                                factorOutputPrefix = paste(outputPrefix,
                                                                           "_3",
                                                                           sep=""))
  } else if (normalization == "house_keeping_genes") {
    factor_table <- BridgeRNormalizationFactorsHK(inputFile = test_table,
                                                  group = group,
                                                  hour = hour,
                                                  inforColumn = inforColumn,
                                                  inforHKGenesRow = inforHKGenesRow,
                                                  HKGenes = HKGenes,
                                                  save = save,
                                                  factorOutputPrefix = paste(outputPrefix,
                                                                             "_3_HK",
                                                                             sep=""))
  }


  # Calc Normalized RPKM values
  normalized_table <- BridgeRNormalization(test_table,
                                           factor_table,
                                           group = group,
                                           hour = hour,
                                           inforColumn = inforColumn,
                                           save = save,
                                           outputPrefix = paste(outputPrefix,
                                                                "_4",
                                                                sep=""))

  # Calc RNA half-life for each gene
  halflife_table <- NULL
  if (method == "default") {
    halflife_table <- BridgeRHalfLifeCalcR2Select(normalized_table,
                                                  group = group,
                                                  hour = hour,
                                                  inforColumn = inforColumn,
                                                  CutoffTimePointNumber = CutoffTimePointNumber,
                                                  R2_criteria = R2_criteria, # 0.90,
                                                  TimePointRemoval1 = TimePointRemoval1,
                                                  TimePointRemoval2 = TimePointRemoval2,
                                                  ThresholdHalfLife1 = ThresholdHalfLife1,
                                                  ThresholdHalfLife2 = ThresholdHalfLife2,
                                                  save = save,
                                                  outputPrefix = paste(outputPrefix,
                                                                       "_5",
                                                                       sep=""))
  } else if (method == "3models") {
    halflife_table <- BridgeRHalfLifeCalc3models(normalized_table,
                                                 group = group,
                                                 hour = hour,
                                                 inforColumn = inforColumn,
                                                 CutoffTimePointNumber = CutoffTimePointNumber,
                                                 save = save,
                                                 outputPrefix = paste(outputPrefix,
                                                                      "_5",
                                                                      sep=""))
  }

  return(halflife_table)
}
