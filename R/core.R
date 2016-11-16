############################################################
#
# BridgeR organization of R files
#
# core.R ......... most of the statistical code
# plots.R ........ all plotting functions
# utils.R ........ utility codes
#



#' BridgeR Core function
#'
#' \code{BridgeRCore} returns RNA half-life for each gene.
#'
#' @param


BridgeRCore <- function(inputFile,
                        inforColumn = 4,
                        group = c("Control","Knockdown"),
                        hour = c(0, 1, 2, 4, 8, 12),
                        RPKMcutoff = 0.1,
                        RelRPKMFig = F,
                        SelectNormFactor = T,
                        CutoffDataPointNumber = 4,
                        CutoffDataPoint1 = c(1, 2),
                        CutoffDataPoint2 = c(8, 12),
                        ThresholdHalfLife = c(8, 12),
                        CutoffRelExp = 0.001,
                        ModelMode = "R2_selection"
){
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


