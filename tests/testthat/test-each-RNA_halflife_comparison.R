# calc-Relative RPKM.R ##########
context("Relative RPKM calculation")

test_that("calculating relative RPKM", {
  group <- c("CTRL","PUM1KD")
  hour <- c(0,1,2,4,8,12)
  test_table <- BridgeRDataSetFromMatrix(inputFile = RNA_halflife_comparison,
                                         group = group,
                                         hour = hour,
                                         cutoff = 0.1,
                                         inforColumn = 4)
  raw_table <- test_table[[1]]
  test_table <- test_table[[2]]
  expect_is(raw_table, "data.table")
  expect_is(test_table, "data.table")
})

# plots-dataset_cheking.R ##########
context("BRIC-esq dataset checking")

test("testing input BRIC-seq dataset", {
  group <- c("CTRL","PUM1KD")
  hour <- c(0,1,2,4,8,12)
  BridgeRDatasetChecker(inputFile = test_table)
})

# calc-Normalization.R ###########
# library(ggplot2)
# factor_table <- BridgeRNormalizationFactors(test_table)
# normalized_table <- BridgeRNormalization(test_table, factor_table)

# plots-dataset_cheking.R ##########
# library(ggplot2)
# BridgeRDatasetChecker(inputFile = normalized_table)

# calc-RNA_halflife_3models.R ##########
# halflife_table <- BridgeRHalfLifeCalc3models(normalized_table)

# calc-RNA_halflife_R2_selection.R ##########
#halflife_table <- BridgeRHalfLifeCalcR2Select(normalized_table)

# calc-pvalue_evaluation.R ##########
# library(BSDA)
# pvalue_table <- BridgeRPvalueEvaluation(halflife_table, calibration = TRUE)

# calc-halflife_SD.R ##########
# library(data.table)
# library(ggplot2)
# library(outliers)
# half_sd_table <- CalcHalflifeDeviation(halflife_table, raw_table,
#                                        outputPrefix = "data/BridgeR_7")
#
# controlFile <- half_sd_table
# grubbs_table <- BridgeRGrubbsTest(half_sd_table,
#                                   compFile = compFile)

# reporting.R ##########
# library(shiny)
# library(shinydashboard)
# library(plotly)
# library(ggplot2)
# library(data.table)
# BridgeReport(pvalue_table)

# plots-RNA_halflife_comparison.R ##########


