context("RNA halflife comparison")

test_that("Comparing RNA halflife", {
  # calc-Relative RPKM.R ##########
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

  # plots-dataset_cheking.R ##########
  BridgeRDatasetChecker(inputFile = test_table)

  # calc-Normalization.R ###########
  factor_table <- BridgeRNormalizationFactors(test_table)
  normalized_table <- BridgeRNormalization(test_table, factor_table)
  expect_is(factor_table, "matrix")
  expect_is(normalized_table, "data.table")

  # plots-dataset_cheking.R ##########
  BridgeRDatasetChecker(inputFile = normalized_table)

  # calc-RNA_halflife_3models.R ##########
  # halflife_table <- BridgeRHalfLifeCalc3models(normalized_table)
  # expect_is(halflife_table, "data.table")

  # calc-RNA_halflife_R2_selection.R ##########
  halflife_table <- BridgeRHalfLifeCalcR2Select(normalized_table)
  expect_is(halflife_table, "data.table")

  # calc-pvalue_evaluation.R ##########
  # pvalue_table <- BridgeRPvalueEvaluation(halflife_table, calibration = TRUE)
})



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


