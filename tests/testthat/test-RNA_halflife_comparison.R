context("RNA halflife comparison")

test_that("Comparing RNA halflife", {
  # calc-Relative RPKM.R ##########
  group <- c("Control","Knockdown")
  hour <- c(0,1,2,4,8,12)
  test_table <- BridgeRDataSetFromMatrix(inputFile = RNA_halflife_comparison,
                                         group = group,
                                         hour = hour,
                                         cutoff = 0.1,
                                         inforColumn = 4,
                                         save = FALSE)
  raw_table <- test_table[[1]]
  test_table <- test_table[[2]]
  expect_is(raw_table, "data.table")
  expect_is(test_table, "data.table")

  # plots-dataset_cheking.R ##########
  fig_list_raw <- BridgeRDatasetChecker(inputFile = test_table,
                                        save = FALSE)
  expect_is(fig_list_raw, "list")

  # calc-Normalization.R ###########
  factor_table <- BridgeRNormalizationFactors(test_table,
                                              save = FALSE)
  normalized_table <- BridgeRNormalization(test_table, factor_table,
                                           save = FALSE)
  expect_is(factor_table, "matrix")
  expect_is(normalized_table, "data.table")

  # plots-dataset_cheking.R ##########
  fig_list_norm <- BridgeRDatasetChecker(inputFile = normalized_table,
                                         save = FALSE)
  expect_is(fig_list_norm, "list")

  # calc-RNA_halflife_3models.R ##########
  halflife_table <- BridgeRHalfLifeCalc3models(normalized_table,
                                               save = FALSE)
  expect_is(halflife_table, "data.table")

  # calc-RNA_halflife_R2_selection.R ##########
  halflife_table <- BridgeRHalfLifeCalcR2Select(normalized_table,
                                                save = FALSE)
  expect_is(halflife_table, "data.table")

  # calc-pvalue_evaluation.R ##########
  pvalue_table <- BridgeRPvalueEvaluation(halflife_table,
                                          calibration = TRUE,
                                          save = FALSE)
  expect_is(pvalue_table, "data.table")

  # reporting.R ##########
  shiny_test <- BridgeReport(pvalue_table)
  expect_is(shiny_test, "shiny.appobj")

  # plots-RNA_halflife_comparison.R ##########
  result_fig <- BridgeRResultChecker(pvalue_table,
                                     save = FALSE)
  expect_is(result_fig, "list")
})


