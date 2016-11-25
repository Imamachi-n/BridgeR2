context("RNA halflife grubbs test")

test_that("performing Grubbs test", {
  # calc-Relative RPKM.R ##########
  group <- c("CTRL_1","CTRL_2","CTRL_3","KD_1")
  hour <- c(0,1,2,4,8,12)
  test_table <- BridgeRDataSetFromMatrix(inputFile = RNA_halflife_grubbs_test,
                                         group = group,
                                         hour = hour,
                                         cutoff = 0.1,
                                         inforColumn = 4,
                                         save = FALSE)
  raw_table <- test_table[[1]]
  test_table <- test_table[[2]]
  expect_is(raw_table, "data.table")
  expect_is(test_table, "data.table")

  # calc-Normalization.R ###########
  factor_table <- BridgeRNormalizationFactors(test_table,
                                              group = group,
                                              hour = hour,
                                              save = FALSE)
  normalized_table <- BridgeRNormalization(test_table,
                                           factor_table,
                                           group = group,
                                           hour = hour,
                                           save = FALSE)
  expect_is(factor_table, "matrix")
  expect_is(normalized_table, "data.table")

  # calc-RNA_halflife_R2_selection.R ##########
  halflife_table <- BridgeRHalfLifeCalcR2Select(normalized_table,
                                                group = group,
                                                hour = hour,
                                                save = FALSE)
  expect_is(halflife_table, "data.table")

  # calc-halflife_SD.R ##########
  half_sd_table <- CalcHalflifeDeviation(halflife_table,
                                         raw_table,
                                         group = c("CTRL_1",
                                                   "CTRL_2",
                                                   "CTRL_3"),
                                         save = FALSE)
  expect_is(half_sd_table, "data.table")

  grubbs_table <- BridgeRGrubbsTest(half_sd_table,
                                    halflife_table,
                                    compIndex = 4,
                                    controlGroup = c("CTRL_1",
                                                     "CTRL_2",
                                                     "CTRL_3"),
                                    save = FALSE)
  expect_is(grubbs_table, "data.table")

})
