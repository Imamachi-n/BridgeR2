context("BridgeRCore testing1")

test_that("testing BridgeRCore", {
  halflife_table <- BridgeRCore(RNA_halflife_comparison,
                                save = FALSE)
  expect_is(halflife_table, "data.table")
})

context("BridgeRCore testing2")

test_that("testing BridgeRCore 2", {
  halflife_table <- BridgeRCore(RNA_halflife_comparison_HK,
                                save = FALSE,
                                normalization = "house_keeping_genes",
                                method = "3models")
  expect_is(halflife_table, "data.table")
})
