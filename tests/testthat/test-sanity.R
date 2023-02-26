library(testthat)

test_that("Make sure ANTs can be loaded", {

  cat("\n\n\n")
  print(ants)
  cat("\n\n\n")

  for(cls in c("ants.proxy", "python.builtin.module", "python.builtin.object"
  )) {
    expect_true(inherits(ants, cls))
  }

  for(nm in names(ants)){
    obj <- do.call(`$`, list(ants, nm))
    cat("Checking if [", nm, "] is ants.proxy\n", sep = "")
    expect_true(inherits(obj, "ants.proxy"))
  }

})


