Sys.setenv("R_TESTS" = "")

library(testthat)
library(ProteoDisco)

testthat::test_check("ProteoDisco")
