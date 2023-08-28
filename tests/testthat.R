Sys.setenv("OMP_THREAD_LIMIT" = 2)

library(testthat)
library(profoc)


test_check("profoc")
