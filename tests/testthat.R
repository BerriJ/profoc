library(testthat)
library(profoc)

Sys.setenv("OMP_THREAD_LIMIT" = 3)

test_check("profoc")
