Rcpp::compileAttributes()
devtools::build()
roxygen2::roxygenize(roclets = "rd")
devtools::check(env_vars = c("NOT_CRAN" = "false"))
devtools::test()

# Release TODO:

# Bump Version Number in Description
# Update the Date field in Description
# Update the Date in NEWS.md
# Create a tagged commit
# Merge develop into main, without squashing commits, wihtout extra merge commit