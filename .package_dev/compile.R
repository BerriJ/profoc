Rcpp::compileAttributes()
devtools::document()
devtools::build()
# devtools::check()
devtools::test()

# Release TODO:

# Bump Version Number in Description
# Update the Date field in Description
# Update the Date in NEWS.md
# Create a tagged commit
# Merge develop into main, without squashing commits, wihtout extra merge commit

devtools::load_all()