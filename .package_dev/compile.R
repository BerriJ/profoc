Rcpp::compileAttributes()
devtools::document()
devtools::build()
devtools::check(cran = TRUE)
devtools::test()

# Release:

# Bump Version Number in Description
# Update the Date field in Description
# Update the Date in NEWS.md

# CRAN Release checklist:

devtools::spell_check()
devtools::check_rhub()
devtools::release()

# Merge submitted_to_cran into main, without squashing commits, wihtout extra merge commit
# Create a tag on main

# Potential library unload problems can be solved by specifying
# CXXFLAGS=-fno-gnu-unique in ~/.R/Makevars
# or
# PKG_CXXFLAGS=-fno-gnu-unique in src/Makevars

devtools::load_all()
devtools::test()

knots <- make_knots2(9, deg = deg) + 5
knots <- seq(-0.2, 1.2, 0.1) + 0.3
order <- 3
devtools::load_all()
devtools::load_all()

P_cpp <- penalty(knots, order)[[2]]
P_cpp
D <- diff(diag(6), differences = 3)
P_r <- t(D) %*% D
P_r

round(P_cpp, 13)
round(P_r, 13)
