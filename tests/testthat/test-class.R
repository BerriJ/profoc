skip_if(debug_mode)

# Test if conline class is exposed:
instance1 <- new(conline)
instance1$trace <- FALSE

# Test if class instance can be passed from R:
expect_false(test_class_input(instance1))

# Test if class instance can be passed to R:
instance2 <- test_class_output()
expect_false(instance2$trace)

# Test if defaults are set correctly:
expect_true(instance2$forget_past_performance == 0)
