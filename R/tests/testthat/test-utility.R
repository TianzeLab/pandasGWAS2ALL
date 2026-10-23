# Test cases ported from py/tests/test_utility.py

test_that("is_API_available returns TRUE when API is reachable", {
  skip_if_offline()

  is_available <- is_API_available()
  expect_true(is_available)
})
