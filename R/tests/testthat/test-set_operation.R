# Test cases ported from py/tests/test_set_operation.py

test_that("bind concatenates two traits", {
  skip_if_offline()

  one <- get_traits(efo_id = "EFO_0009713")
  two <- get_traits(efo_id = "OBA_VT0001321")
  three <- get_traits(efo_id = "EFO_0920093")

  one_and_two <- bind(one, two)
  two_and_three <- bind(two, three)
  expect_equal(length(one_and_two), 2)
  expect_equal(length(two_and_three), 2)
})

test_that("intersect returns common elements", {
  skip_if_offline()

  one <- get_traits(efo_id = "EFO_0009713")
  two <- get_traits(efo_id = "OBA_VT0001321")
  three <- get_traits(efo_id = "EFO_0920093")

  one_and_two <- bind(one, two)
  two_and_three <- bind(two, three)
  result <- intersect(one_and_two, two_and_three)
  expect_equal(length(result), 1)
})

test_that("set_diff returns difference", {
  skip_if_offline()

  one <- get_traits(efo_id = "EFO_0009713")
  two <- get_traits(efo_id = "OBA_VT0001321")
  three <- get_traits(efo_id = "EFO_0920093")

  one_and_two <- bind(one, two)
  two_and_three <- bind(two, three)
  result <- set_diff(one_and_two, two_and_three)
  expect_equal(length(result), 1)
})

test_that("set_xor returns symmetric difference", {
  skip_if_offline()

  one <- get_traits(efo_id = "EFO_0009713")
  two <- get_traits(efo_id = "OBA_VT0001321")
  three <- get_traits(efo_id = "EFO_0920093")

  one_and_two <- bind(one, two)
  two_and_three <- bind(two, three)
  result <- set_xor(one_and_two, two_and_three)
  expect_equal(length(result), 2)
})

test_that("union returns combined unique elements", {
  skip_if_offline()

  one <- get_traits(efo_id = "EFO_0009713")
  two <- get_traits(efo_id = "OBA_VT0001321")
  three <- get_traits(efo_id = "EFO_0920093")

  one_and_two <- bind(one, two)
  two_and_three <- bind(two, three)
  result <- union(one_and_two, two_and_three)
  expect_equal(length(result), 3)
})

test_that("set_equal checks equality", {
  skip_if_offline()

  df1 <- get_traits(efo_id = "EFO_0009713")
  df2 <- get_traits(efo_id = "EFO_0009713")
  expect_true(set_equal(df1, df2))
})
