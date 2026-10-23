# Test cases ported from py/tests/test_get_ancestries.py

# Ancestries by accession_id (Thin fetch is exercised first, then Fat)
test_that("get_ancestries by accession_id", {
  skip_if_offline()

  ancestry_a <- get_ancestries(accession_id = "GCST000854", mode = "Thin")
  expect_s4_class(ancestry_a, "Ancestry")
  expect_equal(ancestry_a@mode, "Thin")

  ancestry_a <- get_ancestries(accession_id = "GCST000854")
  expect_true(length(ancestry_a) >= 1)
  expect_true(nrow(ancestry_a@ancestries) >= 1)
  expect_equal(ancestry_a@ancestries$accession_id[1], "GCST000854")
})

# Single ancestry by accession_id + ancestry_id
test_that("get_ancestries by accession_id and ancestry_id", {
  skip_if_offline()

  ancestry_b <- get_ancestries(accession_id = "GCST000853", ancestry_id = "7529")

  expect_equal(length(ancestry_b), 1)
  expect_equal(nrow(ancestry_b@ancestries), 1)
  expect_equal(ancestry_b@ancestries$accession_id[1], "GCST000853")
  expect_equal(ancestry_b@ancestries$ancestry_id[1], "7529")
})

# Indexing with composite key: Python is 0-based, R is 1-based
test_that("get_ancestries indexing", {
  skip_if_offline()

  ancestry_b <- get_ancestries(accession_id = "GCST000853", ancestry_id = "7529")

  # assert ancestry_b[0] == ancestry_b
  expect_true(ancestry_b[1] == ancestry_b)
  # assert ancestry_b['GCST000853|7529'] == ancestry_b
  expect_true(ancestry_b["GCST000853|7529"] == ancestry_b)
})
