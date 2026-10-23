# Test cases ported from py/tests/test_get_loci.py

# Loci by association_id (Thin fetch is exercised first, then Fat)
test_that("get_loci by association_id", {
  skip_if_offline()

  locus_a <- get_loci(association_id = "5860")
  locus_a <- get_loci(association_id = "169106883", mode = "Thin")
  expect_s4_class(locus_a, "Locus")
  expect_equal(locus_a@mode, "Thin")

  locus_a <- get_loci(association_id = "169106883")
  expect_true(length(locus_a) >= 1)
  expect_true(nrow(locus_a@loci) >= 1)
  expect_equal(locus_a@loci$association_id[1], "169106883")
})

# Single locus by association_id + locus_id
test_that("get_loci by association_id and locus_id", {
  skip_if_offline()

  locus_b <- get_loci(association_id = "169106883", locus_id = "169106882")

  expect_equal(length(locus_b), 1)
  expect_equal(nrow(locus_b@loci), 1)
  expect_equal(locus_b@loci$association_id[1], "169106883")
  expect_equal(locus_b@loci$locus_id[1], "169106882")
  expect_true(nrow(locus_b@strongest_risk_alleles) >= 1)
})

# Indexing with composite key: Python is 0-based, R is 1-based
test_that("get_loci indexing", {
  skip_if_offline()

  locus_b <- get_loci(association_id = "169106883", locus_id = "169106882")

  # assert locus_b[0] == locus_b
  expect_true(locus_b[1] == locus_b)
  # assert locus_b['169106883|169106882'] == locus_b
  expect_true(locus_b["169106883|169106882"] == locus_b)
})
