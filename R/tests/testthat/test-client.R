# Test cases ported from py/tests/test_client.py

# Low-level client fetch of a single trait by its full URL
test_that("get_trait returns one record for a trait URL", {
  skip_if_offline()

  a_one <- get_trait(
    "https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits/EFO_0009713"
  )
  expect_type(a_one, "list")
  expect_equal(length(a_one), 1)
})
