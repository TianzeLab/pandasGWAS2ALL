# Test cases ported from py/tests/test_get_traits.py

# Single trait by efo_id (Thin fetch is exercised first, then Fat)
test_that("get_traits by efo_id", {
  skip_if_offline()

  trait_a <- get_traits(efo_id = "EFO_0009713", mode = "Thin")
  expect_s4_class(trait_a, "Trait")
  expect_equal(trait_a@mode, "Thin")

  trait_a <- get_traits(efo_id = "EFO_0009713")
  expect_equal(length(trait_a), 1)
  expect_equal(nrow(trait_a@traits), 1)
  expect_equal(trait_a@traits$efo_id[1], "EFO_0009713")
})

# Query filters; each must return at least one trait
test_that("get_traits by query filters", {
  skip_if_offline()

  trait_b <- get_traits(pubmed_id = "35241825")
  expect_true(length(trait_b) >= 1)

  trait_c <- get_traits(efo_trait = "celiac disease")
  expect_true(length(trait_c) >= 1)

  trait_d <- get_traits(uri = "http://www.ebi.ac.uk/efo/EFO_0009713")
  expect_true(length(trait_d) >= 1)

  trait_e <- get_traits(mapped_gene = "CCR6")
  expect_true(length(trait_e) >= 1)

  trait_f <- get_traits(mapped_gene = "CCR6", extended_geneset = TRUE)
  expect_true(length(trait_f) >= 1)
})

# Indexing: Python is 0-based, R is 1-based
test_that("get_traits indexing", {
  skip_if_offline()

  trait_a <- get_traits(efo_id = "EFO_0009713")
  trait_e <- get_traits(mapped_gene = "CCR6")

  # assert trait_a[0] == trait_a
  expect_true(trait_a[1] == trait_a)
  # assert trait_a['EFO_0009713'] == trait_a
  expect_true(trait_a["EFO_0009713"] == trait_a)

  # assert len(trait_e[0:2]) == 2
  expect_equal(length(trait_e[1:2]), 2)
  # assert len(trait_e[range(2)]) == 2
  expect_equal(length(trait_e[1:2]), 2)
  # assert trait_e[0] == trait_e[trait_e.raw_data[0]['efo_id']]
  expect_true(
    trait_e[1] == trait_e[trait_e@raw_data[[1]]$efo_id]
  )
})

# Set operations on Trait objects
test_that("get_traits set operations", {
  skip_if_offline()

  trait_a <- get_traits(efo_id = "EFO_0009713")
  trait_b <- get_traits(pubmed_id = "35241825")

  # assert len(trait_a + trait_b[0]) >= 1
  expect_true(length(trait_a + trait_b[1]) >= 1)
  # assert len(trait_b - trait_b[0]) == len(trait_b) - 1
  expect_equal(length(trait_b - trait_b[1]), length(trait_b) - 1)
  # assert len(trait_b & trait_b[0]) == 1
  expect_equal(length(trait_b & trait_b[1]), 1)
  # assert len(trait_b | trait_b[0]) == len(trait_b)
  expect_equal(length(trait_b | trait_b[1]), length(trait_b))
  # assert len(trait_b ^ trait_b[0]) == len(trait_b) - 1
  expect_equal(length(trait_b ^ trait_b[1]), length(trait_b) - 1)
})

# Equality of independently fetched identical objects
test_that("get_traits equality of repeated fetches", {
  skip_if_offline()

  trait_a <- get_traits(efo_id = "EFO_0009713")
  trait_c <- get_traits(efo_id = "EFO_0009713")
  expect_true(trait_c == trait_a)
})
