# Test cases ported from py/tests/test_get_associations.py

# Single association by association_id (Thin fetch is exercised first, then Fat)
test_that("get_associations by association_id", {
  skip_if_offline()

  association_a <- get_associations(association_id = "169106883", mode = "Thin")
  expect_s4_class(association_a, "Association")
  expect_equal(association_a@mode, "Thin")

  association_a <- get_associations(association_id = "169106883")
  expect_equal(length(association_a), 1)
  expect_equal(nrow(association_a@associations), 1)
  expect_equal(association_a@associations$association_id[1], "169106883")
})

# Query filters; each must return at least one association
test_that("get_associations by query filters", {
  skip_if_offline()

  association_b <- get_associations(pubmed_id = "35241825")
  expect_true(length(association_b) >= 1)

  association_c <- get_associations(rs_id = "rs3093017")
  expect_true(length(association_c) >= 1)

  association_d <- get_associations(accession_id = "GCST000854")
  expect_true(length(association_d) >= 1)

  association_e <- get_associations(efo_trait = "celiac disease")
  expect_true(length(association_e) >= 1)

  association_f <- get_associations(efo_id = "MONDO_0005130")
  expect_true(length(association_f) >= 1)

  association_g <- get_associations(mapped_gene = "CCR6")
  expect_true(length(association_g) >= 1)

  association_h <- get_associations(accession_id = "GCST000854",
                                    full_pvalue_set = FALSE)
  expect_true(length(association_h) >= 1)

  association_i <- get_associations(efo_trait = "celiac disease",
                                    show_child_trait = TRUE)
  expect_true(length(association_i) >= 1)

  association_j <- get_associations(mapped_gene = "CCR6",
                                    extended_geneset = TRUE)
  expect_true(length(association_j) >= 1)
})

# Indexing: Python is 0-based, R is 1-based
test_that("get_associations indexing", {
  skip_if_offline()

  association_a <- get_associations(association_id = "169106883")
  association_b <- get_associations(accession_id = "GCST000854")

  # assert association_a[0] == association_a
  expect_true(association_a[1] == association_a)
  # assert association_a['169106883'] == association_a
  expect_true(association_a["169106883"] == association_a)

  # assert len(association_b[0:2]) == 2
  expect_equal(length(association_b[1:2]), 2)
  # assert len(association_b[range(2)]) == 2
  expect_equal(length(association_b[1:2]), 2)
  # assert association_b[0] == association_b[association_b.raw_data[0]['association_id']]
  expect_true(
    association_b[1] ==
      association_b[association_b@raw_data[[1]]$association_id]
  )
})

# Set operations on Association objects
test_that("get_associations set operations", {
  skip_if_offline()

  association_a <- get_associations(association_id = "169106883")
  association_b <- get_associations(accession_id = "GCST000854")

  # assert len(association_a + association_b[0]) >= 1
  expect_true(length(association_a + association_b[1]) >= 1)
  # assert len(association_b - association_b[0]) == len(association_b) - 1
  expect_equal(length(association_b - association_b[1]),
               length(association_b) - 1)
  # assert len(association_b & association_b[0]) == 1
  expect_equal(length(association_b & association_b[1]), 1)
  # assert len(association_b | association_b[0]) == len(association_b)
  expect_equal(length(association_b | association_b[1]),
               length(association_b))
  # assert len(association_b ^ association_b[0]) == len(association_b) - 1
  expect_equal(length(association_b ^ association_b[1]),
               length(association_b) - 1)
})

# Equality of independently fetched identical objects
test_that("get_associations equality of repeated fetches", {
  skip_if_offline()

  association_a <- get_associations(association_id = "169106883")
  association_c <- get_associations(association_id = "169106883")
  expect_true(association_c == association_a)
})
