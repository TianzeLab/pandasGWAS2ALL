# Test cases ported from py/tests/test_get_publications.py

# Single publication by pubmed_id (Thin fetch is exercised first, then Fat)
test_that("get_publications by pubmed_id", {
  skip_if_offline()

  publication_a <- get_publications(pubmed_id = "35241825", mode = "Thin")
  expect_s4_class(publication_a, "Publication")
  expect_equal(publication_a@mode, "Thin")

  publication_a <- get_publications(pubmed_id = "35241825")
  expect_equal(length(publication_a), 1)
  expect_equal(nrow(publication_a@publications), 1)
  expect_equal(publication_a@publications$pubmed_id[1], "35241825")
})

# Query filter by title
test_that("get_publications by title", {
  skip_if_offline()

  publication_b <- get_publications(title = "Genome-wide association study")
  expect_true(length(publication_b) >= 1)
})

# Indexing: Python is 0-based, R is 1-based
test_that("get_publications indexing", {
  skip_if_offline()

  publication_a <- get_publications(pubmed_id = "35241825")
  publication_b <- get_publications(title = "Genome-wide association study")

  # assert publication_a[0] == publication_a
  expect_true(publication_a[1] == publication_a)
  # assert publication_a['35241825'] == publication_a
  expect_true(publication_a["35241825"] == publication_a)

  # assert len(publication_b[0:2]) == 2
  expect_equal(length(publication_b[1:2]), 2)
  # assert len(publication_b[range(2)]) == 2
  expect_equal(length(publication_b[1:2]), 2)
  # assert publication_b[0] == publication_b[publication_b.raw_data[0]['pubmed_id']]
  expect_true(
    publication_b[1] == publication_b[publication_b@raw_data[[1]]$pubmed_id]
  )
})

# Set operations on Publication objects
test_that("get_publications set operations", {
  skip_if_offline()

  publication_a <- get_publications(pubmed_id = "35241825")
  publication_b <- get_publications(title = "Genome-wide association study")

  # assert len(publication_a + publication_b[0]) >= 1
  expect_true(length(publication_a + publication_b[1]) >= 1)
  # assert len(publication_b - publication_b[0]) == len(publication_b) - 1
  expect_equal(length(publication_b - publication_b[1]),
               length(publication_b) - 1)
  # assert len(publication_b & publication_b[0]) == 1
  expect_equal(length(publication_b & publication_b[1]), 1)
  # assert len(publication_b | publication_b[0]) == len(publication_b)
  expect_equal(length(publication_b | publication_b[1]),
               length(publication_b))
  # assert len(publication_b ^ publication_b[0]) == len(publication_b) - 1
  expect_equal(length(publication_b ^ publication_b[1]),
               length(publication_b) - 1)
})

# Equality of independently fetched identical objects
test_that("get_publications equality of repeated fetches", {
  skip_if_offline()

  publication_a <- get_publications(pubmed_id = "35241825")
  publication_d <- get_publications(pubmed_id = "35241825")
  expect_true(publication_d == publication_a)
})
