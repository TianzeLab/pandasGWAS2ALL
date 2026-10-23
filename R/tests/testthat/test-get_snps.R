# Test cases ported from py/tests/test_get_snps.py

# Single SNP by rs_id (Thin fetch is exercised first, then Fat)
test_that("get_snps by rs_id", {
  skip_if_offline()

  snp_a <- get_snps(rs_id = "rs3093017", mode = "Thin")
  expect_s4_class(snp_a, "Snp")
  expect_equal(snp_a@mode, "Thin")

  snp_a <- get_snps(rs_id = "rs3093017")
  expect_equal(length(snp_a), 1)
  expect_equal(nrow(snp_a@snps), 1)
  expect_equal(snp_a@snps$rs_id[1], "rs3093017")
})

# Query filters; each must return at least one SNP
test_that("get_snps by query filters", {
  skip_if_offline()

  snp_d <- get_snps(pubmed_id = "35241825")
  expect_true(length(snp_d) >= 1)

  snp_e <- get_snps(mapped_gene = "CCR6")
  expect_true(length(snp_e) >= 1)

  snp_f <- get_snps(mapped_gene = "CCR6", extended_geneset = TRUE)
  expect_true(length(snp_f) >= 1)
})

# Indexing: Python is 0-based, R is 1-based
test_that("get_snps indexing", {
  skip_if_offline()

  snp_a <- get_snps(rs_id = "rs3093017")
  snp_b <- get_snps(pubmed_id = "35241825")

  # assert snp_a[0] == snp_a
  expect_true(snp_a[1] == snp_a)
  # assert snp_a['rs3093017'] == snp_a
  expect_true(snp_a["rs3093017"] == snp_a)

  # assert len(snp_b[0:2]) == 2
  expect_equal(length(snp_b[1:2]), 2)
  # assert len(snp_b[range(2)]) == 2
  expect_equal(length(snp_b[1:2]), 2)
  # assert snp_b[0] == snp_b[snp_b.raw_data[0]['rs_id']]
  expect_true(
    snp_b[1] == snp_b[snp_b@raw_data[[1]]$rs_id]
  )
})

# Set operations on Snp objects
test_that("get_snps set operations", {
  skip_if_offline()

  snp_a <- get_snps(rs_id = "rs3093017")
  snp_b <- get_snps(pubmed_id = "35241825")

  # assert len(snp_a + snp_b[0]) >= 1
  expect_true(length(snp_a + snp_b[1]) >= 1)
  # assert len(snp_b - snp_b[0]) == len(snp_b) - 1
  expect_equal(length(snp_b - snp_b[1]), length(snp_b) - 1)
  # assert len(snp_b & snp_b[0]) == 1
  expect_equal(length(snp_b & snp_b[1]), 1)
  # assert len(snp_b | snp_b[0]) == len(snp_b)
  expect_equal(length(snp_b | snp_b[1]), length(snp_b))
  # assert len(snp_b ^ snp_b[0]) == len(snp_b) - 1
  expect_equal(length(snp_b ^ snp_b[1]), length(snp_b) - 1)
})

# Equality of independently fetched identical objects
test_that("get_snps equality of repeated fetches", {
  skip_if_offline()

  snp_a <- get_snps(rs_id = "rs3093017")
  snp_c <- get_snps(rs_id = "rs3093017")
  expect_true(snp_c == snp_a)
})
