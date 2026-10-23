# Test cases ported from py/tests/test_get_genes.py

# Single gene by gene_name (Thin fetch is exercised first, then Fat)
test_that("get_genes by gene_name", {
  skip_if_offline()

  gene_a <- get_genes(gene_name = "HBS1L", mode = "Thin")
  expect_s4_class(gene_a, "Gene")
  expect_equal(gene_a@mode, "Thin")

  gene_a <- get_genes(gene_name = "HBS1L")
  expect_equal(length(gene_a), 1)
  expect_equal(nrow(gene_a@genes), 1)
  expect_equal(gene_a@genes$gene_name[1], "HBS1L")
})

# Indexing: Python is 0-based, R is 1-based
test_that("get_genes indexing", {
  skip_if_offline()

  gene_a <- get_genes(gene_name = "HBS1L")
  gene_b <- get_genes(gene_name = "SAMD5")

  # assert gene_a[0] == gene_a
  expect_true(gene_a[1] == gene_a)
  # assert gene_a['HBS1L'] == gene_a
  expect_true(gene_a["HBS1L"] == gene_a)

  # assert gene_b[0] == gene_b[gene_b.raw_data[0]['gene_name']]
  expect_true(
    gene_b[1] == gene_b[gene_b@raw_data[[1]]$gene_name]
  )
})

# Set operations on Gene objects
test_that("get_genes set operations", {
  skip_if_offline()

  gene_a <- get_genes(gene_name = "HBS1L")
  gene_b <- get_genes(gene_name = "SAMD5")

  # assert len(gene_a + gene_b[0]) >= 1
  expect_true(length(gene_a + gene_b[1]) >= 1)
  # assert len(gene_b - gene_b[0]) == len(gene_b) - 1
  expect_equal(length(gene_b - gene_b[1]), length(gene_b) - 1)
  # assert len(gene_b & gene_b[0]) == 1
  expect_equal(length(gene_b & gene_b[1]), 1)
  # assert len(gene_b | gene_b[0]) == len(gene_b)
  expect_equal(length(gene_b | gene_b[1]), length(gene_b))
  # assert len(gene_b ^ gene_b[0]) == len(gene_b) - 1
  expect_equal(length(gene_b ^ gene_b[1]), length(gene_b) - 1)
})

# Equality of independently fetched identical objects
test_that("get_genes equality of repeated fetches", {
  skip_if_offline()

  gene_a <- get_genes(gene_name = "HBS1L")
  gene_c <- get_genes(gene_name = "HBS1L")
  expect_true(gene_c == gene_a)
})
