# Test cases ported from py/tests/test_get_genomic_contexts.py

# Genomic contexts by rs_id (Thin fetch is exercised first, then Fat)
test_that("get_genomic_contexts by rs_id", {
  skip_if_offline()

  genomic_context_a <- get_genomic_contexts(rs_id = "rs3093017", mode = "Thin")
  expect_s4_class(genomic_context_a, "GenomicContext")
  expect_equal(genomic_context_a@mode, "Thin")

  genomic_context_a <- get_genomic_contexts(rs_id = "rs3093017")
  expect_true(length(genomic_context_a) >= 1)
  expect_true(nrow(genomic_context_a@genomic_contexts) >= 1)
  expect_equal(genomic_context_a@genomic_contexts$rs_id[1], "rs3093017")
})

# Single genomic context by rs_id + genomicContext_id
test_that("get_genomic_contexts by rs_id and genomicContext_id", {
  skip_if_offline()

  genomic_context_b <- get_genomic_contexts(rs_id = "rs3093017",
                                            genomicContext_id = "206217817")

  expect_equal(length(genomic_context_b), 1)
  expect_equal(nrow(genomic_context_b@genomic_contexts), 1)
  expect_equal(genomic_context_b@genomic_contexts$rs_id[1], "rs3093017")
  expect_equal(genomic_context_b@genomic_contexts$genomicContext_id[1],
               "206217817")
  expect_true(nrow(genomic_context_b@genes) >= 1)
  expect_true(nrow(genomic_context_b@locations) >= 1)
})

# Indexing with composite key: Python is 0-based, R is 1-based
test_that("get_genomic_contexts indexing", {
  skip_if_offline()

  genomic_context_b <- get_genomic_contexts(rs_id = "rs3093017",
                                            genomicContext_id = "206217817")

  # assert genomic_context_b[0] == genomic_context_b
  expect_true(genomic_context_b[1] == genomic_context_b)
  # assert genomic_context_b['rs3093017|206217817'] == genomic_context_b
  expect_true(genomic_context_b["rs3093017|206217817"] == genomic_context_b)
})
