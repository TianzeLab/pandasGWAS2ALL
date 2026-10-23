# Test cases ported from py/tests/test_get_studies.py

# Single study by accession_id (Thin fetch is exercised first, then Fat)
test_that("get_studies by accession_id", {
  skip_if_offline()

  study_a <- get_studies(accession_id = "GCST000854", mode = "Thin")
  expect_s4_class(study_a, "Study")
  expect_equal(study_a@mode, "Thin")

  study_a <- get_studies(accession_id = "GCST000854")
  expect_equal(length(study_a), 1)
  expect_equal(nrow(study_a@studies), 1)
  expect_equal(study_a@studies$accession_id[1], "GCST000854")
})

# Query filters; each must return at least one study
test_that("get_studies by query filters", {
  skip_if_offline()

  study_b <- get_studies(pubmed_id = "35241825")
  expect_true(length(study_b) >= 1)

  study_c <- get_studies(disease_trait = "Early-onset Parkinson's disease")
  expect_true(length(study_c) >= 1)

  study_d <- get_studies(accession_id = "GCST000854", full_pvalue_set = FALSE)
  expect_true(length(study_d) >= 1)

  study_e <- get_studies(efo_id = "MONDO_0005130")
  expect_true(length(study_e) >= 1)

  study_f <- get_studies(efo_trait = "celiac disease")
  expect_true(length(study_f) >= 1)

  study_g <- get_studies(cohort = "BioImage")
  expect_true(length(study_g) >= 1)

  study_k <- get_studies(efo_id = "MONDO_0005130", show_child_trait = TRUE)
  expect_true(length(study_k) >= 1)

  study_l <- get_studies(mapped_gene = "CCR6")
  expect_true(length(study_l) >= 1)

  study_m <- get_studies(mapped_gene = "CCR6", extended_geneset = TRUE)
  expect_true(length(study_m) >= 1)
})

# Indexing: Python is 0-based, R is 1-based
test_that("get_studies indexing", {
  skip_if_offline()

  study_a <- get_studies(accession_id = "GCST000854")
  study_b <- get_studies(pubmed_id = "35241825")

  # assert study_a[0] == study_a
  expect_true(study_a[1] == study_a)
  # assert study_a['GCST000854'] == study_a
  expect_true(study_a["GCST000854"] == study_a)

  # assert len(study_b[0:2]) == 2
  expect_equal(length(study_b[1:2]), 2)
  # assert len(study_b[range(2)]) == 2
  expect_equal(length(study_b[1:2]), 2)
  # assert study_b[0] == study_b[study_b.raw_data[0]['accession_id']]
  expect_true(
    study_b[1] == study_b[study_b@raw_data[[1]]$accession_id]
  )
})

# Set operations on Study objects
test_that("get_studies set operations", {
  skip_if_offline()

  study_a <- get_studies(accession_id = "GCST000854")
  study_b <- get_studies(pubmed_id = "35241825")

  # assert len(study_a + study_b[0]) >= 1
  expect_true(length(study_a + study_b[1]) >= 1)
  # assert len(study_b - study_b[0]) == len(study_b) - 1
  expect_equal(length(study_b - study_b[1]), length(study_b) - 1)
  # assert len(study_b & study_b[0]) == 1
  expect_equal(length(study_b & study_b[1]), 1)
  # assert len(study_b | study_b[0]) == len(study_b)
  expect_equal(length(study_b | study_b[1]), length(study_b))
  # assert len(study_b ^ study_b[0]) == len(study_b) - 1
  expect_equal(length(study_b ^ study_b[1]), length(study_b) - 1)
})

# Equality of independently fetched identical objects
test_that("get_studies equality of repeated fetches", {
  skip_if_offline()

  study_a <- get_studies(accession_id = "GCST000854")
  study_c <- get_studies(accession_id = "GCST000854")
  expect_true(study_c == study_a)
})
