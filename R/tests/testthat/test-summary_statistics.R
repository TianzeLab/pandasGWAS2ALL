# Test cases ported from py/tests/test_summary_statistics.py

test_that("search returns correct results by PubMed_id and study_accession_id", {
  skip_if_offline()

  search_DF <- summary_statistics$search(
    PubMed_id = "27918534",
    study_accession_id = "GCST003966"
  )
  expect_equal(nrow(search_DF), 1)
})

test_that("search with online_index returns correct results", {
  skip_if_offline()

  search_DF <- summary_statistics$search(
    PubMed_id = "27918534",
    study_accession_id = "GCST003966",
    online_index = TRUE
  )
  expect_equal(nrow(search_DF), 1)
})

test_that("browser opens FTP directories without error", {
  skip_if_offline()

  search_DF <- summary_statistics$search(PubMed_id = "27918534")
  expect_silent(summary_statistics$browser(search_DF, interactive = FALSE))
})

test_that("download retrieves files to pandasgwas_home", {
  skip_if_offline()

  df <- summary_statistics$search(study_accession_id = "GCST001198")
  expect_silent(summary_statistics$download(df))

  search_DF <- summary_statistics$search(PubMed_id = "20081858")
  expect_silent(summary_statistics$download(search_DF))
})

test_that("parse loads downloaded files into a data.frame", {
  skip_if_offline()

  search_DF <- summary_statistics$search(PubMed_id = "20081858")
  summary_statistics$download(search_DF)
  df <- summary_statistics$parse(search_DF, interactive = FALSE)
  expect_s4_class(df, "data.frame")
  expect_true(nrow(df) > 0)
})
