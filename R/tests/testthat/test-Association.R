# Offline unit tests for the Association class.
# py/tests has no dedicated Association test file; class behavior in Python is
# exercised live through test_get_associations.py. This file covers the same
# class behavior offline with crafted records shaped like parsed JSON.

# --- Construction ----------------------------------------------------------

test_that("NULL data yields an empty Fat object", {
  assoc <- Association(NULL, "Fat")
  expect_s4_class(assoc, "Association")
  expect_equal(length(assoc), 0)
  expect_equal(nrow(assoc@associations), 0)
  expect_equal(nrow(assoc@efo_traits), 0)
  expect_equal(nrow(assoc@snp_allele), 0)
  expect_equal(nrow(assoc@bg_efo_traits), 0)
  expect_equal(nrow(assoc@reported_traits), 0)
  expect_equal(nrow(assoc@locations), 0)
  expect_equal(nrow(assoc@mapped_genes), 0)
  expect_equal(nrow(assoc@snp_effect_alleles), 0)
})

test_that("Thin mode stores raw records without building tables", {
  data <- list(
    list(association_id = "12345", pvalue_mantissa = 1.5,
         pvalue_exponent = -8)
  )
  assoc <- Association(data, "Thin")
  expect_equal(assoc@mode, "Thin")
  expect_equal(length(assoc), 1)
  expect_equal(assoc@raw_data[[1]]$association_id, "12345")
  expect_equal(nrow(assoc@associations), 0)
})

test_that("Fat mode builds the main table and nested tables", {
  data <- list(list(
    association_id = "12345",
    pvalue_mantissa = 1.5,
    pvalue_exponent = -8,
    efo_traits = list(list(efo_id = "EFO_0001", efo_trait = "Trait1"))
  ))
  assoc <- Association(data, "Fat")
  expect_equal(assoc@mode, "Fat")
  expect_equal(nrow(assoc@associations), 1)
  expect_equal(nrow(assoc@efo_traits), 1)
})

test_that("invalid mode is rejected", {
  expect_error(Association(NULL, "Invalid"))
})

test_that("non-dict items in data are skipped", {
  data <- list(
    list(association_id = "12345"),
    1, "x", list(1, 2)
  )
  assoc <- Association(data, "Fat")
  expect_equal(length(assoc), 1)
})

# --- Main table columns ----------------------------------------------------

test_that("association_id is kept as character", {
  assoc <- Association(list(list(association_id = 12345)), "Fat")
  expect_equal(assoc@associations$association_id, "12345")
  expect_equal(assoc@raw_data[[1]]$association_id, "12345")
})

test_that("range joins multiple values, single values pass through, NULL is NA", {
  data <- list(
    list(association_id = "1", range = list("chr1", "chr2")),
    list(association_id = "2", range = list("chr3")),
    list(association_id = "3")
  )
  assoc <- Association(data, "Fat")
  expect_equal(assoc@associations$range[1], "chr1; chr2")
  expect_equal(assoc@associations$range[2], "chr3")
  expect_true(is.na(assoc@associations$range[3]))
})

test_that("missing fixed columns are filled with NA, types are preserved", {
  assoc <- Association(list(list(
    association_id = "12345", pvalue_mantissa = 5,
    multi_snp_haplotype = FALSE
  )), "Fat")
  expect_true(is.numeric(assoc@associations$pvalue_mantissa))
  expect_true(is.logical(assoc@associations$multi_snp_haplotype))
  expect_true(is.na(assoc@associations$beta_num))
})

# --- Nested table extraction -----------------------------------------------

test_that("multiple efo_traits and bg_efo_traits expand to rows", {
  data <- list(list(
    association_id = "12345",
    efo_traits = list(
      list(efo_id = "EFO_0001", efo_trait = "Trait1"),
      list(efo_id = "EFO_0002", efo_trait = "Trait2")
    ),
    bg_efo_traits = list(
      list(efo_id = "EFO_0001", efo_trait = "Trait1")
    )
  ))
  assoc <- Association(data, "Fat")
  expect_equal(nrow(assoc@efo_traits), 2)
  expect_equal(nrow(assoc@bg_efo_traits), 1)
  expect_equal(assoc@efo_traits$efo_id, c("EFO_0001", "EFO_0002"))
  expect_true(all(assoc@efo_traits$association_id == "12345"))
})

test_that("single-object snp_allele (named list) yields one row", {
  data <- list(list(
    association_id = "12345",
    snp_allele = list(rs_id = "rs123", effect_allele = "A")
  ))
  assoc <- Association(data, "Fat")
  expect_equal(nrow(assoc@snp_allele), 1)
  expect_equal(assoc@snp_allele$rs_id, "rs123")
  expect_equal(assoc@snp_allele$effect_allele, "A")
  expect_equal(
    colnames(assoc@snp_allele),
    c("association_id", "rs_id", "effect_allele")
  )
})

test_that("scalar reported_trait yields one row", {
  assoc <- Association(list(list(
    association_id = "12345", reported_trait = "Some disease"
  )), "Fat")
  expect_equal(nrow(assoc@reported_traits), 1)
  expect_equal(assoc@reported_traits$reported_trait, "Some disease")
})

test_that("locations become a list column with plural column name", {
  data <- list(list(
    association_id = "12345",
    locations = list(
      list(chromosome_name = "1", chromosome_position = 100L),
      list(chromosome_name = "2", chromosome_position = 200L)
    )
  ))
  assoc <- Association(data, "Fat")
  expect_equal(nrow(assoc@locations), 2)
  expect_equal(colnames(assoc@locations), c("association_id", "locations"))
  expect_true(is.list(assoc@locations$locations))
  expect_equal(
    assoc@locations$locations[[1]],
    list(chromosome_name = "1", chromosome_position = 100L)
  )
})

test_that("mapped_genes flatten into a character column with plural name", {
  data <- list(list(
    association_id = "12345", mapped_genes = list("GENE1", "GENE2")
  ))
  assoc <- Association(data, "Fat")
  expect_equal(nrow(assoc@mapped_genes), 2)
  expect_equal(colnames(assoc@mapped_genes),
               c("association_id", "mapped_genes"))
  expect_equal(assoc@mapped_genes$mapped_genes, c("GENE1", "GENE2"))
})

test_that("snp_effect_allele expands into the plural-named table", {
  data <- list(list(
    association_id = "12345", snp_effect_allele = list("A", "T")
  ))
  assoc <- Association(data, "Fat")
  expect_equal(nrow(assoc@snp_effect_alleles), 2)
  expect_equal(assoc@snp_effect_alleles$snp_effect_allele, c("A", "T"))
})

test_that("NULL nested fields yield zero child rows", {
  assoc <- Association(list(list(association_id = "12345")), "Fat")
  expect_equal(nrow(assoc@efo_traits), 0)
  expect_equal(nrow(assoc@snp_allele), 0)
  expect_equal(nrow(assoc@reported_traits), 0)
  expect_equal(nrow(assoc@locations), 0)
  expect_equal(nrow(assoc@mapped_genes), 0)
})

# --- Indexing --------------------------------------------------------------

test_that("numeric and character indexing return matching subsets", {
  data <- list(
    list(association_id = "12345"),
    list(association_id = "67890"),
    list(association_id = "11111")
  )
  assoc <- Association(data, "Fat")

  expect_equal(length(assoc[1]), 1)
  expect_equal(assoc[1]@raw_data[[1]]$association_id, "12345")
  expect_equal(length(assoc["12345"]), 1)
  expect_equal(assoc["12345"]@raw_data[[1]]$association_id, "12345")
  expect_equal(length(assoc[1:2]), 2)
  expect_equal(length(assoc[c(1, 3)]), 2)
  expect_equal(length(assoc[c("12345", "11111")]), 2)
  expect_equal(length(assoc[]), 3)
})

test_that("unknown id and out-of-bounds position raise errors", {
  assoc <- Association(list(
    list(association_id = "12345"), list(association_id = "67890")
  ), "Fat")
  expect_error(assoc["nope"])
  expect_error(assoc[5])
})

# --- Set operations --------------------------------------------------------

test_that("addition concatenates records", {
  assoc1 <- Association(list(list(association_id = "12345")), "Fat")
  assoc2 <- Association(list(list(association_id = "67890")), "Fat")
  expect_equal(length(assoc1 + assoc2), 2)
})

test_that("subtraction removes matching records", {
  assoc1 <- Association(list(
    list(association_id = "12345"), list(association_id = "67890")
  ), "Fat")
  assoc2 <- Association(list(list(association_id = "12345")), "Fat")
  result <- assoc1 - assoc2
  expect_equal(length(result), 1)
  expect_equal(result@raw_data[[1]]$association_id, "67890")

  assoc_other <- Association(list(list(association_id = "99999")), "Fat")
  unchanged <- assoc1 - assoc_other
  expect_equal(length(unchanged), 2)
  expect_equal(unchanged@raw_data[[1]]$association_id, "12345")
  expect_equal(unchanged@raw_data[[2]]$association_id, "67890")
})

test_that("intersection returns records present in both", {
  assoc1 <- Association(list(
    list(association_id = "12345"), list(association_id = "67890")
  ), "Fat")
  assoc2 <- Association(list(
    list(association_id = "67890"), list(association_id = "11111")
  ), "Fat")
  result <- assoc1 & assoc2
  expect_equal(length(result), 1)
  expect_equal(result@raw_data[[1]]$association_id, "67890")
})

test_that("union keeps unique records and xor keeps symmetric difference", {
  assoc1 <- Association(list(
    list(association_id = "12345"), list(association_id = "67890")
  ), "Fat")
  assoc2 <- Association(list(
    list(association_id = "67890"), list(association_id = "11111")
  ), "Fat")
  expect_equal(length(assoc1 | assoc2), 3)
  expect_equal(length(assoc1 ^ assoc2), 2)
})

test_that("equality reflects raw_data identity and mode", {
  data <- list(list(association_id = "12345"))
  expect_true(Association(data, "Fat") == Association(data, "Fat"))
  expect_false(
    Association(data, "Fat") ==
      Association(list(list(association_id = "67890")), "Fat")
  )
})

test_that("set operations raise errors on mode mismatch", {
  assoc_fat <- Association(list(list(association_id = "12345")), "Fat")
  assoc_thin <- Association(list(list(association_id = "67890")), "Thin")
  expect_error(assoc_fat + assoc_thin)
  expect_error(assoc_fat - assoc_thin)
  expect_error(assoc_fat & assoc_thin)
  expect_error(assoc_fat | assoc_thin)
  expect_error(assoc_fat ^ assoc_thin)
})

# --- Methods: show, summary, names, length ---------------------------------

test_that("show reports mode for Fat and Thin, empty and non-empty", {
  expect_output(show(Association(
    list(list(association_id = "12345")), "Fat"
  )), "fat mode")
  expect_output(show(Association(NULL, "Fat")), "fat mode")
  expect_output(show(Association(
    list(list(association_id = "12345")), "Thin"
  )), "thin mode")
})

test_that("summary lists the object summary", {
  expect_output(summary(Association(
    list(list(association_id = "12345")), "Fat"
  )), "Association Object Summary")
})

test_that("names returns the slot names and length counts raw records", {
  assoc <- Association(list(
    list(association_id = "12345"), list(association_id = "67890")
  ), "Fat")
  expect_equal(names(assoc), c(
    "raw_data", "mode", "associations", "efo_traits", "snp_allele",
    "bg_efo_traits", "reported_traits", "locations", "mapped_genes",
    "snp_effect_alleles"
  ))
  expect_equal(length(assoc), 2)
})
