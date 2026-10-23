#' Locus class
#'
#' An S4 class that stores data of type Locus fetched from the GWAS API.
#' In Fat mode it holds 5 data.frames with hierarchical dependencies; in
#' Thin mode it only stores the raw list data.
#'
#' @slot raw_data list. Converted from obtained JSON data
#' @slot mode character. Fat or Thin
#' @slot loci data.frame
#' @slot strongest_risk_alleles data.frame
#' @slot author_reported_genes data.frame
#' @slot ensembl_gene_ids data.frame
#' @slot entrez_gene_ids data.frame
#'
#' @importFrom methods new
#'
#' @export
#'
setClass(
  "Locus",
  slots = c(
    raw_data = "list",
    mode = "character",
    loci = "data.frame",
    strongest_risk_alleles = "data.frame",
    author_reported_genes = "data.frame",
    ensembl_gene_ids = "data.frame",
    entrez_gene_ids = "data.frame"
  )
)

#' Locus constructor
#'
#' Creates a Locus object from a list of raw records. To understand the
#' significance of each column of the data.frames, visit "Locus" in the
#' GWAS Catalog documentation
#' (https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference).
#'
#' @param data A list of locus records
#' @param mode The mode of the object: "Fat" or "Thin"
#'
#' @return A Locus object
#'
#' @examples
#' \dontrun{
#' ch <- Locus(
#'   list(list(association_id = "169106883", locus_id = "169106882")),
#'   mode = "Fat"
#' )
#' }
#'
#' @export
#'
Locus <- function(data = NULL, mode = "Fat") {
  if (!mode %in% c("Fat", "Thin")) {
    stop("Mode must be Fat or Thin")
  }

  if (is.null(data)) {
    normalized_data <- list()
  } else {
    normalized_data <- list()
    for (item in data) {
      if (!is.list(item) ||
          (length(item) > 0 && is.null(names(item)))) {
        next
      }
      rec <- item
      aliases <- list(
        list("association_id", c("association_id", "associationId")),
        list("locus_id", c("locus_id", "locusId")),
        list("haplotype_snp_count", c("haplotype_snp_count", "haplotypeSnpCount")),
        list("strongest_risk_alleles", c("strongest_risk_alleles", "strongestRiskAlleles")),
        list("author_reported_genes", c("author_reported_genes", "authorReportedGenes", "genes", "gene")),
        list("ensembl_gene_ids", c("ensembl_gene_ids", "ensemblGeneIds")),
        list("entrez_gene_ids", c("entrez_gene_ids", "entrezGeneIds"))
      )
      for (alias_entry in aliases) {
        key <- alias_entry[[1]]
        source_keys <- alias_entry[[2]]
        value <- NULL
        for (source_key in source_keys) {
          if (!is.null(rec[[source_key]])) {
            value <- rec[[source_key]]
            break
          }
        }
        if (!is.null(value)) {
          rec[[key]] <- value
        }
      }
      if (!is.null(rec$association_id)) {
        rec$association_id <- as.character(rec$association_id)
      }
      if (!is.null(rec$locus_id)) {
        rec$locus_id <- as.character(rec$locus_id)
      }
      normalized_data <- c(normalized_data, list(rec))
    }
  }

  LOCUS_COLUMNS <- c(
    "association_id",
    "locus_id",
    "haplotype_snp_count",
    "description"
  )

  STRONGEST_RISK_ALLELES_COLUMNS <- c(
    "association_id",
    "locus_id",
    "risk_allele_name",
    "risk_frequency",
    "genome_wide",
    "limited_list",
    "source",
    "mapping_method",
    "is_closest_gene"
  )

  AUTHOR_REPORTED_GENES_COLUMNS <- c(
    "association_id",
    "locus_id",
    "gene_name",
    "gene_description",
    "location",
    "cytogenic_region",
    "biotype"
  )

  ENSEMBL_GENE_IDS_COLUMNS <- c("association_id", "locus_id", "gene_name", "ensembl_gene_ids")
  ENTREZ_GENE_IDS_COLUMNS <- c("association_id", "locus_id", "gene_name", "entrez_gene_ids")

  NESTED_FIELDS <- c(
    "strongest_risk_alleles",
    "author_reported_genes",
    "ensembl_gene_ids",
    "entrez_gene_ids"
  )

  obj <- new("Locus")
  obj@mode <- mode

  obj@loci <- data.frame(matrix(ncol = length(LOCUS_COLUMNS), nrow = 0))
  colnames(obj@loci) <- LOCUS_COLUMNS

  obj@strongest_risk_alleles <- data.frame(matrix(ncol = length(STRONGEST_RISK_ALLELES_COLUMNS), nrow = 0))
  colnames(obj@strongest_risk_alleles) <- STRONGEST_RISK_ALLELES_COLUMNS

  obj@author_reported_genes <- data.frame(matrix(ncol = length(AUTHOR_REPORTED_GENES_COLUMNS), nrow = 0))
  colnames(obj@author_reported_genes) <- AUTHOR_REPORTED_GENES_COLUMNS

  obj@ensembl_gene_ids <- data.frame(matrix(ncol = length(ENSEMBL_GENE_IDS_COLUMNS), nrow = 0))
  colnames(obj@ensembl_gene_ids) <- ENSEMBL_GENE_IDS_COLUMNS

  obj@entrez_gene_ids <- data.frame(matrix(ncol = length(ENTREZ_GENE_IDS_COLUMNS), nrow = 0))
  colnames(obj@entrez_gene_ids) <- ENTREZ_GENE_IDS_COLUMNS

  if (mode == "Thin") {
    obj@raw_data <- normalized_data
    return(obj)
  }

  standardized_data <- list()
  for (item in normalized_data) {
    rec <- item
    aliases <- list(
      list("association_id", c("association_id", "associationId")),
      list("locus_id", c("locus_id", "locusId")),
      list("haplotype_snp_count", c("haplotype_snp_count", "haplotypeSnpCount")),
      list("strongest_risk_alleles", c("strongest_risk_alleles", "strongestRiskAlleles")),
      list("author_reported_genes", c("author_reported_genes", "authorReportedGenes", "genes", "gene")),
      list("ensembl_gene_ids", c("ensembl_gene_ids", "ensemblGeneIds")),
      list("entrez_gene_ids", c("entrez_gene_ids", "entrezGeneIds"))
    )
    for (alias_entry in aliases) {
      key <- alias_entry[[1]]
      source_keys <- alias_entry[[2]]
      value <- NULL
      for (source_key in source_keys) {
        if (!is.null(rec[[source_key]])) {
          value <- rec[[source_key]]
          break
        }
      }
      if (!is.null(value)) {
        rec[[key]] <- value
      }
    }

    standardized <- list()
    for (col in LOCUS_COLUMNS) {
      standardized[[col]] <- rec[[col]]
    }

    for (field in NESTED_FIELDS) {
      standardized[[field]] <- rec[[field]]
    }

    standardized_data <- c(standardized_data, list(standardized))
  }

  obj@raw_data <- standardized_data

  if (length(standardized_data) == 0) {
    return(obj)
  }

  locus_rows <- list()
  strongest_risk_alleles_rows <- list()
  author_reported_genes_rows <- list()
  ensembl_gene_ids_rows <- list()
  entrez_gene_ids_rows <- list()

  for (record in standardized_data) {
    association_id <- record$association_id
    locus_id <- record$locus_id

    locus_rows <- c(locus_rows, list(record[LOCUS_COLUMNS]))

    raw_strongest_risk_alleles <- record$strongest_risk_alleles
    if (is.null(raw_strongest_risk_alleles)) {
      iterable_strongest_risk_alleles <- list()
    } else if (is.list(raw_strongest_risk_alleles) && is.null(names(raw_strongest_risk_alleles))) {
      iterable_strongest_risk_alleles <- raw_strongest_risk_alleles
    } else {
      iterable_strongest_risk_alleles <- list(raw_strongest_risk_alleles)
    }
    for (item in iterable_strongest_risk_alleles) {
      if (is.list(item) && !is.null(names(item))) {
        risk_allele_name <- item$risk_allele_name
        if (is.null(risk_allele_name)) {
          risk_allele_name <- item$riskAlleleName
        }
        risk_frequency <- item$risk_frequency
        if (is.null(risk_frequency)) {
          risk_frequency <- item$riskFrequency
        }
        genome_wide <- item$genome_wide
        if (is.null(genome_wide)) {
          genome_wide <- item$genomeWide
        }
        limited_list <- item$limited_list
        if (is.null(limited_list)) {
          limited_list <- item$limitedList
        }
        mapping_method <- item$mapping_method
        if (is.null(mapping_method)) {
          mapping_method <- item$mappingMethod
        }
        is_closest_gene <- item$is_closest_gene
        if (is.null(is_closest_gene)) {
          is_closest_gene <- item$isClosestGene
        }
        strongest_risk_alleles_rows <- c(strongest_risk_alleles_rows, list(list(
          association_id = association_id,
          locus_id = locus_id,
          risk_allele_name = risk_allele_name,
          risk_frequency = risk_frequency,
          genome_wide = genome_wide,
          limited_list = limited_list,
          source = item$source,
          mapping_method = mapping_method,
          is_closest_gene = is_closest_gene
        )))
      } else {
        strongest_risk_alleles_rows <- c(strongest_risk_alleles_rows, list(list(
          association_id = association_id,
          locus_id = locus_id,
          risk_allele_name = item,
          risk_frequency = NA,
          genome_wide = NA,
          limited_list = NA,
          source = NA,
          mapping_method = NA,
          is_closest_gene = NA
        )))
      }
    }

    raw_author_reported_genes <- record$author_reported_genes
    if (is.null(raw_author_reported_genes)) {
      iterable_author_reported_genes <- list()
    } else if (is.list(raw_author_reported_genes) && is.null(names(raw_author_reported_genes))) {
      iterable_author_reported_genes <- raw_author_reported_genes
    } else {
      iterable_author_reported_genes <- list(raw_author_reported_genes)
    }
    for (item in iterable_author_reported_genes) {
      if (is.list(item) && !is.null(names(item))) {
        gene_name <- item$gene_name
        if (is.null(gene_name)) {
          gene_name <- item$geneName
        }
        gene_description <- item$gene_description
        if (is.null(gene_description)) {
          gene_description <- item$geneDescription
        }
        cytogenic_region <- item$cytogenic_region
        if (is.null(cytogenic_region)) {
          cytogenic_region <- item$cytogenetic_region
        }
        if (is.null(cytogenic_region)) {
          cytogenic_region <- item$cytogenicRegion
        }
        if (is.null(cytogenic_region)) {
          cytogenic_region <- item$cytogeneticRegion
        }
        author_reported_genes_rows <- c(author_reported_genes_rows, list(list(
          association_id = association_id,
          locus_id = locus_id,
          gene_name = gene_name,
          gene_description = gene_description,
          location = item$location,
          cytogenic_region = cytogenic_region,
          biotype = item$biotype
        )))

        raw_ensembl_gene_ids <- item$ensembl_gene_ids
        if (is.null(raw_ensembl_gene_ids)) {
          raw_ensembl_gene_ids <- item$ensemblGeneIds
        }
        if (is.null(raw_ensembl_gene_ids)) {
          iterable_ensembl_gene_ids <- list()
        } else if (is.list(raw_ensembl_gene_ids) && is.null(names(raw_ensembl_gene_ids))) {
          iterable_ensembl_gene_ids <- raw_ensembl_gene_ids
        } else {
          iterable_ensembl_gene_ids <- list(raw_ensembl_gene_ids)
        }
        for (ensembl_gene_id in iterable_ensembl_gene_ids) {
          ensembl_gene_ids_rows <- c(ensembl_gene_ids_rows, list(list(
            association_id = association_id,
            locus_id = locus_id,
            gene_name = gene_name,
            ensembl_gene_ids = ensembl_gene_id
          )))
        }

        raw_entrez_gene_ids <- item$entrez_gene_ids
        if (is.null(raw_entrez_gene_ids)) {
          raw_entrez_gene_ids <- item$entrezGeneIds
        }
        if (is.null(raw_entrez_gene_ids)) {
          iterable_entrez_gene_ids <- list()
        } else if (is.list(raw_entrez_gene_ids) && is.null(names(raw_entrez_gene_ids))) {
          iterable_entrez_gene_ids <- raw_entrez_gene_ids
        } else {
          iterable_entrez_gene_ids <- list(raw_entrez_gene_ids)
        }
        for (entrez_gene_id in iterable_entrez_gene_ids) {
          entrez_gene_ids_rows <- c(entrez_gene_ids_rows, list(list(
            association_id = association_id,
            locus_id = locus_id,
            gene_name = gene_name,
            entrez_gene_ids = entrez_gene_id
          )))
        }
      } else {
        author_reported_genes_rows <- c(author_reported_genes_rows, list(list(
          association_id = association_id,
          locus_id = locus_id,
          gene_name = item,
          gene_description = NA,
          location = NA,
          cytogenic_region = NA,
          biotype = NA
        )))
      }
    }

    raw_ensembl_gene_ids <- record$ensembl_gene_ids
    if (is.null(raw_ensembl_gene_ids)) {
      iterable_ensembl_gene_ids <- list()
    } else if (is.list(raw_ensembl_gene_ids) && is.null(names(raw_ensembl_gene_ids))) {
      iterable_ensembl_gene_ids <- raw_ensembl_gene_ids
    } else {
      iterable_ensembl_gene_ids <- list(raw_ensembl_gene_ids)
    }
    for (item in iterable_ensembl_gene_ids) {
      if (is.list(item) && !is.null(names(item))) {
        gene_name <- item$gene_name
        if (is.null(gene_name)) {
          gene_name <- item$geneName
        }
        ensembl_gene_ids <- item$ensembl_gene_ids
        if (is.null(ensembl_gene_ids)) {
          ensembl_gene_ids <- item$ensemblGeneIds
        }
        ensembl_gene_ids_rows <- c(ensembl_gene_ids_rows, list(list(
          association_id = association_id,
          locus_id = locus_id,
          gene_name = gene_name,
          ensembl_gene_ids = ensembl_gene_ids
        )))
      } else {
        ensembl_gene_ids_rows <- c(ensembl_gene_ids_rows, list(list(
          association_id = association_id,
          locus_id = locus_id,
          gene_name = NA,
          ensembl_gene_ids = item
        )))
      }
    }

    raw_entrez_gene_ids <- record$entrez_gene_ids
    if (is.null(raw_entrez_gene_ids)) {
      iterable_entrez_gene_ids <- list()
    } else if (is.list(raw_entrez_gene_ids) && is.null(names(raw_entrez_gene_ids))) {
      iterable_entrez_gene_ids <- raw_entrez_gene_ids
    } else {
      iterable_entrez_gene_ids <- list(raw_entrez_gene_ids)
    }
    for (item in iterable_entrez_gene_ids) {
      if (is.list(item) && !is.null(names(item))) {
        gene_name <- item$gene_name
        if (is.null(gene_name)) {
          gene_name <- item$geneName
        }
        entrez_gene_ids <- item$entrez_gene_ids
        if (is.null(entrez_gene_ids)) {
          entrez_gene_ids <- item$entrezGeneIds
        }
        entrez_gene_ids_rows <- c(entrez_gene_ids_rows, list(list(
          association_id = association_id,
          locus_id = locus_id,
          gene_name = gene_name,
          entrez_gene_ids = entrez_gene_ids
        )))
      } else {
        entrez_gene_ids_rows <- c(entrez_gene_ids_rows, list(list(
          association_id = association_id,
          locus_id = locus_id,
          gene_name = NA,
          entrez_gene_ids = item
        )))
      }
    }
  }

  if (length(locus_rows) > 0) {
    built <- data.frame(row.names = seq_along(locus_rows))
    for (col in LOCUS_COLUMNS) {
      col_values <- vector("list", length(locus_rows))
      for (r_idx in seq_along(locus_rows)) {
        cell <- locus_rows[[r_idx]][[col]]
        if (is.null(cell)) {
          cell <- NA
        }
        col_values[[r_idx]] <- cell
      }
      atomic_col <- TRUE
      for (cell in col_values) {
        if (!is.atomic(cell) || length(cell) != 1) {
          atomic_col <- FALSE
          break
        }
      }
      if (atomic_col) {
        built[[col]] <- unlist(col_values)
      } else {
        built[[col]] <- col_values
      }
    }
    obj@loci <- built
  }

  if (length(strongest_risk_alleles_rows) > 0) {
    built <- data.frame(row.names = seq_along(strongest_risk_alleles_rows))
    for (col in STRONGEST_RISK_ALLELES_COLUMNS) {
      col_values <- vector("list", length(strongest_risk_alleles_rows))
      for (r_idx in seq_along(strongest_risk_alleles_rows)) {
        cell <- strongest_risk_alleles_rows[[r_idx]][[col]]
        if (is.null(cell)) {
          cell <- NA
        }
        col_values[[r_idx]] <- cell
      }
      atomic_col <- TRUE
      for (cell in col_values) {
        if (!is.atomic(cell) || length(cell) != 1) {
          atomic_col <- FALSE
          break
        }
      }
      if (atomic_col) {
        built[[col]] <- unlist(col_values)
      } else {
        built[[col]] <- col_values
      }
    }
    obj@strongest_risk_alleles <- built
  }

  if (length(author_reported_genes_rows) > 0) {
    built <- data.frame(row.names = seq_along(author_reported_genes_rows))
    for (col in AUTHOR_REPORTED_GENES_COLUMNS) {
      col_values <- vector("list", length(author_reported_genes_rows))
      for (r_idx in seq_along(author_reported_genes_rows)) {
        cell <- author_reported_genes_rows[[r_idx]][[col]]
        if (is.null(cell)) {
          cell <- NA
        }
        col_values[[r_idx]] <- cell
      }
      atomic_col <- TRUE
      for (cell in col_values) {
        if (!is.atomic(cell) || length(cell) != 1) {
          atomic_col <- FALSE
          break
        }
      }
      if (atomic_col) {
        built[[col]] <- unlist(col_values)
      } else {
        built[[col]] <- col_values
      }
    }
    obj@author_reported_genes <- built
  }

  if (length(ensembl_gene_ids_rows) > 0) {
    built <- data.frame(row.names = seq_along(ensembl_gene_ids_rows))
    for (col in ENSEMBL_GENE_IDS_COLUMNS) {
      col_values <- vector("list", length(ensembl_gene_ids_rows))
      for (r_idx in seq_along(ensembl_gene_ids_rows)) {
        cell <- ensembl_gene_ids_rows[[r_idx]][[col]]
        if (is.null(cell)) {
          cell <- NA
        }
        col_values[[r_idx]] <- cell
      }
      atomic_col <- TRUE
      for (cell in col_values) {
        if (!is.atomic(cell) || length(cell) != 1) {
          atomic_col <- FALSE
          break
        }
      }
      if (atomic_col) {
        built[[col]] <- unlist(col_values)
      } else {
        built[[col]] <- col_values
      }
    }
    obj@ensembl_gene_ids <- built
  }

  if (length(entrez_gene_ids_rows) > 0) {
    built <- data.frame(row.names = seq_along(entrez_gene_ids_rows))
    for (col in ENTREZ_GENE_IDS_COLUMNS) {
      col_values <- vector("list", length(entrez_gene_ids_rows))
      for (r_idx in seq_along(entrez_gene_ids_rows)) {
        cell <- entrez_gene_ids_rows[[r_idx]][[col]]
        if (is.null(cell)) {
          cell <- NA
        }
        col_values[[r_idx]] <- cell
      }
      atomic_col <- TRUE
      for (cell in col_values) {
        if (!is.atomic(cell) || length(cell) != 1) {
          atomic_col <- FALSE
          break
        }
      }
      if (atomic_col) {
        built[[col]] <- unlist(col_values)
      } else {
        built[[col]] <- col_values
      }
    }
    obj@entrez_gene_ids <- built
  }

  return(obj)
}

#' Show a Locus object
#'
#' Print a concise overview of a Locus object: the current mode and the
#' size of each data component.
#'
#' @param object A Locus object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "show",
  signature = "Locus",
  function(object) {
    if (object@mode == "Fat") {
      cat(
        "Locus is running in fat mode. It has 5 DataFrames with hierarchical dependencies.\n",
        "loci: ", nrow(object@loci), " rows\n",
        "|\n -strongest_risk_alleles: ", nrow(object@strongest_risk_alleles), " rows\n",
        "|\n -author_reported_genes: ", nrow(object@author_reported_genes), " rows\n",
        "|\n -ensembl_gene_ids: ", nrow(object@ensembl_gene_ids), " rows\n",
        "|\n -entrez_gene_ids: ", nrow(object@entrez_gene_ids), " rows\n",
        sep = ""
      )
    } else if (object@mode == "Thin") {
      cat(
        "Locus is running in thin mode. It has 1 list that contains the raw data.\n",
        "raw_data: a list of size ", length(object@raw_data), ".\n",
        sep = ""
      )
    }
  }
)

#' Length of a Locus object
#'
#' Get the number of locus records stored in the object.
#'
#' @param x A Locus object
#'
#' @return An integer, the number of records in the raw_data slot
#'
#' @export
setMethod(
  "length",
  signature = "Locus",
  function(x) {
    return(length(x@raw_data))
  }
)

#' Component names of a Locus object
#'
#' Get the names of the components (slots) of a Locus object.
#'
#' @param x A Locus object
#'
#' @return A character vector of component names
#'
#' @export
setMethod(
  "names",
  signature = "Locus",
  function(x) {
    return(c(
      "raw_data",
      "mode",
      "loci",
      "strongest_risk_alleles",
      "author_reported_genes",
      "ensembl_gene_ids",
      "entrez_gene_ids"
    ))
  }
)

#' Subset a Locus object
#'
#' Extract locus records by position or by the composite identifier
#' "association_id|locus_id" and return a new Locus object of the same
#' mode containing only the selected records.
#'
#' @param x A Locus object
#' @param i A numeric vector of positions (1-based) or a character vector of
#'   composite identifiers in the form "association_id|locus_id". A list
#'   is unlisted before use. Leaving i missing returns x unchanged.
#'
#' @return A Locus object containing the selected records
#'
#' @examples
#' \dontrun{
#' a <- get_loci(association_id = "169106883", locus_id = "169106882")
#' a[1]
#' a["169106883|169106882"]
#' }
#'
#' @export
setMethod(
  "[",
  signature = "Locus",
  function(x, i) {
    if (missing(i)) {
      return(x)
    }

    if (is.character(i)) {
      arr <- i
    } else if (is.numeric(i)) {
      arr <- i
    } else if (is.list(i)) {
      arr <- unlist(i)
    } else {
      stop("Invalid argument type: ", class(i))
    }

    raw_data <- x@raw_data
    raw_data_dict <- list()
    for (j in seq_along(raw_data)) {
      if (!is.null(raw_data[[j]]$association_id) &&
          !is.null(raw_data[[j]]$locus_id)) {
        key <- paste0(raw_data[[j]]$association_id, "|", raw_data[[j]]$locus_id)
        raw_data_dict[[key]] <- raw_data[[j]]
      }
    }

    sub_set <- list()
    for (i_val in arr) {
      if (is.character(i_val)) {
        if (!i_val %in% names(raw_data_dict)) {
          stop("Unknown association_id|locus_id: ", i_val)
        }
        sub_set <- c(sub_set, list(raw_data_dict[[i_val]]))
      } else if (is.numeric(i_val)) {
        position <- as.integer(i_val)
        if (is.na(position) || position < 1 ||
            position > length(raw_data)) {
          stop("subscript out of bounds: ", i_val)
        }
        sub_set <- c(sub_set, list(raw_data[[position]]))
      } else {
        stop("Invalid item type: ", class(i_val))
      }
    }

    return(Locus(sub_set, x@mode))
  }
)

#' Set operations on Locus objects
#'
#' Combine or compare two Locus objects of the same mode as sets of locus
#' records, identified and deduplicated by the composite key
#' "association_id|locus_id".
#'
#' @param e1 A Locus object
#' @param e2 A Locus object
#'
#' @return A Locus object of the same mode
#'
#' @examples
#' \dontrun{
#' a <- get_loci(association_id = "169106883", locus_id = "169106882")
#' b <- get_loci(association_id = "218958492", locus_id = "218958491")
#' a + b
#' a - b
#' a & b
#' a | b
#' a ^ b
#' }
#'
#' @rdname Locus_set_operations
#' @export
setMethod(
  "+",
  signature = c("Locus", "Locus"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      return(Locus(c(e1@raw_data, e2@raw_data), e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Locus_set_operations
#' @export
setMethod(
  "-",
  signature = c("Locus", "Locus"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$association_id) && !is.null(i$locus_id)) {
          key <- paste0(i$association_id, "|", i$locus_id)
          self_key_set <- unique(c(self_key_set, key))
          self_dict[[key]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$association_id) && !is.null(j$locus_id)) {
          other_key_set <- unique(c(other_key_set, paste0(j$association_id, "|", j$locus_id)))
        }
      }

      sub_key <- self_key_set[!self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Locus(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Locus_set_operations
#' @export
setMethod(
  "&",
  signature = c("Locus", "Locus"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$association_id) && !is.null(i$locus_id)) {
          key <- paste0(i$association_id, "|", i$locus_id)
          self_key_set <- unique(c(self_key_set, key))
          self_dict[[key]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$association_id) && !is.null(j$locus_id)) {
          other_key_set <- unique(c(other_key_set, paste0(j$association_id, "|", j$locus_id)))
        }
      }

      sub_key <- self_key_set[self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Locus(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Locus_set_operations
#' @export
setMethod(
  "|",
  signature = c("Locus", "Locus"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      and_dict <- list()

      for (i in e1@raw_data) {
        if (!is.null(i$association_id) && !is.null(i$locus_id)) {
          and_dict[[paste0(i$association_id, "|", i$locus_id)]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$association_id) && !is.null(j$locus_id)) {
          and_dict[[paste0(j$association_id, "|", j$locus_id)]] <- j
        }
      }

      data <- unname(and_dict)
      return(Locus(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Locus_set_operations
#' @export
setMethod(
  "^",
  signature = c("Locus", "Locus"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      and_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$association_id) && !is.null(i$locus_id)) {
          key <- paste0(i$association_id, "|", i$locus_id)
          self_key_set <- unique(c(self_key_set, key))
          and_dict[[key]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$association_id) && !is.null(j$locus_id)) {
          key <- paste0(j$association_id, "|", j$locus_id)
          other_key_set <- unique(c(other_key_set, key))
          and_dict[[key]] <- j
        }
      }

      all_keys <- unique(c(self_key_set, other_key_set))
      common_keys <- self_key_set[self_key_set %in% other_key_set]
      sub_key <- all_keys[!all_keys %in% common_keys]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(and_dict[[k]]))
      }

      return(Locus(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' Compare Locus objects for equality
#'
#' Test whether two Locus objects are equal. Two objects are equal when
#' they run in the same mode and store identical records.
#'
#' @param e1 A Locus object
#' @param e2 A Locus object
#'
#' @return A logical value
#'
#' @export
setMethod(
  "==",
  signature = c("Locus", "Locus"),
  function(e1, e2) {
    if (is.null(e1) || is.null(e2)) {
      return(is.null(e1) && is.null(e2))
    }
    return(identical(e1@raw_data, e2@raw_data) && e1@mode == e2@mode)
  }
)

#' Summarize a Locus object
#'
#' Print a summary of a Locus object: the mode, the number of raw records,
#' the number of loci and, in Fat mode, the row counts of all 5
#' data.frames.
#'
#' @param object A Locus object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "summary",
  signature = "Locus",
  function(object) {
    cat("Locus Object Summary\n")
    cat("====================\n")
    cat("Mode:", object@mode, "\n")
    cat("Number of raw records:", length(object@raw_data), "\n")
    cat("Number of loci:", nrow(object@loci), "\n")
    if (object@mode == "Fat") {
      cat("\nDataFrames:\n")
      cat("  loci:", nrow(object@loci), "rows\n")
      cat("  strongest_risk_alleles:", nrow(object@strongest_risk_alleles), "rows\n")
      cat("  author_reported_genes:", nrow(object@author_reported_genes), "rows\n")
      cat("  ensembl_gene_ids:", nrow(object@ensembl_gene_ids), "rows\n")
      cat("  entrez_gene_ids:", nrow(object@entrez_gene_ids), "rows\n")
    }
  }
)
