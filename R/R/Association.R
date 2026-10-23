#' Association class
#'
#' An S4 class that stores data of type Association fetched from the GWAS API.
#' In Fat mode it holds 8 data.frames with hierarchical dependencies; in Thin
#' mode it only stores the raw list data.
#'
#' @slot raw_data list. Converted from obtained JSON data
#' @slot mode character. Fat or Thin
#' @slot associations data.frame
#' @slot efo_traits data.frame
#' @slot snp_allele data.frame
#' @slot bg_efo_traits data.frame
#' @slot reported_traits data.frame
#' @slot locations data.frame
#' @slot mapped_genes data.frame
#' @slot snp_effect_alleles data.frame
#'
#' @importFrom methods new
#'
#' @export
#'
setClass(
  "Association",
  slots = c(
    raw_data = "list",
    mode = "character",
    associations = "data.frame",
    efo_traits = "data.frame",
    snp_allele = "data.frame",
    bg_efo_traits = "data.frame",
    reported_traits = "data.frame",
    locations = "data.frame",
    mapped_genes = "data.frame",
    snp_effect_alleles = "data.frame"
  )
)

#' Association constructor
#'
#' Creates an Association object from a list of raw records. To understand the
#' significance of each column of the data.frames, visit "Association" in the
#' GWAS Catalog documentation
#' (https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference).
#'
#' @param data A list of association records
#' @param mode The mode of the object: "Fat" or "Thin"
#'
#' @return An Association object
#'
#' @examples
#' \dontrun{
#' ch <- Association(
#'   list(list(association_id = "169106883", pvalue_mantissa = 5)),
#'   mode = "Fat"
#' )
#' }
#'
#' @export
#'
Association <- function(data = NULL, mode = "Fat") {
  if (!mode %in% c("Fat", "Thin")) {
    stop("Mode must be Fat or Thin")
  }

  if (is.null(data)) {
    normalized_data <- list()
  } else {
    normalized_data <- list()
    for (item in data) {
      # JSON objects parse to named lists; arrays parse to unnamed lists
      if (!is.list(item) ||
          (length(item) > 0 && is.null(names(item)))) {
        next
      }
      rec <- item
      if (!is.null(rec$association_id)) {
        rec$association_id <- as.character(rec$association_id)
      }
      normalized_data <- c(normalized_data, list(rec))
    }
  }

  ASSOCIATION_COLUMNS <- c(
    "association_id",
    "risk_frequency",
    "pvalue_description",
    "pvalue_mantissa",
    "pvalue_exponent",
    "multi_snp_haplotype",
    "snp_interaction",
    "snp_type",
    "range",
    "description",
    "or_per_copy_num",
    "or_value",
    "beta_num",
    "beta_unit",
    "beta_direction",
    "beta",
    "last_mapping_date",
    "last_update_date",
    "p_value",
    "pubmed_id",
    "first_author",
    "ci_lower",
    "ci_upper"
  )

  EFO_TRAITS_COLUMNS <- c("association_id", "efo_id", "efo_trait")
  SNP_ALLELE_COLUMNS <- c("association_id", "rs_id", "effect_allele")
  BG_EFO_TRAITS_COLUMNS <- c("association_id", "efo_id", "efo_trait")
  REPORTED_TRAITS_COLUMNS <- c("association_id", "reported_trait")
  LOCATIONS_COLUMNS <- c("association_id", "locations")
  MAPPED_GENES_COLUMNS <- c("association_id", "mapped_genes")
  SNP_EFFECT_ALLELES_COLUMNS <- c("association_id", "snp_effect_allele")

  NESTED_FIELDS <- c(
    "efo_traits",
    "snp_allele",
    "bg_efo_traits",
    "reported_trait",
    "locations",
    "mapped_genes",
    "snp_effect_allele"
  )

  obj <- new("Association")
  obj@mode <- mode

  obj@associations <- data.frame(matrix(ncol = length(ASSOCIATION_COLUMNS), nrow = 0))
  colnames(obj@associations) <- ASSOCIATION_COLUMNS

  obj@efo_traits <- data.frame(matrix(ncol = length(EFO_TRAITS_COLUMNS), nrow = 0))
  colnames(obj@efo_traits) <- EFO_TRAITS_COLUMNS

  obj@snp_allele <- data.frame(matrix(ncol = length(SNP_ALLELE_COLUMNS), nrow = 0))
  colnames(obj@snp_allele) <- SNP_ALLELE_COLUMNS

  obj@bg_efo_traits <- data.frame(matrix(ncol = length(BG_EFO_TRAITS_COLUMNS), nrow = 0))
  colnames(obj@bg_efo_traits) <- BG_EFO_TRAITS_COLUMNS

  obj@reported_traits <- data.frame(matrix(ncol = length(REPORTED_TRAITS_COLUMNS), nrow = 0))
  colnames(obj@reported_traits) <- REPORTED_TRAITS_COLUMNS

  obj@locations <- data.frame(matrix(ncol = length(LOCATIONS_COLUMNS), nrow = 0))
  colnames(obj@locations) <- LOCATIONS_COLUMNS

  obj@mapped_genes <- data.frame(matrix(ncol = length(MAPPED_GENES_COLUMNS), nrow = 0))
  colnames(obj@mapped_genes) <- MAPPED_GENES_COLUMNS

  obj@snp_effect_alleles <- data.frame(matrix(ncol = length(SNP_EFFECT_ALLELES_COLUMNS), nrow = 0))
  colnames(obj@snp_effect_alleles) <- SNP_EFFECT_ALLELES_COLUMNS

  if (mode == "Thin") {
    obj@raw_data <- normalized_data
    return(obj)
  }

  # Fat mode: standardize fixed columns, keep nested fields, fill missing with NA
  standardized_data <- list()
  for (item in normalized_data) {
    standardized <- list()
    for (col in ASSOCIATION_COLUMNS) {
      value <- item[[col]]
      if (col == "range") {
        if (is.null(value)) {
          value <- NA
        } else if (is.list(value)) {
          if (length(value) == 0) {
            value <- NA
          } else if (length(value) == 1) {
            value <- if (is.null(value[[1]])) NA else as.character(value[[1]])
          } else {
            pieces <- character(length(value))
            for (v_idx in seq_along(value)) {
              pieces[v_idx] <- if (is.null(value[[v_idx]])) {
                ""
              } else {
                as.character(value[[v_idx]])
              }
            }
            value <- paste(pieces, collapse = "; ")
          }
        } else {
          value <- as.character(value)
        }
      }
      standardized[[col]] <- value
    }

    for (field in NESTED_FIELDS) {
      standardized[[field]] <- item[[field]]
    }

    standardized_data <- c(standardized_data, list(standardized))
  }

  obj@raw_data <- standardized_data

  if (length(standardized_data) == 0) {
    return(obj)
  }

  association_rows <- list()
  efo_traits_rows <- list()
  snp_allele_rows <- list()
  bg_efo_traits_rows <- list()
  reported_traits_rows <- list()
  locations_rows <- list()
  mapped_genes_rows <- list()
  snp_effect_alleles_rows <- list()

  for (record in standardized_data) {
    association_id <- record$association_id

    # Main table: fixed columns, missing filled with NA
    association_rows <- c(
      association_rows,
      list(record[ASSOCIATION_COLUMNS])
    )

    # efo_traits
    raw_efo_traits <- record$efo_traits
    if (is.null(raw_efo_traits)) {
      iterable_efo_traits <- list()
    } else if (is.list(raw_efo_traits) && is.null(names(raw_efo_traits))) {
      iterable_efo_traits <- raw_efo_traits
    } else {
      iterable_efo_traits <- list(raw_efo_traits)
    }
    for (item in iterable_efo_traits) {
      if (is.list(item) && !is.null(names(item))) {
        efo_traits_rows <- c(efo_traits_rows, list(list(
          association_id = association_id,
          efo_id = item$efo_id,
          efo_trait = item$efo_trait
        )))
      } else {
        efo_traits_rows <- c(efo_traits_rows, list(list(
          association_id = association_id,
          efo_id = NA,
          efo_trait = item
        )))
      }
    }

    # snp_allele
    raw_snp_allele <- record$snp_allele
    if (is.null(raw_snp_allele)) {
      iterable_snp_allele <- list()
    } else if (is.list(raw_snp_allele) && is.null(names(raw_snp_allele))) {
      iterable_snp_allele <- raw_snp_allele
    } else {
      iterable_snp_allele <- list(raw_snp_allele)
    }
    for (item in iterable_snp_allele) {
      if (is.list(item) && !is.null(names(item))) {
        snp_allele_rows <- c(snp_allele_rows, list(list(
          association_id = association_id,
          rs_id = item$rs_id,
          effect_allele = item$effect_allele
        )))
      } else {
        snp_allele_rows <- c(snp_allele_rows, list(list(
          association_id = association_id,
          rs_id = item,
          effect_allele = NA
        )))
      }
    }

    # bg_efo_traits
    raw_bg_efo_traits <- record$bg_efo_traits
    if (is.null(raw_bg_efo_traits)) {
      iterable_bg_efo_traits <- list()
    } else if (is.list(raw_bg_efo_traits) && is.null(names(raw_bg_efo_traits))) {
      iterable_bg_efo_traits <- raw_bg_efo_traits
    } else {
      iterable_bg_efo_traits <- list(raw_bg_efo_traits)
    }
    for (item in iterable_bg_efo_traits) {
      if (is.list(item) && !is.null(names(item))) {
        bg_efo_traits_rows <- c(bg_efo_traits_rows, list(list(
          association_id = association_id,
          efo_id = item$efo_id,
          efo_trait = item$efo_trait
        )))
      } else {
        bg_efo_traits_rows <- c(bg_efo_traits_rows, list(list(
          association_id = association_id,
          efo_id = NA,
          efo_trait = item
        )))
      }
    }

    # reported_traits (attribute stays plural, data field is reported_trait)
    raw_reported_trait <- record$reported_trait
    if (is.null(raw_reported_trait)) {
      iterable_reported_trait <- list()
    } else if (is.list(raw_reported_trait) && is.null(names(raw_reported_trait))) {
      iterable_reported_trait <- raw_reported_trait
    } else {
      iterable_reported_trait <- list(raw_reported_trait)
    }
    for (item in iterable_reported_trait) {
      reported_traits_rows <- c(reported_traits_rows, list(list(
        association_id = association_id,
        reported_trait = item
      )))
    }

    # locations
    raw_locations <- record$locations
    if (is.null(raw_locations)) {
      iterable_locations <- list()
    } else if (is.list(raw_locations) && is.null(names(raw_locations))) {
      iterable_locations <- raw_locations
    } else {
      iterable_locations <- list(raw_locations)
    }
    for (item in iterable_locations) {
      locations_rows <- c(locations_rows, list(list(
        association_id = association_id,
        locations = item
      )))
    }

    # mapped_genes
    raw_mapped_genes <- record$mapped_genes
    if (is.null(raw_mapped_genes)) {
      iterable_mapped_genes <- list()
    } else if (is.list(raw_mapped_genes) && is.null(names(raw_mapped_genes))) {
      iterable_mapped_genes <- raw_mapped_genes
    } else {
      iterable_mapped_genes <- list(raw_mapped_genes)
    }
    for (item in iterable_mapped_genes) {
      mapped_genes_rows <- c(mapped_genes_rows, list(list(
        association_id = association_id,
        mapped_genes = item
      )))
    }

    # snp_effect_alleles (attribute stays plural, data field is snp_effect_allele)
    raw_snp_effect_allele <- record$snp_effect_allele
    if (is.null(raw_snp_effect_allele)) {
      iterable_snp_effect_allele <- list()
    } else if (is.list(raw_snp_effect_allele) && is.null(names(raw_snp_effect_allele))) {
      iterable_snp_effect_allele <- raw_snp_effect_allele
    } else {
      iterable_snp_effect_allele <- list(raw_snp_effect_allele)
    }
    for (item in iterable_snp_effect_allele) {
      snp_effect_alleles_rows <- c(snp_effect_alleles_rows, list(list(
        association_id = association_id,
        snp_effect_allele = item
      )))
    }
  }

  # Build associations data.frame
  if (length(association_rows) > 0) {
    built <- data.frame(row.names = seq_along(association_rows))
    for (col in ASSOCIATION_COLUMNS) {
      col_values <- vector("list", length(association_rows))
      for (r_idx in seq_along(association_rows)) {
        cell <- association_rows[[r_idx]][[col]]
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
    obj@associations <- built
  }

  # Build efo_traits data.frame
  if (length(efo_traits_rows) > 0) {
    built <- data.frame(row.names = seq_along(efo_traits_rows))
    for (col in EFO_TRAITS_COLUMNS) {
      col_values <- vector("list", length(efo_traits_rows))
      for (r_idx in seq_along(efo_traits_rows)) {
        cell <- efo_traits_rows[[r_idx]][[col]]
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
    obj@efo_traits <- built
  }

  # Build snp_allele data.frame
  if (length(snp_allele_rows) > 0) {
    built <- data.frame(row.names = seq_along(snp_allele_rows))
    for (col in SNP_ALLELE_COLUMNS) {
      col_values <- vector("list", length(snp_allele_rows))
      for (r_idx in seq_along(snp_allele_rows)) {
        cell <- snp_allele_rows[[r_idx]][[col]]
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
    obj@snp_allele <- built
  }

  # Build bg_efo_traits data.frame
  if (length(bg_efo_traits_rows) > 0) {
    built <- data.frame(row.names = seq_along(bg_efo_traits_rows))
    for (col in BG_EFO_TRAITS_COLUMNS) {
      col_values <- vector("list", length(bg_efo_traits_rows))
      for (r_idx in seq_along(bg_efo_traits_rows)) {
        cell <- bg_efo_traits_rows[[r_idx]][[col]]
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
    obj@bg_efo_traits <- built
  }

  # Build reported_traits data.frame
  if (length(reported_traits_rows) > 0) {
    built <- data.frame(row.names = seq_along(reported_traits_rows))
    for (col in REPORTED_TRAITS_COLUMNS) {
      col_values <- vector("list", length(reported_traits_rows))
      for (r_idx in seq_along(reported_traits_rows)) {
        cell <- reported_traits_rows[[r_idx]][[col]]
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
    obj@reported_traits <- built
  }

  # Build locations data.frame
  if (length(locations_rows) > 0) {
    built <- data.frame(row.names = seq_along(locations_rows))
    for (col in LOCATIONS_COLUMNS) {
      col_values <- vector("list", length(locations_rows))
      for (r_idx in seq_along(locations_rows)) {
        cell <- locations_rows[[r_idx]][[col]]
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
    obj@locations <- built
  }

  # Build mapped_genes data.frame
  if (length(mapped_genes_rows) > 0) {
    built <- data.frame(row.names = seq_along(mapped_genes_rows))
    for (col in MAPPED_GENES_COLUMNS) {
      col_values <- vector("list", length(mapped_genes_rows))
      for (r_idx in seq_along(mapped_genes_rows)) {
        cell <- mapped_genes_rows[[r_idx]][[col]]
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
    obj@mapped_genes <- built
  }

  # Build snp_effect_alleles data.frame
  if (length(snp_effect_alleles_rows) > 0) {
    built <- data.frame(row.names = seq_along(snp_effect_alleles_rows))
    for (col in SNP_EFFECT_ALLELES_COLUMNS) {
      col_values <- vector("list", length(snp_effect_alleles_rows))
      for (r_idx in seq_along(snp_effect_alleles_rows)) {
        cell <- snp_effect_alleles_rows[[r_idx]][[col]]
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
    obj@snp_effect_alleles <- built
  }

  return(obj)
}

#' Show an Association object
#'
#' Print a concise overview of an Association object: the current mode and
#' the size of each data component. In Fat mode the row counts of the 8
#' data.frames are listed; in Thin mode the size of the raw record list is
#' shown.
#'
#' @param object An Association object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "show",
  signature = "Association",
  function(object) {
    if (object@mode == "Fat") {
      cat(
        "Association is running in fat mode. It has 8 DataFrames with hierarchical dependencies.\n",
        "associations: ", nrow(object@associations), " rows\n",
        "|\n -efo_traits: ", nrow(object@efo_traits), " rows\n",
        "|\n -snp_allele: ", nrow(object@snp_allele), " rows\n",
        "|\n -bg_efo_traits: ", nrow(object@bg_efo_traits), " rows\n",
        "|\n -reported_traits: ", nrow(object@reported_traits), " rows\n",
        "|\n -locations: ", nrow(object@locations), " rows\n",
        "|\n -mapped_genes: ", nrow(object@mapped_genes), " rows\n",
        "|\n -snp_effect_alleles: ", nrow(object@snp_effect_alleles), " rows\n",
        sep = ""
      )
    } else if (object@mode == "Thin") {
      cat(
        "Association is running in thin mode. It has 1 list that contains the raw data.\n",
        "raw_data: a list of size ", length(object@raw_data), ".\n",
        sep = ""
      )
    }
  }
)

#' Length of an Association object
#'
#' Get the number of association records stored in the object.
#'
#' @param x An Association object
#'
#' @return An integer, the number of records in the raw_data slot
#'
#' @export
setMethod(
  "length",
  signature = "Association",
  function(x) {
    return(length(x@raw_data))
  }
)

#' Component names of an Association object
#'
#' Get the names of the components (slots) of an Association object.
#'
#' @param x An Association object
#'
#' @return A character vector of component names
#'
#' @export
setMethod(
  "names",
  signature = "Association",
  function(x) {
    return(c(
      "raw_data",
      "mode",
      "associations",
      "efo_traits",
      "snp_allele",
      "bg_efo_traits",
      "reported_traits",
      "locations",
      "mapped_genes",
      "snp_effect_alleles"
    ))
  }
)

#' Subset an Association object
#'
#' Extract association records by position or by association identifier and
#' return a new Association object of the same mode containing only the
#' selected records.
#'
#' @param x An Association object
#' @param i A numeric vector of positions (1-based) or a character vector of
#'   association_ids. A list is unlisted before use, so a list of positions
#'   or identifiers is also accepted. Leaving i missing returns x unchanged.
#'
#' @return An Association object containing the selected records
#'
#' @examples
#' \dontrun{
#' a <- get_associations(accession_id = "GCST000854")
#' a[1]
#' a["169106883"]
#' }
#'
#' @export
setMethod(
  "[",
  signature = "Association",
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
      if (!is.null(raw_data[[j]]$association_id)) {
        raw_data_dict[[raw_data[[j]]$association_id]] <- raw_data[[j]]
      }
    }

    sub_set <- list()
    for (i_val in arr) {
      if (is.character(i_val)) {
        if (!i_val %in% names(raw_data_dict)) {
          stop("Unknown association_id: ", i_val)
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

    return(Association(sub_set, x@mode))
  }
)

#' Set operations on Association objects
#'
#' Combine or compare two Association objects of the same mode as sets of
#' association records, identified and deduplicated by association_id.
#'
#' `+` concatenates the records of both objects, in the order of e1 then
#' e2. `-` returns the records of e1 whose association_id does not appear
#' in e2 (difference). `&` returns the records of e1 whose association_id
#' also appears in e2 (intersection). `|` returns the records of both
#' objects deduplicated by association_id, keeping the record of e1 when
#' both objects contain the same id (union). `^` returns the records whose
#' association_id appears in exactly one of the two objects (symmetric
#' difference). All operators raise an error when the two objects do not
#' run in the same mode.
#'
#' @param e1 An Association object
#' @param e2 An Association object
#'
#' @return An Association object of the same mode
#'
#' @examples
#' \dontrun{
#' a <- get_associations(efo_trait = "celiac disease")
#' b <- get_associations(mapped_gene = "CCR6")
#' a + b
#' a - b
#' a & b
#' a | b
#' a ^ b
#' }
#'
#' @rdname set_operations
#' @export
setMethod(
  "+",
  signature = c("Association", "Association"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      return(Association(c(e1@raw_data, e2@raw_data), e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname set_operations
#' @export
setMethod(
  "-",
  signature = c("Association", "Association"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$association_id)) {
          self_key_set <- unique(c(self_key_set, i$association_id))
          self_dict[[i$association_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$association_id)) {
          other_key_set <- unique(c(other_key_set, j$association_id))
        }
      }

      sub_key <- self_key_set[!self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Association(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname set_operations
#' @export
setMethod(
  "&",
  signature = c("Association", "Association"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$association_id)) {
          self_key_set <- unique(c(self_key_set, i$association_id))
          self_dict[[i$association_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$association_id)) {
          other_key_set <- unique(c(other_key_set, j$association_id))
        }
      }

      sub_key <- self_key_set[self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Association(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname set_operations
#' @export
setMethod(
  "|",
  signature = c("Association", "Association"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      and_dict <- list()

      for (i in e1@raw_data) {
        if (!is.null(i$association_id)) {
          and_dict[[i$association_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$association_id)) {
          and_dict[[j$association_id]] <- j
        }
      }

      data <- unname(and_dict)
      return(Association(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname set_operations
#' @export
setMethod(
  "^",
  signature = c("Association", "Association"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      and_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$association_id)) {
          self_key_set <- unique(c(self_key_set, i$association_id))
          and_dict[[i$association_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$association_id)) {
          other_key_set <- unique(c(other_key_set, j$association_id))
          and_dict[[j$association_id]] <- j
        }
      }

      all_keys <- unique(c(self_key_set, other_key_set))
      common_keys <- self_key_set[self_key_set %in% other_key_set]
      sub_key <- all_keys[!all_keys %in% common_keys]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(and_dict[[k]]))
      }

      return(Association(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' Compare Association objects for equality
#'
#' Test whether two Association objects are equal. Two objects are equal
#' when they run in the same mode and store identical records.
#'
#' @param e1 An Association object
#' @param e2 An Association object
#'
#' @return A logical value
#'
#' @export
setMethod(
  "==",
  signature = c("Association", "Association"),
  function(e1, e2) {
    if (is.null(e1) || is.null(e2)) {
      return(is.null(e1) && is.null(e2))
    }
    return(identical(e1@raw_data, e2@raw_data) && e1@mode == e2@mode)
  }
)

#' Summarize an Association object
#'
#' Print a summary of an Association object: the mode, the number of raw
#' records, the number of associations and, in Fat mode, the row counts of
#' all 8 data.frames.
#'
#' @param object An Association object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "summary",
  signature = "Association",
  function(object) {
    cat("Association Object Summary\n")
    cat("==========================\n")
    cat("Mode:", object@mode, "\n")
    cat("Number of raw records:", length(object@raw_data), "\n")
    cat("Number of associations:", nrow(object@associations), "\n")
    if (object@mode == "Fat") {
      cat("\nDataFrames:\n")
      cat("  associations:", nrow(object@associations), "rows\n")
      cat("  efo_traits:", nrow(object@efo_traits), "rows\n")
      cat("  snp_allele:", nrow(object@snp_allele), "rows\n")
      cat("  bg_efo_traits:", nrow(object@bg_efo_traits), "rows\n")
      cat("  reported_traits:", nrow(object@reported_traits), "rows\n")
      cat("  locations:", nrow(object@locations), "rows\n")
      cat("  mapped_genes:", nrow(object@mapped_genes), "rows\n")
      cat("  snp_effect_alleles:", nrow(object@snp_effect_alleles), "rows\n")
    }
  }
)
