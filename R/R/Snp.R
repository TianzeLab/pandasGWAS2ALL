#' Snp class
#'
#' An S4 class that stores data of type Snp fetched from the GWAS API.
#' In Fat mode it holds 3 data.frames with hierarchical dependencies; in
#' Thin mode it only stores the raw list data.
#'
#' @slot raw_data list. Converted from obtained JSON data
#' @slot mode character. Fat or Thin
#' @slot snps data.frame
#' @slot locations data.frame
#' @slot mapped_genes data.frame
#'
#' @importFrom methods new
#'
#' @export
#'
setClass(
  "Snp",
  slots = c(
    raw_data = "list",
    mode = "character",
    snps = "data.frame",
    locations = "data.frame",
    mapped_genes = "data.frame"
  )
)

#' Snp constructor
#'
#' Creates a Snp object from a list of raw records. To understand the
#' significance of each column of the data.frames, visit "Snp" in the
#' GWAS Catalog documentation
#' (https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference).
#'
#' @param data A list of snp records
#' @param mode The mode of the object: "Fat" or "Thin"
#'
#' @return A Snp object
#'
#' @examples
#' \dontrun{
#' ch <- Snp(
#'   list(list(rs_id = "rs3093017", merged = TRUE)),
#'   mode = "Fat"
#' )
#' }
#'
#' @export
#'
Snp <- function(data = NULL, mode = "Fat") {
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
      if (!is.null(rec$rs_id)) {
        rec$rs_id <- as.character(rec$rs_id)
      }
      normalized_data <- c(normalized_data, list(rec))
    }
  }

  SNP_COLUMNS <- c(
    "rs_id",
    "merged",
    "functional_class",
    "current_snp",
    "last_update_date",
    "maf",
    "minor_allele",
    "alleles",
    "most_severe_consequence"
  )

  LOCATIONS_COLUMNS <- c("rs_id", "chromosome_name", "chromosome_position", "region.name")
  MAPPED_GENES_COLUMNS <- c("rs_id", "mapped_genes")

  NESTED_FIELDS <- c(
    "locations",
    "mapped_genes"
  )

  obj <- new("Snp")
  obj@mode <- mode

  obj@snps <- data.frame(matrix(ncol = length(SNP_COLUMNS), nrow = 0))
  colnames(obj@snps) <- SNP_COLUMNS

  obj@locations <- data.frame(matrix(ncol = length(LOCATIONS_COLUMNS), nrow = 0))
  colnames(obj@locations) <- LOCATIONS_COLUMNS

  obj@mapped_genes <- data.frame(matrix(ncol = length(MAPPED_GENES_COLUMNS), nrow = 0))
  colnames(obj@mapped_genes) <- MAPPED_GENES_COLUMNS

  if (mode == "Thin") {
    obj@raw_data <- normalized_data
    return(obj)
  }

  standardized_data <- list()
  for (item in normalized_data) {
    standardized <- list()
    for (col in SNP_COLUMNS) {
      standardized[[col]] <- item[[col]]
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

  snp_rows <- list()
  locations_rows <- list()
  mapped_genes_rows <- list()

  for (record in standardized_data) {
    rs_id <- record$rs_id

    snp_rows <- c(snp_rows, list(record[SNP_COLUMNS]))

    raw_locations <- record$locations
    if (is.null(raw_locations)) {
      iterable_locations <- list()
    } else if (is.list(raw_locations) && is.null(names(raw_locations))) {
      iterable_locations <- raw_locations
    } else {
      iterable_locations <- list(raw_locations)
    }
    for (item in iterable_locations) {
      if (is.list(item) && !is.null(names(item))) {
        raw_region <- item$region
        if (is.null(raw_region)) {
          region_name <- NA
        } else if (is.list(raw_region)) {
          region_name <- raw_region$name
        } else {
          region_name <- raw_region
        }
        locations_rows <- c(locations_rows, list(list(
          rs_id = rs_id,
          chromosome_name = item$chromosome_name,
          chromosome_position = item$chromosome_position,
          region.name = region_name
        )))
      } else {
        locations_rows <- c(locations_rows, list(list(
          rs_id = rs_id,
          chromosome_name = item,
          chromosome_position = NA,
          region.name = NA
        )))
      }
    }

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
        rs_id = rs_id,
        mapped_genes = item
      )))
    }
  }

  if (length(snp_rows) > 0) {
    built <- data.frame(row.names = seq_along(snp_rows))
    for (col in SNP_COLUMNS) {
      col_values <- vector("list", length(snp_rows))
      for (r_idx in seq_along(snp_rows)) {
        cell <- snp_rows[[r_idx]][[col]]
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
    obj@snps <- built
  }

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

  return(obj)
}

#' Show a Snp object
#'
#' Print a concise overview of a Snp object: the current mode and the
#' size of each data component.
#'
#' @param object A Snp object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "show",
  signature = "Snp",
  function(object) {
    if (object@mode == "Fat") {
      cat(
        "Snp is running in fat mode. It has 3 DataFrames with hierarchical dependencies.\n",
        "snps: ", nrow(object@snps), " rows\n",
        "|\n -locations: ", nrow(object@locations), " rows\n",
        "|\n -mapped_genes: ", nrow(object@mapped_genes), " rows\n",
        sep = ""
      )
    } else if (object@mode == "Thin") {
      cat(
        "Snp is running in thin mode. It has 1 list that contains the raw data.\n",
        "raw_data: a list of size ", length(object@raw_data), ".\n",
        sep = ""
      )
    }
  }
)

#' Length of a Snp object
#'
#' Get the number of snp records stored in the object.
#'
#' @param x A Snp object
#'
#' @return An integer, the number of records in the raw_data slot
#'
#' @export
setMethod(
  "length",
  signature = "Snp",
  function(x) {
    return(length(x@raw_data))
  }
)

#' Component names of a Snp object
#'
#' Get the names of the components (slots) of a Snp object.
#'
#' @param x A Snp object
#'
#' @return A character vector of component names
#'
#' @export
setMethod(
  "names",
  signature = "Snp",
  function(x) {
    return(c("raw_data", "mode", "snps", "locations", "mapped_genes"))
  }
)

#' Subset a Snp object
#'
#' Extract snp records by position or by rs identifier and return a new
#' Snp object of the same mode containing only the selected records.
#'
#' @param x A Snp object
#' @param i A numeric vector of positions (1-based) or a character vector of
#'   rs_ids. A list is unlisted before use. Leaving i missing returns x
#'   unchanged.
#'
#' @return A Snp object containing the selected records
#'
#' @examples
#' \dontrun{
#' a <- get_snps(rs_id = "rs3093017")
#' a[1]
#' a["rs3093017"]
#' }
#'
#' @export
setMethod(
  "[",
  signature = "Snp",
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
      if (!is.null(raw_data[[j]]$rs_id)) {
        raw_data_dict[[raw_data[[j]]$rs_id]] <- raw_data[[j]]
      }
    }

    sub_set <- list()
    for (i_val in arr) {
      if (is.character(i_val)) {
        if (!i_val %in% names(raw_data_dict)) {
          stop("Unknown rs_id: ", i_val)
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

    return(Snp(sub_set, x@mode))
  }
)

#' Set operations on Snp objects
#'
#' Combine or compare two Snp objects of the same mode as sets of snp
#' records, identified and deduplicated by rs_id.
#'
#' @param e1 A Snp object
#' @param e2 A Snp object
#'
#' @return A Snp object of the same mode
#'
#' @examples
#' \dontrun{
#' a <- get_snps(rs_id = "rs3093017")
#' b <- get_snps(rs_id = "rs7671745")
#' a + b
#' a - b
#' a & b
#' a | b
#' a ^ b
#' }
#'
#' @rdname Snp_set_operations
#' @export
setMethod(
  "+",
  signature = c("Snp", "Snp"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      return(Snp(c(e1@raw_data, e2@raw_data), e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Snp_set_operations
#' @export
setMethod(
  "-",
  signature = c("Snp", "Snp"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$rs_id)) {
          self_key_set <- unique(c(self_key_set, i$rs_id))
          self_dict[[i$rs_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$rs_id)) {
          other_key_set <- unique(c(other_key_set, j$rs_id))
        }
      }

      sub_key <- self_key_set[!self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Snp(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Snp_set_operations
#' @export
setMethod(
  "&",
  signature = c("Snp", "Snp"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$rs_id)) {
          self_key_set <- unique(c(self_key_set, i$rs_id))
          self_dict[[i$rs_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$rs_id)) {
          other_key_set <- unique(c(other_key_set, j$rs_id))
        }
      }

      sub_key <- self_key_set[self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Snp(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Snp_set_operations
#' @export
setMethod(
  "|",
  signature = c("Snp", "Snp"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      and_dict <- list()

      for (i in e1@raw_data) {
        if (!is.null(i$rs_id)) {
          and_dict[[i$rs_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$rs_id)) {
          and_dict[[j$rs_id]] <- j
        }
      }

      data <- unname(and_dict)
      return(Snp(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Snp_set_operations
#' @export
setMethod(
  "^",
  signature = c("Snp", "Snp"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      and_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$rs_id)) {
          self_key_set <- unique(c(self_key_set, i$rs_id))
          and_dict[[i$rs_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$rs_id)) {
          other_key_set <- unique(c(other_key_set, j$rs_id))
          and_dict[[j$rs_id]] <- j
        }
      }

      all_keys <- unique(c(self_key_set, other_key_set))
      common_keys <- self_key_set[self_key_set %in% other_key_set]
      sub_key <- all_keys[!all_keys %in% common_keys]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(and_dict[[k]]))
      }

      return(Snp(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' Compare Snp objects for equality
#'
#' Test whether two Snp objects are equal. Two objects are equal when
#' they run in the same mode and store identical records.
#'
#' @param e1 A Snp object
#' @param e2 A Snp object
#'
#' @return A logical value
#'
#' @export
setMethod(
  "==",
  signature = c("Snp", "Snp"),
  function(e1, e2) {
    if (is.null(e1) || is.null(e2)) {
      return(is.null(e1) && is.null(e2))
    }
    return(identical(e1@raw_data, e2@raw_data) && e1@mode == e2@mode)
  }
)

#' Summarize a Snp object
#'
#' Print a summary of a Snp object: the mode, the number of raw records,
#' the number of snps and, in Fat mode, the row counts of all 3
#' data.frames.
#'
#' @param object A Snp object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "summary",
  signature = "Snp",
  function(object) {
    cat("Snp Object Summary\n")
    cat("==================\n")
    cat("Mode:", object@mode, "\n")
    cat("Number of raw records:", length(object@raw_data), "\n")
    cat("Number of snps:", nrow(object@snps), "\n")
    if (object@mode == "Fat") {
      cat("\nDataFrames:\n")
      cat("  snps:", nrow(object@snps), "rows\n")
      cat("  locations:", nrow(object@locations), "rows\n")
      cat("  mapped_genes:", nrow(object@mapped_genes), "rows\n")
    }
  }
)
