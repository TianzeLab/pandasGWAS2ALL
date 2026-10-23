#' Ancestry class
#'
#' An S4 class that stores data of type Ancestry fetched from the GWAS
#' API. In Fat mode it holds 4 data.frames with hierarchical dependencies;
#' in Thin mode it only stores the raw list data.
#'
#' @slot raw_data list. Converted from obtained JSON data
#' @slot mode character. Fat or Thin
#' @slot ancestries data.frame
#' @slot ancestral_groups data.frame
#' @slot country_of_origins data.frame
#' @slot country_of_recruitments data.frame
#'
#' @importFrom methods new
#'
#' @export
#'
setClass(
  "Ancestry",
  slots = c(
    raw_data = "list",
    mode = "character",
    ancestries = "data.frame",
    ancestral_groups = "data.frame",
    country_of_origins = "data.frame",
    country_of_recruitments = "data.frame"
  )
)

#' Ancestry constructor
#'
#' Creates an Ancestry object from a list of raw records. To understand
#' the significance of each column of the data.frames, visit "Ancestry"
#' in the GWAS Catalog documentation
#' (https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference).
#'
#' @param data A list of ancestry records
#' @param mode The mode of the object: "Fat" or "Thin"
#'
#' @return An Ancestry object
#'
#' @examples
#' \dontrun{
#' ch <- Ancestry(
#'   list(list(accession_id = "GCST000854", ancestry_id = "4595")),
#'   mode = "Fat"
#' )
#' }
#'
#' @export
#'
Ancestry <- function(data = NULL, mode = "Fat") {
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
      if (!is.null(rec$accession_id)) {
        rec$accession_id <- as.character(rec$accession_id)
      }
      if (!is.null(rec$ancestry_id)) {
        rec$ancestry_id <- as.character(rec$ancestry_id)
      }
      normalized_data <- c(normalized_data, list(rec))
    }
  }

  ANCESTRY_COLUMNS <- c(
    "accession_id",
    "ancestry_id",
    "type",
    "number_of_individuals"
  )

  ANCESTRAL_GROUPS_COLUMNS <- c("accession_id", "ancestry_id", "ancestral_group")
  COUNTRY_OF_ORIGINS_COLUMNS <- c("accession_id", "ancestry_id", "major_area", "region", "country_name")
  COUNTRY_OF_RECRUITMENTS_COLUMNS <- c("accession_id", "ancestry_id", "major_area", "region", "country_name")

  NESTED_FIELDS <- c(
    "ancestral_groups",
    "country_of_origin",
    "country_of_recruitment"
  )

  obj <- new("Ancestry")
  obj@mode <- mode

  obj@ancestries <- data.frame(matrix(ncol = length(ANCESTRY_COLUMNS), nrow = 0))
  colnames(obj@ancestries) <- ANCESTRY_COLUMNS

  obj@ancestral_groups <- data.frame(matrix(ncol = length(ANCESTRAL_GROUPS_COLUMNS), nrow = 0))
  colnames(obj@ancestral_groups) <- ANCESTRAL_GROUPS_COLUMNS

  obj@country_of_origins <- data.frame(matrix(ncol = length(COUNTRY_OF_ORIGINS_COLUMNS), nrow = 0))
  colnames(obj@country_of_origins) <- COUNTRY_OF_ORIGINS_COLUMNS

  obj@country_of_recruitments <- data.frame(matrix(ncol = length(COUNTRY_OF_RECRUITMENTS_COLUMNS), nrow = 0))
  colnames(obj@country_of_recruitments) <- COUNTRY_OF_RECRUITMENTS_COLUMNS

  if (mode == "Thin") {
    obj@raw_data <- normalized_data
    return(obj)
  }

  standardized_data <- list()
  for (item in normalized_data) {
    standardized <- list()
    for (col in ANCESTRY_COLUMNS) {
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

  ancestry_rows <- list()
  ancestral_groups_rows <- list()
  country_of_origins_rows <- list()
  country_of_recruitments_rows <- list()

  for (record in standardized_data) {
    accession_id <- record$accession_id
    ancestry_id <- record$ancestry_id

    ancestry_rows <- c(ancestry_rows, list(record[ANCESTRY_COLUMNS]))

    raw_ancestral_groups <- record$ancestral_groups
    if (is.null(raw_ancestral_groups)) {
      iterable_ancestral_groups <- list()
    } else if (is.list(raw_ancestral_groups) && is.null(names(raw_ancestral_groups))) {
      iterable_ancestral_groups <- raw_ancestral_groups
    } else {
      iterable_ancestral_groups <- list(raw_ancestral_groups)
    }
    for (item in iterable_ancestral_groups) {
      if (is.list(item) && !is.null(names(item))) {
        ancestral_groups_rows <- c(ancestral_groups_rows, list(list(
          accession_id = accession_id,
          ancestry_id = ancestry_id,
          ancestral_group = item$ancestral_group
        )))
      } else {
        ancestral_groups_rows <- c(ancestral_groups_rows, list(list(
          accession_id = accession_id,
          ancestry_id = ancestry_id,
          ancestral_group = item
        )))
      }
    }

    raw_country_of_origins <- record$country_of_origin
    if (is.null(raw_country_of_origins)) {
      iterable_country_of_origins <- list()
    } else if (is.list(raw_country_of_origins) && is.null(names(raw_country_of_origins))) {
      iterable_country_of_origins <- raw_country_of_origins
    } else {
      iterable_country_of_origins <- list(raw_country_of_origins)
    }
    for (item in iterable_country_of_origins) {
      if (is.list(item) && !is.null(names(item))) {
        country_of_origins_rows <- c(country_of_origins_rows, list(list(
          accession_id = accession_id,
          ancestry_id = ancestry_id,
          major_area = item$major_area,
          region = item$region,
          country_name = item$country_name
        )))
      } else {
        country_of_origins_rows <- c(country_of_origins_rows, list(list(
          accession_id = accession_id,
          ancestry_id = ancestry_id,
          major_area = NA,
          region = NA,
          country_name = item
        )))
      }
    }

    raw_country_of_recruitments <- record$country_of_recruitment
    if (is.null(raw_country_of_recruitments)) {
      iterable_country_of_recruitments <- list()
    } else if (is.list(raw_country_of_recruitments) && is.null(names(raw_country_of_recruitments))) {
      iterable_country_of_recruitments <- raw_country_of_recruitments
    } else {
      iterable_country_of_recruitments <- list(raw_country_of_recruitments)
    }
    for (item in iterable_country_of_recruitments) {
      if (is.list(item) && !is.null(names(item))) {
        country_of_recruitments_rows <- c(country_of_recruitments_rows, list(list(
          accession_id = accession_id,
          ancestry_id = ancestry_id,
          major_area = item$major_area,
          region = item$region,
          country_name = item$country_name
        )))
      } else {
        country_of_recruitments_rows <- c(country_of_recruitments_rows, list(list(
          accession_id = accession_id,
          ancestry_id = ancestry_id,
          major_area = NA,
          region = NA,
          country_name = item
        )))
      }
    }
  }

  if (length(ancestry_rows) > 0) {
    built <- data.frame(row.names = seq_along(ancestry_rows))
    for (col in ANCESTRY_COLUMNS) {
      col_values <- vector("list", length(ancestry_rows))
      for (r_idx in seq_along(ancestry_rows)) {
        cell <- ancestry_rows[[r_idx]][[col]]
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
    obj@ancestries <- built
  }

  if (length(ancestral_groups_rows) > 0) {
    built <- data.frame(row.names = seq_along(ancestral_groups_rows))
    for (col in ANCESTRAL_GROUPS_COLUMNS) {
      col_values <- vector("list", length(ancestral_groups_rows))
      for (r_idx in seq_along(ancestral_groups_rows)) {
        cell <- ancestral_groups_rows[[r_idx]][[col]]
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
    obj@ancestral_groups <- built
  }

  if (length(country_of_origins_rows) > 0) {
    built <- data.frame(row.names = seq_along(country_of_origins_rows))
    for (col in COUNTRY_OF_ORIGINS_COLUMNS) {
      col_values <- vector("list", length(country_of_origins_rows))
      for (r_idx in seq_along(country_of_origins_rows)) {
        cell <- country_of_origins_rows[[r_idx]][[col]]
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
    obj@country_of_origins <- built
  }

  if (length(country_of_recruitments_rows) > 0) {
    built <- data.frame(row.names = seq_along(country_of_recruitments_rows))
    for (col in COUNTRY_OF_RECRUITMENTS_COLUMNS) {
      col_values <- vector("list", length(country_of_recruitments_rows))
      for (r_idx in seq_along(country_of_recruitments_rows)) {
        cell <- country_of_recruitments_rows[[r_idx]][[col]]
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
    obj@country_of_recruitments <- built
  }

  return(obj)
}

#' Show an Ancestry object
#'
#' Print a concise overview of an Ancestry object: the current mode and
#' the size of each data component.
#'
#' @param object An Ancestry object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "show",
  signature = "Ancestry",
  function(object) {
    if (object@mode == "Fat") {
      cat(
        "Ancestry is running in fat mode. It has 4 DataFrames with hierarchical dependencies.\n",
        "ancestries: ", nrow(object@ancestries), " rows\n",
        "|\n -ancestral_groups: ", nrow(object@ancestral_groups), " rows\n",
        "|\n -country_of_origins: ", nrow(object@country_of_origins), " rows\n",
        "|\n -country_of_recruitments: ", nrow(object@country_of_recruitments), " rows\n",
        sep = ""
      )
    } else if (object@mode == "Thin") {
      cat(
        "Ancestry is running in thin mode. It has 1 list that contains the raw data.\n",
        "raw_data: a list of size ", length(object@raw_data), ".\n",
        sep = ""
      )
    }
  }
)

#' Length of an Ancestry object
#'
#' Get the number of ancestry records stored in the object.
#'
#' @param x An Ancestry object
#'
#' @return An integer, the number of records in the raw_data slot
#'
#' @export
setMethod(
  "length",
  signature = "Ancestry",
  function(x) {
    return(length(x@raw_data))
  }
)

#' Component names of an Ancestry object
#'
#' Get the names of the components (slots) of an Ancestry object.
#'
#' @param x An Ancestry object
#'
#' @return A character vector of component names
#'
#' @export
setMethod(
  "names",
  signature = "Ancestry",
  function(x) {
    return(c(
      "raw_data",
      "mode",
      "ancestries",
      "ancestral_groups",
      "country_of_origins",
      "country_of_recruitments"
    ))
  }
)

#' Subset an Ancestry object
#'
#' Extract ancestry records by position or by the composite identifier
#' "accession_id|ancestry_id" and return a new Ancestry object of the same
#' mode containing only the selected records.
#'
#' @param x An Ancestry object
#' @param i A numeric vector of positions (1-based) or a character vector of
#'   composite identifiers in the form "accession_id|ancestry_id". A list
#'   is unlisted before use. Leaving i missing returns x unchanged.
#'
#' @return An Ancestry object containing the selected records
#'
#' @examples
#' \dontrun{
#' a <- get_ancestries(accession_id = "GCST000854", ancestry_id = "4595")
#' a[1]
#' a["GCST000854|4595"]
#' }
#'
#' @export
setMethod(
  "[",
  signature = "Ancestry",
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
      if (!is.null(raw_data[[j]]$accession_id) &&
          !is.null(raw_data[[j]]$ancestry_id)) {
        key <- paste0(raw_data[[j]]$accession_id, "|", raw_data[[j]]$ancestry_id)
        raw_data_dict[[key]] <- raw_data[[j]]
      }
    }

    sub_set <- list()
    for (i_val in arr) {
      if (is.character(i_val)) {
        if (!i_val %in% names(raw_data_dict)) {
          stop("Unknown accession_id|ancestry_id: ", i_val)
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

    return(Ancestry(sub_set, x@mode))
  }
)

#' Set operations on Ancestry objects
#'
#' Combine or compare two Ancestry objects of the same mode as sets of
#' ancestry records, identified and deduplicated by the composite key
#' "accession_id|ancestry_id".
#'
#' @param e1 An Ancestry object
#' @param e2 An Ancestry object
#'
#' @return An Ancestry object of the same mode
#'
#' @examples
#' \dontrun{
#' a <- get_ancestries(accession_id = "GCST000854", ancestry_id = "4595")
#' b <- get_ancestries(accession_id = "GCST000854", ancestry_id = "7528")
#' a + b
#' a - b
#' a & b
#' a | b
#' a ^ b
#' }
#'
#' @rdname Ancestry_set_operations
#' @export
setMethod(
  "+",
  signature = c("Ancestry", "Ancestry"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      return(Ancestry(c(e1@raw_data, e2@raw_data), e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Ancestry_set_operations
#' @export
setMethod(
  "-",
  signature = c("Ancestry", "Ancestry"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$accession_id) && !is.null(i$ancestry_id)) {
          key <- paste0(i$accession_id, "|", i$ancestry_id)
          self_key_set <- unique(c(self_key_set, key))
          self_dict[[key]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$accession_id) && !is.null(j$ancestry_id)) {
          other_key_set <- unique(c(other_key_set, paste0(j$accession_id, "|", j$ancestry_id)))
        }
      }

      sub_key <- self_key_set[!self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Ancestry(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Ancestry_set_operations
#' @export
setMethod(
  "&",
  signature = c("Ancestry", "Ancestry"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$accession_id) && !is.null(i$ancestry_id)) {
          key <- paste0(i$accession_id, "|", i$ancestry_id)
          self_key_set <- unique(c(self_key_set, key))
          self_dict[[key]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$accession_id) && !is.null(j$ancestry_id)) {
          other_key_set <- unique(c(other_key_set, paste0(j$accession_id, "|", j$ancestry_id)))
        }
      }

      sub_key <- self_key_set[self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Ancestry(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Ancestry_set_operations
#' @export
setMethod(
  "|",
  signature = c("Ancestry", "Ancestry"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      and_dict <- list()

      for (i in e1@raw_data) {
        if (!is.null(i$accession_id) && !is.null(i$ancestry_id)) {
          and_dict[[paste0(i$accession_id, "|", i$ancestry_id)]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$accession_id) && !is.null(j$ancestry_id)) {
          and_dict[[paste0(j$accession_id, "|", j$ancestry_id)]] <- j
        }
      }

      data <- unname(and_dict)
      return(Ancestry(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Ancestry_set_operations
#' @export
setMethod(
  "^",
  signature = c("Ancestry", "Ancestry"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      and_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$accession_id) && !is.null(i$ancestry_id)) {
          key <- paste0(i$accession_id, "|", i$ancestry_id)
          self_key_set <- unique(c(self_key_set, key))
          and_dict[[key]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$accession_id) && !is.null(j$ancestry_id)) {
          key <- paste0(j$accession_id, "|", j$ancestry_id)
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

      return(Ancestry(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' Compare Ancestry objects for equality
#'
#' Test whether two Ancestry objects are equal. Two objects are equal
#' when they run in the same mode and store identical records.
#'
#' @param e1 An Ancestry object
#' @param e2 An Ancestry object
#'
#' @return A logical value
#'
#' @export
setMethod(
  "==",
  signature = c("Ancestry", "Ancestry"),
  function(e1, e2) {
    if (is.null(e1) || is.null(e2)) {
      return(is.null(e1) && is.null(e2))
    }
    return(identical(e1@raw_data, e2@raw_data) && e1@mode == e2@mode)
  }
)

#' Summarize an Ancestry object
#'
#' Print a summary of an Ancestry object: the mode, the number of raw
#' records, the number of ancestries and, in Fat mode, the row counts of
#' all 4 data.frames.
#'
#' @param object An Ancestry object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "summary",
  signature = "Ancestry",
  function(object) {
    cat("Ancestry Object Summary\n")
    cat("=======================\n")
    cat("Mode:", object@mode, "\n")
    cat("Number of raw records:", length(object@raw_data), "\n")
    cat("Number of ancestries:", nrow(object@ancestries), "\n")
    if (object@mode == "Fat") {
      cat("\nDataFrames:\n")
      cat("  ancestries:", nrow(object@ancestries), "rows\n")
      cat("  ancestral_groups:", nrow(object@ancestral_groups), "rows\n")
      cat("  country_of_origins:", nrow(object@country_of_origins), "rows\n")
      cat("  country_of_recruitments:", nrow(object@country_of_recruitments), "rows\n")
    }
  }
)
