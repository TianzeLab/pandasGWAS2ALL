#' Trait class
#'
#' An S4 class that stores data of type Trait fetched from the GWAS API.
#' In Fat mode it holds 1 data.frame; in Thin mode it only stores the raw
#' list data.
#'
#' @slot raw_data list. Converted from obtained JSON data
#' @slot mode character. Fat or Thin
#' @slot traits data.frame
#'
#' @importFrom methods new
#'
#' @export
#'
setClass(
  "Trait",
  slots = c(
    raw_data = "list",
    mode = "character",
    traits = "data.frame"
  )
)

#' Trait constructor
#'
#' Creates a Trait object from a list of raw records. To understand the
#' significance of each column of the data.frame, visit "Trait" in the
#' GWAS Catalog documentation
#' (https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference).
#'
#' @param data A list of trait records
#' @param mode The mode of the object: "Fat" or "Thin"
#'
#' @return A Trait object
#'
#' @examples
#' \dontrun{
#' ch <- Trait(
#'   list(list(efo_id = "EFO_0009713", efo_trait = "asthma")),
#'   mode = "Fat"
#' )
#' }
#'
#' @export
#'
Trait <- function(data = NULL, mode = "Fat") {
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
      if (!is.null(rec$efo_id)) {
        rec$efo_id <- as.character(rec$efo_id)
      }
      normalized_data <- c(normalized_data, list(rec))
    }
  }

  TRAIT_COLUMNS <- c(
    "efo_trait",
    "uri",
    "efo_id"
  )

  obj <- new("Trait")
  obj@mode <- mode

  obj@traits <- data.frame(matrix(ncol = length(TRAIT_COLUMNS), nrow = 0))
  colnames(obj@traits) <- TRAIT_COLUMNS

  if (mode == "Thin") {
    obj@raw_data <- normalized_data
    return(obj)
  }

  standardized_data <- list()
  for (item in normalized_data) {
    standardized <- list()
    for (col in TRAIT_COLUMNS) {
      standardized[[col]] <- item[[col]]
    }
    standardized_data <- c(standardized_data, list(standardized))
  }

  obj@raw_data <- standardized_data

  if (length(standardized_data) == 0) {
    return(obj)
  }

  trait_rows <- list()

  for (record in standardized_data) {
    trait_rows <- c(trait_rows, list(record[TRAIT_COLUMNS]))
  }

  if (length(trait_rows) > 0) {
    built <- data.frame(row.names = seq_along(trait_rows))
    for (col in TRAIT_COLUMNS) {
      col_values <- vector("list", length(trait_rows))
      for (r_idx in seq_along(trait_rows)) {
        cell <- trait_rows[[r_idx]][[col]]
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
    obj@traits <- built
  }

  return(obj)
}

#' Show a Trait object
#'
#' Print a concise overview of a Trait object: the current mode and the
#' size of each data component. In Fat mode the row count of the
#' data.frame is listed; in Thin mode the size of the raw record list is
#' shown.
#'
#' @param object A Trait object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "show",
  signature = "Trait",
  function(object) {
    if (object@mode == "Fat") {
      cat(
        "Trait is running in fat mode. It has 1 DataFrame with hierarchical dependencies.\n",
        "traits: ", nrow(object@traits), " rows\n",
        sep = ""
      )
    } else if (object@mode == "Thin") {
      cat(
        "Trait is running in thin mode. It has 1 list that contains the raw data.\n",
        "raw_data: a list of size ", length(object@raw_data), ".\n",
        sep = ""
      )
    }
  }
)

#' Length of a Trait object
#'
#' Get the number of trait records stored in the object.
#'
#' @param x A Trait object
#'
#' @return An integer, the number of records in the raw_data slot
#'
#' @export
setMethod(
  "length",
  signature = "Trait",
  function(x) {
    return(length(x@raw_data))
  }
)

#' Component names of a Trait object
#'
#' Get the names of the components (slots) of a Trait object.
#'
#' @param x A Trait object
#'
#' @return A character vector of component names
#'
#' @export
setMethod(
  "names",
  signature = "Trait",
  function(x) {
    return(c("raw_data", "mode", "traits"))
  }
)

#' Subset a Trait object
#'
#' Extract trait records by position or by efo identifier and return a new
#' Trait object of the same mode containing only the selected records.
#'
#' @param x A Trait object
#' @param i A numeric vector of positions (1-based) or a character vector of
#'   efo_ids. A list is unlisted before use, so a list of positions or
#'   identifiers is also accepted. Leaving i missing returns x unchanged.
#'
#' @return A Trait object containing the selected records
#'
#' @examples
#' \dontrun{
#' a <- get_traits(efo_id = "EFO_0009713")
#' a[1]
#' a["EFO_0009713"]
#' }
#'
#' @export
setMethod(
  "[",
  signature = "Trait",
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
      if (!is.null(raw_data[[j]]$efo_id)) {
        raw_data_dict[[raw_data[[j]]$efo_id]] <- raw_data[[j]]
      }
    }

    sub_set <- list()
    for (i_val in arr) {
      if (is.character(i_val)) {
        if (!i_val %in% names(raw_data_dict)) {
          stop("Unknown efo_id: ", i_val)
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

    return(Trait(sub_set, x@mode))
  }
)

#' Set operations on Trait objects
#'
#' Combine or compare two Trait objects of the same mode as sets of trait
#' records, identified and deduplicated by efo_id.
#'
#' `+` concatenates the records of both objects, in the order of e1 then
#' e2. `-` returns the records of e1 whose efo_id does not appear in e2
#' (difference). `&` returns the records of e1 whose efo_id also appears
#' in e2 (intersection). `|` returns the records of both objects
#' deduplicated by efo_id, keeping the record of e1 when both objects
#' contain the same id (union). `^` returns the records whose efo_id
#' appears in exactly one of the two objects (symmetric difference). All
#' operators raise an error when the two objects do not run in the same
#' mode.
#'
#' @param e1 A Trait object
#' @param e2 A Trait object
#'
#' @return A Trait object of the same mode
#'
#' @examples
#' \dontrun{
#' a <- get_traits(efo_id = "EFO_0009713")
#' b <- get_traits(efo_id = "OBA_VT0001321")
#' a + b
#' a - b
#' a & b
#' a | b
#' a ^ b
#' }
#'
#' @rdname Trait_set_operations
#' @export
setMethod(
  "+",
  signature = c("Trait", "Trait"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      return(Trait(c(e1@raw_data, e2@raw_data), e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Trait_set_operations
#' @export
setMethod(
  "-",
  signature = c("Trait", "Trait"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$efo_id)) {
          self_key_set <- unique(c(self_key_set, i$efo_id))
          self_dict[[i$efo_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$efo_id)) {
          other_key_set <- unique(c(other_key_set, j$efo_id))
        }
      }

      sub_key <- self_key_set[!self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Trait(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Trait_set_operations
#' @export
setMethod(
  "&",
  signature = c("Trait", "Trait"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$efo_id)) {
          self_key_set <- unique(c(self_key_set, i$efo_id))
          self_dict[[i$efo_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$efo_id)) {
          other_key_set <- unique(c(other_key_set, j$efo_id))
        }
      }

      sub_key <- self_key_set[self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Trait(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Trait_set_operations
#' @export
setMethod(
  "|",
  signature = c("Trait", "Trait"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      and_dict <- list()

      for (i in e1@raw_data) {
        if (!is.null(i$efo_id)) {
          and_dict[[i$efo_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$efo_id)) {
          and_dict[[j$efo_id]] <- j
        }
      }

      data <- unname(and_dict)
      return(Trait(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Trait_set_operations
#' @export
setMethod(
  "^",
  signature = c("Trait", "Trait"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      and_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$efo_id)) {
          self_key_set <- unique(c(self_key_set, i$efo_id))
          and_dict[[i$efo_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$efo_id)) {
          other_key_set <- unique(c(other_key_set, j$efo_id))
          and_dict[[j$efo_id]] <- j
        }
      }

      all_keys <- unique(c(self_key_set, other_key_set))
      common_keys <- self_key_set[self_key_set %in% other_key_set]
      sub_key <- all_keys[!all_keys %in% common_keys]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(and_dict[[k]]))
      }

      return(Trait(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' Compare Trait objects for equality
#'
#' Test whether two Trait objects are equal. Two objects are equal when
#' they run in the same mode and store identical records.
#'
#' @param e1 A Trait object
#' @param e2 A Trait object
#'
#' @return A logical value
#'
#' @export
setMethod(
  "==",
  signature = c("Trait", "Trait"),
  function(e1, e2) {
    if (is.null(e1) || is.null(e2)) {
      return(is.null(e1) && is.null(e2))
    }
    return(identical(e1@raw_data, e2@raw_data) && e1@mode == e2@mode)
  }
)

#' Summarize a Trait object
#'
#' Print a summary of a Trait object: the mode, the number of raw records,
#' the number of traits and, in Fat mode, the row count of the data.frame.
#'
#' @param object A Trait object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "summary",
  signature = "Trait",
  function(object) {
    cat("Trait Object Summary\n")
    cat("====================\n")
    cat("Mode:", object@mode, "\n")
    cat("Number of raw records:", length(object@raw_data), "\n")
    cat("Number of traits:", nrow(object@traits), "\n")
    if (object@mode == "Fat") {
      cat("\nDataFrames:\n")
      cat("  traits:", nrow(object@traits), "rows\n")
    }
  }
)
