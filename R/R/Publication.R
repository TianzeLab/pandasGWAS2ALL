#' Publication class
#'
#' An S4 class that stores data of type Publication fetched from the GWAS
#' API. In Fat mode it holds 2 data.frames with hierarchical dependencies;
#' in Thin mode it only stores the raw list data.
#'
#' @slot raw_data list. Converted from obtained JSON data
#' @slot mode character. Fat or Thin
#' @slot publications data.frame
#' @slot authors data.frame
#'
#' @importFrom methods new
#'
#' @export
#'
setClass(
  "Publication",
  slots = c(
    raw_data = "list",
    mode = "character",
    publications = "data.frame",
    authors = "data.frame"
  )
)

#' Publication constructor
#'
#' Creates a Publication object from a list of raw records. To understand
#' the significance of each column of the data.frames, visit "Publication"
#' in the GWAS Catalog documentation
#' (https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference).
#'
#' @param data A list of publication records
#' @param mode The mode of the object: "Fat" or "Thin"
#'
#' @return A Publication object
#'
#' @examples
#' \dontrun{
#' ch <- Publication(
#'   list(list(pubmed_id = "35241825", journal = "Nature")),
#'   mode = "Fat"
#' )
#' }
#'
#' @export
#'
Publication <- function(data = NULL, mode = "Fat") {
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
      if (!is.null(rec$pubmed_id)) {
        rec$pubmed_id <- as.character(rec$pubmed_id)
      }
      normalized_data <- c(normalized_data, list(rec))
    }
  }

  PUBLICATION_COLUMNS <- c(
    "pubmed_id",
    "publication_date",
    "journal",
    "title",
    "first_author.full_name",
    "first_author.full_name_standard",
    "first_author.first_name",
    "first_author.last_name",
    "first_author.initials",
    "first_author.affiliation",
    "first_author.orcid"
  )

  AUTHORS_COLUMNS <- c(
    "pubmed_id",
    "full_name",
    "full_name_standard",
    "first_name",
    "last_name",
    "initials",
    "affiliation",
    "orcid"
  )

  NESTED_FIELDS <- c(
    "authors"
  )

  obj <- new("Publication")
  obj@mode <- mode

  obj@publications <- data.frame(matrix(ncol = length(PUBLICATION_COLUMNS), nrow = 0))
  colnames(obj@publications) <- PUBLICATION_COLUMNS

  obj@authors <- data.frame(matrix(ncol = length(AUTHORS_COLUMNS), nrow = 0))
  colnames(obj@authors) <- AUTHORS_COLUMNS

  if (mode == "Thin") {
    obj@raw_data <- normalized_data
    return(obj)
  }

  standardized_data <- list()
  for (item in normalized_data) {
    raw_first_author <- item$first_author
    if (is.list(raw_first_author) && is.null(names(raw_first_author))) {
      if (length(raw_first_author) == 0) {
        first_author <- NULL
      } else {
        first_author <- raw_first_author[[1]]
      }
    } else {
      first_author <- raw_first_author
    }

    standardized <- list()
    for (col in PUBLICATION_COLUMNS) {
      if (startsWith(col, "first_author.")) {
        author_col <- sub("^first_author\\.", "", col)
        if (col %in% names(item)) {
          value <- item[[col]]
        } else if (is.list(first_author) && !is.null(names(first_author))) {
          value <- first_author[[author_col]]
        } else if (author_col == "full_name") {
          value <- first_author
        } else {
          value <- NULL
        }
      } else {
        value <- item[[col]]
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

  publication_rows <- list()
  authors_rows <- list()

  for (record in standardized_data) {
    pubmed_id <- record$pubmed_id

    publication_rows <- c(publication_rows, list(record[PUBLICATION_COLUMNS]))

    raw_authors <- record$authors
    if (is.null(raw_authors)) {
      iterable_authors <- list()
    } else if (is.list(raw_authors) && is.null(names(raw_authors))) {
      iterable_authors <- raw_authors
    } else {
      iterable_authors <- list(raw_authors)
    }
    for (item in iterable_authors) {
      if (is.list(item) && !is.null(names(item))) {
        authors_rows <- c(authors_rows, list(list(
          pubmed_id = pubmed_id,
          full_name = item$full_name,
          full_name_standard = item$full_name_standard,
          first_name = item$first_name,
          last_name = item$last_name,
          initials = item$initials,
          affiliation = item$affiliation,
          orcid = item$orcid
        )))
      } else {
        authors_rows <- c(authors_rows, list(list(
          pubmed_id = pubmed_id,
          full_name = item,
          full_name_standard = NA,
          first_name = NA,
          last_name = NA,
          initials = NA,
          affiliation = NA,
          orcid = NA
        )))
      }
    }
  }

  if (length(publication_rows) > 0) {
    built <- data.frame(row.names = seq_along(publication_rows))
    for (col in PUBLICATION_COLUMNS) {
      col_values <- vector("list", length(publication_rows))
      for (r_idx in seq_along(publication_rows)) {
        cell <- publication_rows[[r_idx]][[col]]
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
    obj@publications <- built
  }

  if (length(authors_rows) > 0) {
    built <- data.frame(row.names = seq_along(authors_rows))
    for (col in AUTHORS_COLUMNS) {
      col_values <- vector("list", length(authors_rows))
      for (r_idx in seq_along(authors_rows)) {
        cell <- authors_rows[[r_idx]][[col]]
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
    obj@authors <- built
  }

  return(obj)
}

#' Show a Publication object
#'
#' Print a concise overview of a Publication object: the current mode and
#' the size of each data component.
#'
#' @param object A Publication object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "show",
  signature = "Publication",
  function(object) {
    if (object@mode == "Fat") {
      cat(
        "Publication is running in fat mode. It has 2 DataFrames with hierarchical dependencies.\n",
        "publications: ", nrow(object@publications), " rows\n",
        "|\n -authors: ", nrow(object@authors), " rows\n",
        sep = ""
      )
    } else if (object@mode == "Thin") {
      cat(
        "Publication is running in thin mode. It has 1 list that contains the raw data.\n",
        "raw_data: a list of size ", length(object@raw_data), ".\n",
        sep = ""
      )
    }
  }
)

#' Length of a Publication object
#'
#' Get the number of publication records stored in the object.
#'
#' @param x A Publication object
#'
#' @return An integer, the number of records in the raw_data slot
#'
#' @export
setMethod(
  "length",
  signature = "Publication",
  function(x) {
    return(length(x@raw_data))
  }
)

#' Component names of a Publication object
#'
#' Get the names of the components (slots) of a Publication object.
#'
#' @param x A Publication object
#'
#' @return A character vector of component names
#'
#' @export
setMethod(
  "names",
  signature = "Publication",
  function(x) {
    return(c("raw_data", "mode", "publications", "authors"))
  }
)

#' Subset a Publication object
#'
#' Extract publication records by position or by pubmed identifier and
#' return a new Publication object of the same mode containing only the
#' selected records.
#'
#' @param x A Publication object
#' @param i A numeric vector of positions (1-based) or a character vector of
#'   pubmed_ids. A list is unlisted before use. Leaving i missing returns x
#'   unchanged.
#'
#' @return A Publication object containing the selected records
#'
#' @examples
#' \dontrun{
#' a <- get_publications(pubmed_id = "35241825")
#' a[1]
#' a["35241825"]
#' }
#'
#' @export
setMethod(
  "[",
  signature = "Publication",
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
      if (!is.null(raw_data[[j]]$pubmed_id)) {
        raw_data_dict[[raw_data[[j]]$pubmed_id]] <- raw_data[[j]]
      }
    }

    sub_set <- list()
    for (i_val in arr) {
      if (is.character(i_val)) {
        if (!i_val %in% names(raw_data_dict)) {
          stop("Unknown pubmed_id: ", i_val)
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

    return(Publication(sub_set, x@mode))
  }
)

#' Set operations on Publication objects
#'
#' Combine or compare two Publication objects of the same mode as sets of
#' publication records, identified and deduplicated by pubmed_id.
#'
#' @param e1 A Publication object
#' @param e2 A Publication object
#'
#' @return A Publication object of the same mode
#'
#' @examples
#' \dontrun{
#' a <- get_publications(pubmed_id = "35241825")
#' b <- get_publications(pubmed_id = "42127933")
#' a + b
#' a - b
#' a & b
#' a | b
#' a ^ b
#' }
#'
#' @rdname Publication_set_operations
#' @export
setMethod(
  "+",
  signature = c("Publication", "Publication"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      return(Publication(c(e1@raw_data, e2@raw_data), e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Publication_set_operations
#' @export
setMethod(
  "-",
  signature = c("Publication", "Publication"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$pubmed_id)) {
          self_key_set <- unique(c(self_key_set, i$pubmed_id))
          self_dict[[i$pubmed_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$pubmed_id)) {
          other_key_set <- unique(c(other_key_set, j$pubmed_id))
        }
      }

      sub_key <- self_key_set[!self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Publication(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Publication_set_operations
#' @export
setMethod(
  "&",
  signature = c("Publication", "Publication"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$pubmed_id)) {
          self_key_set <- unique(c(self_key_set, i$pubmed_id))
          self_dict[[i$pubmed_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$pubmed_id)) {
          other_key_set <- unique(c(other_key_set, j$pubmed_id))
        }
      }

      sub_key <- self_key_set[self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Publication(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Publication_set_operations
#' @export
setMethod(
  "|",
  signature = c("Publication", "Publication"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      and_dict <- list()

      for (i in e1@raw_data) {
        if (!is.null(i$pubmed_id)) {
          and_dict[[i$pubmed_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$pubmed_id)) {
          and_dict[[j$pubmed_id]] <- j
        }
      }

      data <- unname(and_dict)
      return(Publication(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Publication_set_operations
#' @export
setMethod(
  "^",
  signature = c("Publication", "Publication"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      and_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$pubmed_id)) {
          self_key_set <- unique(c(self_key_set, i$pubmed_id))
          and_dict[[i$pubmed_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$pubmed_id)) {
          other_key_set <- unique(c(other_key_set, j$pubmed_id))
          and_dict[[j$pubmed_id]] <- j
        }
      }

      all_keys <- unique(c(self_key_set, other_key_set))
      common_keys <- self_key_set[self_key_set %in% other_key_set]
      sub_key <- all_keys[!all_keys %in% common_keys]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(and_dict[[k]]))
      }

      return(Publication(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' Compare Publication objects for equality
#'
#' Test whether two Publication objects are equal. Two objects are equal
#' when they run in the same mode and store identical records.
#'
#' @param e1 A Publication object
#' @param e2 A Publication object
#'
#' @return A logical value
#'
#' @export
setMethod(
  "==",
  signature = c("Publication", "Publication"),
  function(e1, e2) {
    if (is.null(e1) || is.null(e2)) {
      return(is.null(e1) && is.null(e2))
    }
    return(identical(e1@raw_data, e2@raw_data) && e1@mode == e2@mode)
  }
)

#' Summarize a Publication object
#'
#' Print a summary of a Publication object: the mode, the number of raw
#' records, the number of publications and, in Fat mode, the row counts
#' of both data.frames.
#'
#' @param object A Publication object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "summary",
  signature = "Publication",
  function(object) {
    cat("Publication Object Summary\n")
    cat("==========================\n")
    cat("Mode:", object@mode, "\n")
    cat("Number of raw records:", length(object@raw_data), "\n")
    cat("Number of publications:", nrow(object@publications), "\n")
    if (object@mode == "Fat") {
      cat("\nDataFrames:\n")
      cat("  publications:", nrow(object@publications), "rows\n")
      cat("  authors:", nrow(object@authors), "rows\n")
    }
  }
)
