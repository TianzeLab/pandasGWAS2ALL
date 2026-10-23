#' Study class
#'
#' An S4 class that stores data of type Study fetched from the GWAS API.
#' In Fat mode it holds 8 data.frames with hierarchical dependencies; in
#' Thin mode it only stores the raw list data.
#'
#' @slot raw_data list. Converted from obtained JSON data
#' @slot mode character. Fat or Thin
#' @slot studies data.frame
#' @slot efo_traits data.frame
#' @slot genotyping_technologies data.frame
#' @slot discovery_ancestries data.frame
#' @slot replication_ancestries data.frame
#' @slot bg_efo_traits data.frame
#' @slot cohorts data.frame
#' @slot array_manufacturers data.frame
#'
#' @importFrom methods new
#'
#' @export
#'
setClass(
  "Study",
  slots = c(
    raw_data = "list",
    mode = "character",
    studies = "data.frame",
    efo_traits = "data.frame",
    genotyping_technologies = "data.frame",
    discovery_ancestries = "data.frame",
    replication_ancestries = "data.frame",
    bg_efo_traits = "data.frame",
    cohorts = "data.frame",
    array_manufacturers = "data.frame"
  )
)

#' Study constructor
#'
#' Creates a Study object from a list of raw records. To understand the
#' significance of each column of the data.frames, visit "Study" in the
#' GWAS Catalog documentation
#' (https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference).
#'
#' @param data A list of study records
#' @param mode The mode of the object: "Fat" or "Thin"
#'
#' @return A Study object
#'
#' @examples
#' \dontrun{
#' ch <- Study(
#'   list(list(accession_id = "GCST000854", study_id = "GCST000854")),
#'   mode = "Fat"
#' )
#' }
#'
#' @export
#'
Study <- function(data = NULL, mode = "Fat") {
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
      normalized_data <- c(normalized_data, list(rec))
    }
  }

  STUDY_COLUMNS <- c(
    "initial_sample_size",
    "replication_sample_size",
    "study_id",
    "gxe",
    "gxg",
    "snp_count",
    "qualifier",
    "imputed",
    "pooled",
    "accession_id",
    "user_requested",
    "pubmed_id",
    "platforms",
    "disease_trait",
    "terms_of_license",
    "full_summary_stats",
    "full_summary_stats_available"
  )

  EFO_TRAITS_COLUMNS <- c("accession_id", "efo_id", "efo_trait")
  GENOTYPING_TECHNOLOGIES_COLUMNS <- c("accession_id", "genotyping_technologies")
  DISCOVERY_ANCESTRIES_COLUMNS <- c("accession_id", "discovery_ancestry")
  REPLICATION_ANCESTRIES_COLUMNS <- c("accession_id", "replication_ancestry")
  BG_EFO_TRAITS_COLUMNS <- c("accession_id", "efo_id", "efo_trait")
  COHORTS_COLUMNS <- c("accession_id", "cohort")
  ARRAY_MANUFACTURERS_COLUMNS <- c("accession_id", "array_manufacturer")

  NESTED_FIELDS <- c(
    "efo_traits",
    "genotyping_technologies",
    "discovery_ancestry",
    "replication_ancestry",
    "bg_efo_traits",
    "cohort",
    "array_manufacturer"
  )

  obj <- new("Study")
  obj@mode <- mode

  obj@studies <- data.frame(matrix(ncol = length(STUDY_COLUMNS), nrow = 0))
  colnames(obj@studies) <- STUDY_COLUMNS

  obj@efo_traits <- data.frame(matrix(ncol = length(EFO_TRAITS_COLUMNS), nrow = 0))
  colnames(obj@efo_traits) <- EFO_TRAITS_COLUMNS

  obj@genotyping_technologies <- data.frame(matrix(ncol = length(GENOTYPING_TECHNOLOGIES_COLUMNS), nrow = 0))
  colnames(obj@genotyping_technologies) <- GENOTYPING_TECHNOLOGIES_COLUMNS

  obj@discovery_ancestries <- data.frame(matrix(ncol = length(DISCOVERY_ANCESTRIES_COLUMNS), nrow = 0))
  colnames(obj@discovery_ancestries) <- DISCOVERY_ANCESTRIES_COLUMNS

  obj@replication_ancestries <- data.frame(matrix(ncol = length(REPLICATION_ANCESTRIES_COLUMNS), nrow = 0))
  colnames(obj@replication_ancestries) <- REPLICATION_ANCESTRIES_COLUMNS

  obj@bg_efo_traits <- data.frame(matrix(ncol = length(BG_EFO_TRAITS_COLUMNS), nrow = 0))
  colnames(obj@bg_efo_traits) <- BG_EFO_TRAITS_COLUMNS

  obj@cohorts <- data.frame(matrix(ncol = length(COHORTS_COLUMNS), nrow = 0))
  colnames(obj@cohorts) <- COHORTS_COLUMNS

  obj@array_manufacturers <- data.frame(matrix(ncol = length(ARRAY_MANUFACTURERS_COLUMNS), nrow = 0))
  colnames(obj@array_manufacturers) <- ARRAY_MANUFACTURERS_COLUMNS

  if (mode == "Thin") {
    obj@raw_data <- normalized_data
    return(obj)
  }

  standardized_data <- list()
  for (item in normalized_data) {
    standardized <- list()
    for (col in STUDY_COLUMNS) {
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

  study_rows <- list()
  efo_traits_rows <- list()
  genotyping_technologies_rows <- list()
  discovery_ancestries_rows <- list()
  replication_ancestries_rows <- list()
  bg_efo_traits_rows <- list()
  cohorts_rows <- list()
  array_manufacturers_rows <- list()

  for (record in standardized_data) {
    accession_id <- record$accession_id

    study_rows <- c(study_rows, list(record[STUDY_COLUMNS]))

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
          accession_id = accession_id,
          efo_id = item$efo_id,
          efo_trait = item$efo_trait
        )))
      } else {
        efo_traits_rows <- c(efo_traits_rows, list(list(
          accession_id = accession_id,
          efo_id = NA,
          efo_trait = item
        )))
      }
    }

    raw_genotyping_technologies <- record$genotyping_technologies
    if (is.null(raw_genotyping_technologies)) {
      iterable_genotyping_technologies <- list()
    } else if (is.list(raw_genotyping_technologies) && is.null(names(raw_genotyping_technologies))) {
      iterable_genotyping_technologies <- raw_genotyping_technologies
    } else {
      iterable_genotyping_technologies <- list(raw_genotyping_technologies)
    }
    for (item in iterable_genotyping_technologies) {
      genotyping_technologies_rows <- c(genotyping_technologies_rows, list(list(
        accession_id = accession_id,
        genotyping_technologies = item
      )))
    }

    raw_discovery_ancestry <- record$discovery_ancestry
    if (is.null(raw_discovery_ancestry)) {
      iterable_discovery_ancestry <- list()
    } else if (is.list(raw_discovery_ancestry) && is.null(names(raw_discovery_ancestry))) {
      iterable_discovery_ancestry <- raw_discovery_ancestry
    } else {
      iterable_discovery_ancestry <- list(raw_discovery_ancestry)
    }
    for (item in iterable_discovery_ancestry) {
      discovery_ancestries_rows <- c(discovery_ancestries_rows, list(list(
        accession_id = accession_id,
        discovery_ancestry = item
      )))
    }

    raw_replication_ancestry <- record$replication_ancestry
    if (is.null(raw_replication_ancestry)) {
      iterable_replication_ancestry <- list()
    } else if (is.list(raw_replication_ancestry) && is.null(names(raw_replication_ancestry))) {
      iterable_replication_ancestry <- raw_replication_ancestry
    } else {
      iterable_replication_ancestry <- list(raw_replication_ancestry)
    }
    for (item in iterable_replication_ancestry) {
      replication_ancestries_rows <- c(replication_ancestries_rows, list(list(
        accession_id = accession_id,
        replication_ancestry = item
      )))
    }

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
          accession_id = accession_id,
          efo_id = item$efo_id,
          efo_trait = item$efo_trait
        )))
      } else {
        bg_efo_traits_rows <- c(bg_efo_traits_rows, list(list(
          accession_id = accession_id,
          efo_id = NA,
          efo_trait = item
        )))
      }
    }

    raw_cohort <- record$cohort
    if (is.null(raw_cohort)) {
      iterable_cohort <- list()
    } else if (is.list(raw_cohort) && is.null(names(raw_cohort))) {
      iterable_cohort <- raw_cohort
    } else {
      iterable_cohort <- list(raw_cohort)
    }
    for (item in iterable_cohort) {
      cohorts_rows <- c(cohorts_rows, list(list(
        accession_id = accession_id,
        cohort = item
      )))
    }

    raw_array_manufacturer <- record$array_manufacturer
    if (is.null(raw_array_manufacturer)) {
      iterable_array_manufacturer <- list()
    } else if (is.list(raw_array_manufacturer) && is.null(names(raw_array_manufacturer))) {
      iterable_array_manufacturer <- raw_array_manufacturer
    } else {
      iterable_array_manufacturer <- list(raw_array_manufacturer)
    }
    for (item in iterable_array_manufacturer) {
      array_manufacturers_rows <- c(array_manufacturers_rows, list(list(
        accession_id = accession_id,
        array_manufacturer = item
      )))
    }
  }

  if (length(study_rows) > 0) {
    built <- data.frame(row.names = seq_along(study_rows))
    for (col in STUDY_COLUMNS) {
      col_values <- vector("list", length(study_rows))
      for (r_idx in seq_along(study_rows)) {
        cell <- study_rows[[r_idx]][[col]]
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
    obj@studies <- built
  }

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

  if (length(genotyping_technologies_rows) > 0) {
    built <- data.frame(row.names = seq_along(genotyping_technologies_rows))
    for (col in GENOTYPING_TECHNOLOGIES_COLUMNS) {
      col_values <- vector("list", length(genotyping_technologies_rows))
      for (r_idx in seq_along(genotyping_technologies_rows)) {
        cell <- genotyping_technologies_rows[[r_idx]][[col]]
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
    obj@genotyping_technologies <- built
  }

  if (length(discovery_ancestries_rows) > 0) {
    built <- data.frame(row.names = seq_along(discovery_ancestries_rows))
    for (col in DISCOVERY_ANCESTRIES_COLUMNS) {
      col_values <- vector("list", length(discovery_ancestries_rows))
      for (r_idx in seq_along(discovery_ancestries_rows)) {
        cell <- discovery_ancestries_rows[[r_idx]][[col]]
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
    obj@discovery_ancestries <- built
  }

  if (length(replication_ancestries_rows) > 0) {
    built <- data.frame(row.names = seq_along(replication_ancestries_rows))
    for (col in REPLICATION_ANCESTRIES_COLUMNS) {
      col_values <- vector("list", length(replication_ancestries_rows))
      for (r_idx in seq_along(replication_ancestries_rows)) {
        cell <- replication_ancestries_rows[[r_idx]][[col]]
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
    obj@replication_ancestries <- built
  }

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

  if (length(cohorts_rows) > 0) {
    built <- data.frame(row.names = seq_along(cohorts_rows))
    for (col in COHORTS_COLUMNS) {
      col_values <- vector("list", length(cohorts_rows))
      for (r_idx in seq_along(cohorts_rows)) {
        cell <- cohorts_rows[[r_idx]][[col]]
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
    obj@cohorts <- built
  }

  if (length(array_manufacturers_rows) > 0) {
    built <- data.frame(row.names = seq_along(array_manufacturers_rows))
    for (col in ARRAY_MANUFACTURERS_COLUMNS) {
      col_values <- vector("list", length(array_manufacturers_rows))
      for (r_idx in seq_along(array_manufacturers_rows)) {
        cell <- array_manufacturers_rows[[r_idx]][[col]]
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
    obj@array_manufacturers <- built
  }

  return(obj)
}

#' Show a Study object
#'
#' Print a concise overview of a Study object: the current mode and the
#' size of each data component.
#'
#' @param object A Study object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "show",
  signature = "Study",
  function(object) {
    if (object@mode == "Fat") {
      cat(
        "Study is running in fat mode. It has 8 DataFrames with hierarchical dependencies.\n",
        "studies: ", nrow(object@studies), " rows\n",
        "|\n -efo_traits: ", nrow(object@efo_traits), " rows\n",
        "|\n -genotyping_technologies: ", nrow(object@genotyping_technologies), " rows\n",
        "|\n -discovery_ancestries: ", nrow(object@discovery_ancestries), " rows\n",
        "|\n -replication_ancestries: ", nrow(object@replication_ancestries), " rows\n",
        "|\n -bg_efo_traits: ", nrow(object@bg_efo_traits), " rows\n",
        "|\n -cohorts: ", nrow(object@cohorts), " rows\n",
        "|\n -array_manufacturers: ", nrow(object@array_manufacturers), " rows\n",
        sep = ""
      )
    } else if (object@mode == "Thin") {
      cat(
        "Study is running in thin mode. It has 1 list that contains the raw data.\n",
        "raw_data: a list of size ", length(object@raw_data), ".\n",
        sep = ""
      )
    }
  }
)

#' Length of a Study object
#'
#' Get the number of study records stored in the object.
#'
#' @param x A Study object
#'
#' @return An integer, the number of records in the raw_data slot
#'
#' @export
setMethod(
  "length",
  signature = "Study",
  function(x) {
    return(length(x@raw_data))
  }
)

#' Component names of a Study object
#'
#' Get the names of the components (slots) of a Study object.
#'
#' @param x A Study object
#'
#' @return A character vector of component names
#'
#' @export
setMethod(
  "names",
  signature = "Study",
  function(x) {
    return(c(
      "raw_data",
      "mode",
      "studies",
      "efo_traits",
      "genotyping_technologies",
      "discovery_ancestries",
      "replication_ancestries",
      "bg_efo_traits",
      "cohorts",
      "array_manufacturers"
    ))
  }
)

#' Subset a Study object
#'
#' Extract study records by position or by accession identifier and
#' return a new Study object of the same mode containing only the
#' selected records.
#'
#' @param x A Study object
#' @param i A numeric vector of positions (1-based) or a character vector of
#'   accession_ids. A list is unlisted before use. Leaving i missing
#'   returns x unchanged.
#'
#' @return A Study object containing the selected records
#'
#' @examples
#' \dontrun{
#' a <- get_studies(accession_id = "GCST000854")
#' a[1]
#' a["GCST000854"]
#' }
#'
#' @export
setMethod(
  "[",
  signature = "Study",
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
      if (!is.null(raw_data[[j]]$accession_id)) {
        raw_data_dict[[raw_data[[j]]$accession_id]] <- raw_data[[j]]
      }
    }

    sub_set <- list()
    for (i_val in arr) {
      if (is.character(i_val)) {
        if (!i_val %in% names(raw_data_dict)) {
          stop("Unknown accession_id: ", i_val)
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

    return(Study(sub_set, x@mode))
  }
)

#' Set operations on Study objects
#'
#' Combine or compare two Study objects of the same mode as sets of study
#' records, identified and deduplicated by accession_id.
#'
#' @param e1 A Study object
#' @param e2 A Study object
#'
#' @return A Study object of the same mode
#'
#' @examples
#' \dontrun{
#' a <- get_studies(accession_id = "GCST000854")
#' b <- get_studies(accession_id = "GCST90841394")
#' a + b
#' a - b
#' a & b
#' a | b
#' a ^ b
#' }
#'
#' @rdname Study_set_operations
#' @export
setMethod(
  "+",
  signature = c("Study", "Study"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      return(Study(c(e1@raw_data, e2@raw_data), e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Study_set_operations
#' @export
setMethod(
  "-",
  signature = c("Study", "Study"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$accession_id)) {
          self_key_set <- unique(c(self_key_set, i$accession_id))
          self_dict[[i$accession_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$accession_id)) {
          other_key_set <- unique(c(other_key_set, j$accession_id))
        }
      }

      sub_key <- self_key_set[!self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Study(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Study_set_operations
#' @export
setMethod(
  "&",
  signature = c("Study", "Study"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$accession_id)) {
          self_key_set <- unique(c(self_key_set, i$accession_id))
          self_dict[[i$accession_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$accession_id)) {
          other_key_set <- unique(c(other_key_set, j$accession_id))
        }
      }

      sub_key <- self_key_set[self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Study(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Study_set_operations
#' @export
setMethod(
  "|",
  signature = c("Study", "Study"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      and_dict <- list()

      for (i in e1@raw_data) {
        if (!is.null(i$accession_id)) {
          and_dict[[i$accession_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$accession_id)) {
          and_dict[[j$accession_id]] <- j
        }
      }

      data <- unname(and_dict)
      return(Study(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Study_set_operations
#' @export
setMethod(
  "^",
  signature = c("Study", "Study"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      and_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$accession_id)) {
          self_key_set <- unique(c(self_key_set, i$accession_id))
          and_dict[[i$accession_id]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$accession_id)) {
          other_key_set <- unique(c(other_key_set, j$accession_id))
          and_dict[[j$accession_id]] <- j
        }
      }

      all_keys <- unique(c(self_key_set, other_key_set))
      common_keys <- self_key_set[self_key_set %in% other_key_set]
      sub_key <- all_keys[!all_keys %in% common_keys]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(and_dict[[k]]))
      }

      return(Study(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' Compare Study objects for equality
#'
#' Test whether two Study objects are equal. Two objects are equal when
#' they run in the same mode and store identical records.
#'
#' @param e1 A Study object
#' @param e2 A Study object
#'
#' @return A logical value
#'
#' @export
setMethod(
  "==",
  signature = c("Study", "Study"),
  function(e1, e2) {
    if (is.null(e1) || is.null(e2)) {
      return(is.null(e1) && is.null(e2))
    }
    return(identical(e1@raw_data, e2@raw_data) && e1@mode == e2@mode)
  }
)

#' Summarize a Study object
#'
#' Print a summary of a Study object: the mode, the number of raw records,
#' the number of studies and, in Fat mode, the row counts of all 8
#' data.frames.
#'
#' @param object A Study object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "summary",
  signature = "Study",
  function(object) {
    cat("Study Object Summary\n")
    cat("====================\n")
    cat("Mode:", object@mode, "\n")
    cat("Number of raw records:", length(object@raw_data), "\n")
    cat("Number of studies:", nrow(object@studies), "\n")
    if (object@mode == "Fat") {
      cat("\nDataFrames:\n")
      cat("  studies:", nrow(object@studies), "rows\n")
      cat("  efo_traits:", nrow(object@efo_traits), "rows\n")
      cat("  genotyping_technologies:", nrow(object@genotyping_technologies), "rows\n")
      cat("  discovery_ancestries:", nrow(object@discovery_ancestries), "rows\n")
      cat("  replication_ancestries:", nrow(object@replication_ancestries), "rows\n")
      cat("  bg_efo_traits:", nrow(object@bg_efo_traits), "rows\n")
      cat("  cohorts:", nrow(object@cohorts), "rows\n")
      cat("  array_manufacturers:", nrow(object@array_manufacturers), "rows\n")
    }
  }
)
