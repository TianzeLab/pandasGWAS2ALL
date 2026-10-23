#' Gene class
#'
#' An S4 class that stores data of type Gene fetched from the GWAS API.
#' In Fat mode it holds 3 data.frames with hierarchical dependencies; in
#' Thin mode it only stores the raw list data.
#'
#' @slot raw_data list. Converted from obtained JSON data
#' @slot mode character. Fat or Thin
#' @slot genes data.frame
#' @slot ensembl_gene_ids data.frame
#' @slot entrez_gene_ids data.frame
#'
#' @importFrom methods new
#'
#' @export
#'
setClass(
  "Gene",
  slots = c(
    raw_data = "list",
    mode = "character",
    genes = "data.frame",
    ensembl_gene_ids = "data.frame",
    entrez_gene_ids = "data.frame"
  )
)

#' Gene constructor
#'
#' Creates a Gene object from a list of raw records. To understand the
#' significance of each column of the data.frames, visit "Gene" in the
#' GWAS Catalog documentation
#' (https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference).
#'
#' @param data A list of gene records
#' @param mode The mode of the object: "Fat" or "Thin"
#'
#' @return A Gene object
#'
#' @examples
#' \dontrun{
#' ch <- Gene(
#'   list(list(gene_name = "HBS1L", biotype = "protein_coding")),
#'   mode = "Fat"
#' )
#' }
#'
#' @export
#'
Gene <- function(data = NULL, mode = "Fat") {
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
      if (!is.null(rec$gene_name)) {
        rec$gene_name <- as.character(rec$gene_name)
      }
      normalized_data <- c(normalized_data, list(rec))
    }
  }

  GENE_COLUMNS <- c(
    "gene_name",
    "gene_description",
    "location",
    "cytogenic_region",
    "biotype"
  )

  ENSEMBL_GENE_IDS_COLUMNS <- c("gene_name", "ensembl_gene_ids")
  ENTREZ_GENE_IDS_COLUMNS <- c("gene_name", "entrez_gene_ids")

  NESTED_FIELDS <- c(
    "ensembl_gene_ids",
    "entrez_gene_ids"
  )

  obj <- new("Gene")
  obj@mode <- mode

  obj@genes <- data.frame(matrix(ncol = length(GENE_COLUMNS), nrow = 0))
  colnames(obj@genes) <- GENE_COLUMNS

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
    standardized <- list()
    for (col in GENE_COLUMNS) {
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

  gene_rows <- list()
  ensembl_gene_ids_rows <- list()
  entrez_gene_ids_rows <- list()

  for (record in standardized_data) {
    gene_name <- record$gene_name

    gene_rows <- c(gene_rows, list(record[GENE_COLUMNS]))

    raw_ensembl_gene_ids <- record$ensembl_gene_ids
    if (is.null(raw_ensembl_gene_ids)) {
      iterable_ensembl_gene_ids <- list()
    } else if (is.list(raw_ensembl_gene_ids) && is.null(names(raw_ensembl_gene_ids))) {
      iterable_ensembl_gene_ids <- raw_ensembl_gene_ids
    } else {
      iterable_ensembl_gene_ids <- list(raw_ensembl_gene_ids)
    }
    for (item in iterable_ensembl_gene_ids) {
      ensembl_gene_ids_rows <- c(ensembl_gene_ids_rows, list(list(
        gene_name = gene_name,
        ensembl_gene_ids = item
      )))
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
      entrez_gene_ids_rows <- c(entrez_gene_ids_rows, list(list(
        gene_name = gene_name,
        entrez_gene_ids = item
      )))
    }
  }

  if (length(gene_rows) > 0) {
    built <- data.frame(row.names = seq_along(gene_rows))
    for (col in GENE_COLUMNS) {
      col_values <- vector("list", length(gene_rows))
      for (r_idx in seq_along(gene_rows)) {
        cell <- gene_rows[[r_idx]][[col]]
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
    obj@genes <- built
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

#' Show a Gene object
#'
#' Print a concise overview of a Gene object: the current mode and the
#' size of each data component.
#'
#' @param object A Gene object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "show",
  signature = "Gene",
  function(object) {
    if (object@mode == "Fat") {
      cat(
        "Gene is running in fat mode. It has 3 DataFrames with hierarchical dependencies.\n",
        "genes: ", nrow(object@genes), " rows\n",
        "|\n -ensembl_gene_ids: ", nrow(object@ensembl_gene_ids), " rows\n",
        "|\n -entrez_gene_ids: ", nrow(object@entrez_gene_ids), " rows\n",
        sep = ""
      )
    } else if (object@mode == "Thin") {
      cat(
        "Gene is running in thin mode. It has 1 list that contains the raw data.\n",
        "raw_data: a list of size ", length(object@raw_data), ".\n",
        sep = ""
      )
    }
  }
)

#' Length of a Gene object
#'
#' Get the number of gene records stored in the object.
#'
#' @param x A Gene object
#'
#' @return An integer, the number of records in the raw_data slot
#'
#' @export
setMethod(
  "length",
  signature = "Gene",
  function(x) {
    return(length(x@raw_data))
  }
)

#' Component names of a Gene object
#'
#' Get the names of the components (slots) of a Gene object.
#'
#' @param x A Gene object
#'
#' @return A character vector of component names
#'
#' @export
setMethod(
  "names",
  signature = "Gene",
  function(x) {
    return(c("raw_data", "mode", "genes", "ensembl_gene_ids", "entrez_gene_ids"))
  }
)

#' Subset a Gene object
#'
#' Extract gene records by position or by gene name and return a new Gene
#' object of the same mode containing only the selected records.
#'
#' @param x A Gene object
#' @param i A numeric vector of positions (1-based) or a character vector of
#'   gene_names. A list is unlisted before use. Leaving i missing returns x
#'   unchanged.
#'
#' @return A Gene object containing the selected records
#'
#' @examples
#' \dontrun{
#' a <- get_genes(gene_name = "HBS1L")
#' a[1]
#' a["HBS1L"]
#' }
#'
#' @export
setMethod(
  "[",
  signature = "Gene",
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
      if (!is.null(raw_data[[j]]$gene_name)) {
        raw_data_dict[[raw_data[[j]]$gene_name]] <- raw_data[[j]]
      }
    }

    sub_set <- list()
    for (i_val in arr) {
      if (is.character(i_val)) {
        if (!i_val %in% names(raw_data_dict)) {
          stop("Unknown gene_name: ", i_val)
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

    return(Gene(sub_set, x@mode))
  }
)

#' Set operations on Gene objects
#'
#' Combine or compare two Gene objects of the same mode as sets of gene
#' records, identified and deduplicated by gene_name.
#'
#' @param e1 A Gene object
#' @param e2 A Gene object
#'
#' @return A Gene object of the same mode
#'
#' @examples
#' \dontrun{
#' a <- get_genes(gene_name = "HBS1L")
#' b <- get_genes(gene_name = "SAMD5")
#' a + b
#' a - b
#' a & b
#' a | b
#' a ^ b
#' }
#'
#' @rdname Gene_set_operations
#' @export
setMethod(
  "+",
  signature = c("Gene", "Gene"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      return(Gene(c(e1@raw_data, e2@raw_data), e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Gene_set_operations
#' @export
setMethod(
  "-",
  signature = c("Gene", "Gene"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$gene_name)) {
          self_key_set <- unique(c(self_key_set, i$gene_name))
          self_dict[[i$gene_name]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$gene_name)) {
          other_key_set <- unique(c(other_key_set, j$gene_name))
        }
      }

      sub_key <- self_key_set[!self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Gene(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Gene_set_operations
#' @export
setMethod(
  "&",
  signature = c("Gene", "Gene"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$gene_name)) {
          self_key_set <- unique(c(self_key_set, i$gene_name))
          self_dict[[i$gene_name]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$gene_name)) {
          other_key_set <- unique(c(other_key_set, j$gene_name))
        }
      }

      sub_key <- self_key_set[self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(Gene(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Gene_set_operations
#' @export
setMethod(
  "|",
  signature = c("Gene", "Gene"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      and_dict <- list()

      for (i in e1@raw_data) {
        if (!is.null(i$gene_name)) {
          and_dict[[i$gene_name]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$gene_name)) {
          and_dict[[j$gene_name]] <- j
        }
      }

      data <- unname(and_dict)
      return(Gene(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname Gene_set_operations
#' @export
setMethod(
  "^",
  signature = c("Gene", "Gene"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      and_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$gene_name)) {
          self_key_set <- unique(c(self_key_set, i$gene_name))
          and_dict[[i$gene_name]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$gene_name)) {
          other_key_set <- unique(c(other_key_set, j$gene_name))
          and_dict[[j$gene_name]] <- j
        }
      }

      all_keys <- unique(c(self_key_set, other_key_set))
      common_keys <- self_key_set[self_key_set %in% other_key_set]
      sub_key <- all_keys[!all_keys %in% common_keys]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(and_dict[[k]]))
      }

      return(Gene(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' Compare Gene objects for equality
#'
#' Test whether two Gene objects are equal. Two objects are equal when
#' they run in the same mode and store identical records.
#'
#' @param e1 A Gene object
#' @param e2 A Gene object
#'
#' @return A logical value
#'
#' @export
setMethod(
  "==",
  signature = c("Gene", "Gene"),
  function(e1, e2) {
    if (is.null(e1) || is.null(e2)) {
      return(is.null(e1) && is.null(e2))
    }
    return(identical(e1@raw_data, e2@raw_data) && e1@mode == e2@mode)
  }
)

#' Summarize a Gene object
#'
#' Print a summary of a Gene object: the mode, the number of raw records,
#' the number of genes and, in Fat mode, the row counts of all 3
#' data.frames.
#'
#' @param object A Gene object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "summary",
  signature = "Gene",
  function(object) {
    cat("Gene Object Summary\n")
    cat("===================\n")
    cat("Mode:", object@mode, "\n")
    cat("Number of raw records:", length(object@raw_data), "\n")
    cat("Number of genes:", nrow(object@genes), "\n")
    if (object@mode == "Fat") {
      cat("\nDataFrames:\n")
      cat("  genes:", nrow(object@genes), "rows\n")
      cat("  ensembl_gene_ids:", nrow(object@ensembl_gene_ids), "rows\n")
      cat("  entrez_gene_ids:", nrow(object@entrez_gene_ids), "rows\n")
    }
  }
)
