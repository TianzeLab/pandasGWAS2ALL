#' GenomicContext class
#'
#' An S4 class that stores data of type GenomicContext fetched from the
#' GWAS API. In Fat mode it holds 5 data.frames with hierarchical
#' dependencies; in Thin mode it only stores the raw list data.
#'
#' @slot raw_data list. Converted from obtained JSON data
#' @slot mode character. Fat or Thin
#' @slot genomic_contexts data.frame
#' @slot genes data.frame
#' @slot locations data.frame
#' @slot ensembl_gene_ids data.frame
#' @slot entrez_gene_ids data.frame
#'
#' @importFrom methods new
#'
#' @export
#'
setClass(
  "GenomicContext",
  slots = c(
    raw_data = "list",
    mode = "character",
    genomic_contexts = "data.frame",
    genes = "data.frame",
    locations = "data.frame",
    ensembl_gene_ids = "data.frame",
    entrez_gene_ids = "data.frame"
  )
)

#' GenomicContext constructor
#'
#' Creates a GenomicContext object from a list of raw records. To
#' understand the significance of each column of the data.frames, visit
#' "GenomicContext" in the GWAS Catalog documentation
#' (https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference).
#'
#' @param data A list of genomic context records
#' @param mode The mode of the object: "Fat" or "Thin"
#'
#' @return A GenomicContext object
#'
#' @examples
#' \dontrun{
#' ch <- GenomicContext(
#'   list(list(rs_id = "rs3093017", genomicContext_id = "206217817")),
#'   mode = "Fat"
#' )
#' }
#'
#' @export
#'
GenomicContext <- function(data = NULL, mode = "Fat") {
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
        list("rs_id", c("rs_id", "rsId")),
        list("genomicContext_id", c("genomicContext_id", "genomicContextId", "genomic_context_id")),
        list("is_intergenic", c("is_intergenic", "isIntergenic")),
        list("is_upstream", c("is_upstream", "isUpstream")),
        list("is_downstream", c("is_downstream", "isDownstream")),
        list("mapping_method", c("mapping_method", "mappingMethod")),
        list("is_closest_gene", c("is_closest_gene", "isClosestGene")),
        list("genes", c("genes", "gene")),
        list("locations", c("locations", "location")),
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
      if (!is.null(rec$rs_id)) {
        rec$rs_id <- as.character(rec$rs_id)
      }
      if (!is.null(rec$genomicContext_id)) {
        rec$genomicContext_id <- as.character(rec$genomicContext_id)
      }
      normalized_data <- c(normalized_data, list(rec))
    }
  }

  GENOMIC_CONTEXT_COLUMNS <- c(
    "rs_id",
    "genomicContext_id",
    "is_intergenic",
    "is_upstream",
    "is_downstream",
    "distance",
    "source",
    "mapping_method",
    "is_closest_gene"
  )

  GENES_COLUMNS <- c(
    "rs_id",
    "genomicContext_id",
    "gene_name",
    "gene_description",
    "location",
    "cytogenic_region",
    "biotype"
  )

  LOCATIONS_COLUMNS <- c("rs_id", "genomicContext_id", "chromosome_name", "chromosome_position", "region.name")
  ENSEMBL_GENE_IDS_COLUMNS <- c("rs_id", "genomicContext_id", "gene_name", "ensembl_gene_ids")
  ENTREZ_GENE_IDS_COLUMNS <- c("rs_id", "genomicContext_id", "gene_name", "entrez_gene_ids")

  NESTED_FIELDS <- c(
    "genes",
    "locations",
    "ensembl_gene_ids",
    "entrez_gene_ids"
  )

  obj <- new("GenomicContext")
  obj@mode <- mode

  obj@genomic_contexts <- data.frame(matrix(ncol = length(GENOMIC_CONTEXT_COLUMNS), nrow = 0))
  colnames(obj@genomic_contexts) <- GENOMIC_CONTEXT_COLUMNS

  obj@genes <- data.frame(matrix(ncol = length(GENES_COLUMNS), nrow = 0))
  colnames(obj@genes) <- GENES_COLUMNS

  obj@locations <- data.frame(matrix(ncol = length(LOCATIONS_COLUMNS), nrow = 0))
  colnames(obj@locations) <- LOCATIONS_COLUMNS

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
      list("rs_id", c("rs_id", "rsId")),
      list("genomicContext_id", c("genomicContext_id", "genomicContextId", "genomic_context_id")),
      list("is_intergenic", c("is_intergenic", "isIntergenic")),
      list("is_upstream", c("is_upstream", "isUpstream")),
      list("is_downstream", c("is_downstream", "isDownstream")),
      list("mapping_method", c("mapping_method", "mappingMethod")),
      list("is_closest_gene", c("is_closest_gene", "isClosestGene")),
      list("genes", c("genes", "gene")),
      list("locations", c("locations", "location")),
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
    for (col in GENOMIC_CONTEXT_COLUMNS) {
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

  genomic_context_rows <- list()
  genes_rows <- list()
  locations_rows <- list()
  ensembl_gene_ids_rows <- list()
  entrez_gene_ids_rows <- list()

  for (record in standardized_data) {
    rs_id <- record$rs_id
    genomicContext_id <- record$genomicContext_id

    genomic_context_rows <- c(genomic_context_rows, list(record[GENOMIC_CONTEXT_COLUMNS]))

    raw_genes <- record$genes
    if (is.null(raw_genes)) {
      iterable_genes <- list()
    } else if (is.list(raw_genes) && is.null(names(raw_genes))) {
      iterable_genes <- raw_genes
    } else {
      iterable_genes <- list(raw_genes)
    }
    for (item in iterable_genes) {
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
        genes_rows <- c(genes_rows, list(list(
          rs_id = rs_id,
          genomicContext_id = genomicContext_id,
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
            rs_id = rs_id,
            genomicContext_id = genomicContext_id,
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
            rs_id = rs_id,
            genomicContext_id = genomicContext_id,
            gene_name = gene_name,
            entrez_gene_ids = entrez_gene_id
          )))
        }
      } else {
        genes_rows <- c(genes_rows, list(list(
          rs_id = rs_id,
          genomicContext_id = genomicContext_id,
          gene_name = item,
          gene_description = NA,
          location = NA,
          cytogenic_region = NA,
          biotype = NA
        )))
      }
    }

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
        chromosome_name <- item$chromosome_name
        if (is.null(chromosome_name)) {
          chromosome_name <- item$chromosomeName
        }
        chromosome_position <- item$chromosome_position
        if (is.null(chromosome_position)) {
          chromosome_position <- item$chromosomePosition
        }
        region_name <- NULL
        raw_region <- item$region
        if (is.list(raw_region) && !is.null(names(raw_region))) {
          region_name <- raw_region$name
        } else if (!is.null(raw_region)) {
          region_name <- raw_region
        }
        locations_rows <- c(locations_rows, list(list(
          rs_id = rs_id,
          genomicContext_id = genomicContext_id,
          chromosome_name = chromosome_name,
          chromosome_position = chromosome_position,
          region.name = region_name
        )))
      } else {
        locations_rows <- c(locations_rows, list(list(
          rs_id = rs_id,
          genomicContext_id = genomicContext_id,
          chromosome_name = NA,
          chromosome_position = item,
          region.name = NA
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
          rs_id = rs_id,
          genomicContext_id = genomicContext_id,
          gene_name = gene_name,
          ensembl_gene_ids = ensembl_gene_ids
        )))
      } else {
        ensembl_gene_ids_rows <- c(ensembl_gene_ids_rows, list(list(
          rs_id = rs_id,
          genomicContext_id = genomicContext_id,
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
          rs_id = rs_id,
          genomicContext_id = genomicContext_id,
          gene_name = gene_name,
          entrez_gene_ids = entrez_gene_ids
        )))
      } else {
        entrez_gene_ids_rows <- c(entrez_gene_ids_rows, list(list(
          rs_id = rs_id,
          genomicContext_id = genomicContext_id,
          gene_name = NA,
          entrez_gene_ids = item
        )))
      }
    }
  }

  if (length(genomic_context_rows) > 0) {
    built <- data.frame(row.names = seq_along(genomic_context_rows))
    for (col in GENOMIC_CONTEXT_COLUMNS) {
      col_values <- vector("list", length(genomic_context_rows))
      for (r_idx in seq_along(genomic_context_rows)) {
        cell <- genomic_context_rows[[r_idx]][[col]]
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
    obj@genomic_contexts <- built
  }

  if (length(genes_rows) > 0) {
    built <- data.frame(row.names = seq_along(genes_rows))
    for (col in GENES_COLUMNS) {
      col_values <- vector("list", length(genes_rows))
      for (r_idx in seq_along(genes_rows)) {
        cell <- genes_rows[[r_idx]][[col]]
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

#' Show a GenomicContext object
#'
#' Print a concise overview of a GenomicContext object: the current mode
#' and the size of each data component.
#'
#' @param object A GenomicContext object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "show",
  signature = "GenomicContext",
  function(object) {
    if (object@mode == "Fat") {
      cat(
        "GenomicContext is running in fat mode. It has 5 DataFrames with hierarchical dependencies.\n",
        "genomic_contexts: ", nrow(object@genomic_contexts), " rows\n",
        "|\n -genes: ", nrow(object@genes), " rows\n",
        "|\n -locations: ", nrow(object@locations), " rows\n",
        "|\n -ensembl_gene_ids: ", nrow(object@ensembl_gene_ids), " rows\n",
        "|\n -entrez_gene_ids: ", nrow(object@entrez_gene_ids), " rows\n",
        sep = ""
      )
    } else if (object@mode == "Thin") {
      cat(
        "GenomicContext is running in thin mode. It has 1 list that contains the raw data.\n",
        "raw_data: a list of size ", length(object@raw_data), ".\n",
        sep = ""
      )
    }
  }
)

#' Length of a GenomicContext object
#'
#' Get the number of genomic context records stored in the object.
#'
#' @param x A GenomicContext object
#'
#' @return An integer, the number of records in the raw_data slot
#'
#' @export
setMethod(
  "length",
  signature = "GenomicContext",
  function(x) {
    return(length(x@raw_data))
  }
)

#' Component names of a GenomicContext object
#'
#' Get the names of the components (slots) of a GenomicContext object.
#'
#' @param x A GenomicContext object
#'
#' @return A character vector of component names
#'
#' @export
setMethod(
  "names",
  signature = "GenomicContext",
  function(x) {
    return(c(
      "raw_data",
      "mode",
      "genomic_contexts",
      "genes",
      "locations",
      "ensembl_gene_ids",
      "entrez_gene_ids"
    ))
  }
)

#' Subset a GenomicContext object
#'
#' Extract genomic context records by position or by the composite
#' identifier "rs_id|genomicContext_id" and return a new GenomicContext
#' object of the same mode containing only the selected records.
#'
#' @param x A GenomicContext object
#' @param i A numeric vector of positions (1-based) or a character vector of
#'   composite identifiers in the form "rs_id|genomicContext_id". A list
#'   is unlisted before use. Leaving i missing returns x unchanged.
#'
#' @return A GenomicContext object containing the selected records
#'
#' @examples
#' \dontrun{
#' a <- get_genomic_contexts(rs_id = "rs3093017", genomicContext_id = "206217817")
#' a[1]
#' a["rs3093017|206217817"]
#' }
#'
#' @export
setMethod(
  "[",
  signature = "GenomicContext",
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
      if (!is.null(raw_data[[j]]$rs_id) &&
          !is.null(raw_data[[j]]$genomicContext_id)) {
        key <- paste0(raw_data[[j]]$rs_id, "|", raw_data[[j]]$genomicContext_id)
        raw_data_dict[[key]] <- raw_data[[j]]
      }
    }

    sub_set <- list()
    for (i_val in arr) {
      if (is.character(i_val)) {
        if (!i_val %in% names(raw_data_dict)) {
          stop("Unknown rs_id|genomicContext_id: ", i_val)
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

    return(GenomicContext(sub_set, x@mode))
  }
)

#' Set operations on GenomicContext objects
#'
#' Combine or compare two GenomicContext objects of the same mode as sets
#' of genomic context records, identified and deduplicated by the
#' composite key "rs_id|genomicContext_id".
#'
#' @param e1 A GenomicContext object
#' @param e2 A GenomicContext object
#'
#' @return A GenomicContext object of the same mode
#'
#' @examples
#' \dontrun{
#' a <- get_genomic_contexts(rs_id = "rs3093017", genomicContext_id = "206217817")
#' b <- get_genomic_contexts(rs_id = "rs3093017", genomicContext_id = "206217818")
#' a + b
#' a - b
#' a & b
#' a | b
#' a ^ b
#' }
#'
#' @rdname GenomicContext_set_operations
#' @export
setMethod(
  "+",
  signature = c("GenomicContext", "GenomicContext"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      return(GenomicContext(c(e1@raw_data, e2@raw_data), e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname GenomicContext_set_operations
#' @export
setMethod(
  "-",
  signature = c("GenomicContext", "GenomicContext"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$rs_id) && !is.null(i$genomicContext_id)) {
          key <- paste0(i$rs_id, "|", i$genomicContext_id)
          self_key_set <- unique(c(self_key_set, key))
          self_dict[[key]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$rs_id) && !is.null(j$genomicContext_id)) {
          other_key_set <- unique(c(other_key_set, paste0(j$rs_id, "|", j$genomicContext_id)))
        }
      }

      sub_key <- self_key_set[!self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(GenomicContext(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname GenomicContext_set_operations
#' @export
setMethod(
  "&",
  signature = c("GenomicContext", "GenomicContext"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      self_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$rs_id) && !is.null(i$genomicContext_id)) {
          key <- paste0(i$rs_id, "|", i$genomicContext_id)
          self_key_set <- unique(c(self_key_set, key))
          self_dict[[key]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$rs_id) && !is.null(j$genomicContext_id)) {
          other_key_set <- unique(c(other_key_set, paste0(j$rs_id, "|", j$genomicContext_id)))
        }
      }

      sub_key <- self_key_set[self_key_set %in% other_key_set]
      data <- list()
      for (k in sub_key) {
        data <- c(data, list(self_dict[[k]]))
      }

      return(GenomicContext(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname GenomicContext_set_operations
#' @export
setMethod(
  "|",
  signature = c("GenomicContext", "GenomicContext"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      and_dict <- list()

      for (i in e1@raw_data) {
        if (!is.null(i$rs_id) && !is.null(i$genomicContext_id)) {
          and_dict[[paste0(i$rs_id, "|", i$genomicContext_id)]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$rs_id) && !is.null(j$genomicContext_id)) {
          and_dict[[paste0(j$rs_id, "|", j$genomicContext_id)]] <- j
        }
      }

      data <- unname(and_dict)
      return(GenomicContext(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' @rdname GenomicContext_set_operations
#' @export
setMethod(
  "^",
  signature = c("GenomicContext", "GenomicContext"),
  function(e1, e2) {
    if (e1@mode == e2@mode) {
      self_key_set <- character()
      and_dict <- list()
      other_key_set <- character()

      for (i in e1@raw_data) {
        if (!is.null(i$rs_id) && !is.null(i$genomicContext_id)) {
          key <- paste0(i$rs_id, "|", i$genomicContext_id)
          self_key_set <- unique(c(self_key_set, key))
          and_dict[[key]] <- i
        }
      }

      for (j in e2@raw_data) {
        if (!is.null(j$rs_id) && !is.null(j$genomicContext_id)) {
          key <- paste0(j$rs_id, "|", j$genomicContext_id)
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

      return(GenomicContext(data, e1@mode))
    } else {
      stop("Please input the same mode")
    }
  }
)

#' Compare GenomicContext objects for equality
#'
#' Test whether two GenomicContext objects are equal. Two objects are
#' equal when they run in the same mode and store identical records.
#'
#' @param e1 A GenomicContext object
#' @param e2 A GenomicContext object
#'
#' @return A logical value
#'
#' @export
setMethod(
  "==",
  signature = c("GenomicContext", "GenomicContext"),
  function(e1, e2) {
    if (is.null(e1) || is.null(e2)) {
      return(is.null(e1) && is.null(e2))
    }
    return(identical(e1@raw_data, e2@raw_data) && e1@mode == e2@mode)
  }
)

#' Summarize a GenomicContext object
#'
#' Print a summary of a GenomicContext object: the mode, the number of raw
#' records, the number of genomic contexts and, in Fat mode, the row
#' counts of all 5 data.frames.
#'
#' @param object A GenomicContext object
#'
#' @return None (called for its side effect of printing to the console)
#'
#' @export
setMethod(
  "summary",
  signature = "GenomicContext",
  function(object) {
    cat("GenomicContext Object Summary\n")
    cat("=============================\n")
    cat("Mode:", object@mode, "\n")
    cat("Number of raw records:", length(object@raw_data), "\n")
    cat("Number of genomic_contexts:", nrow(object@genomic_contexts), "\n")
    if (object@mode == "Fat") {
      cat("\nDataFrames:\n")
      cat("  genomic_contexts:", nrow(object@genomic_contexts), "rows\n")
      cat("  genes:", nrow(object@genes), "rows\n")
      cat("  locations:", nrow(object@locations), "rows\n")
      cat("  ensembl_gene_ids:", nrow(object@ensembl_gene_ids), "rows\n")
      cat("  entrez_gene_ids:", nrow(object@entrez_gene_ids), "rows\n")
    }
  }
)
