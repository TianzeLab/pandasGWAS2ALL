#' Get Association data from the server
#'
#' This function retrieves association data from the GWAS API based on various
#' search criteria.
#'
#' @param association_id Primary identifier of association table
#' @param pubmed_id pubmedId of the publication
#' @param rs_id SNP rsID; if a haplotype it may include more than one rs number
#'   (multiple SNPs comprising the haplotype)
#' @param full_pvalue_set Whether full summary statistics are available for this
#'   study
#' @param accession_id The study's GWAS Catalog accession ID
#' @param efo_trait The trait name or label
#' @param efo_id The trait URI shortform
#' @param show_child_trait display entities for descendants of a parent Efo Trait
#'   Term
#' @param mapped_gene Gene(s) overlapping the variant. If a variant is
#'   intergenic, the closest 5' and 3' genes are listed
#' @param extended_geneset Show extended matching genes in addition to the mapped
#'   genes
#' @param cached Whether or not to try to get data from the cache
#' @param mode Fat or Thin. Specifies the mode of the returned object
#'
#' @return An Association object. Attributes of type data.frame have hierarchical
#'   dependencies
#'
#' @examples
#' \dontrun{
#' ch <- get_associations(association_id = "169106883")
#' }
#'
#' @importFrom utils URLencode
#'
#' @export
#'
get_associations <- function(
    association_id = NULL,
    pubmed_id = NULL,
    rs_id = NULL,
    full_pvalue_set = NULL,
    accession_id = NULL,
    efo_trait = NULL,
    efo_id = NULL,
    show_child_trait = NULL,
    mapped_gene = NULL,
    extended_geneset = NULL,
    cached = TRUE,
    mode = "Fat"
) {
  BASE_URL <- "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"

  params <- list(
    pubmed_id = pubmed_id,
    rs_id = rs_id,
    full_pvalue_set = full_pvalue_set,
    accession_id = accession_id,
    efo_trait = efo_trait,
    efo_id = efo_id,
    show_child_trait = show_child_trait,
    mapped_gene = mapped_gene,
    extended_geneset = extended_geneset
  )
  params <- params[!vapply(params, is.null, logical(1))]

  if (length(params) > 0) {
    encoded_pieces <- character(0)
    for (param_name in names(params)) {
      raw_value <- params[[param_name]]
      if (is.logical(raw_value)) {
        value_str <- if (raw_value) "True" else "False"
      } else {
        value_str <- as.character(raw_value)
      }
      encoded_value <- gsub(
        "%20", "+",
        URLencode(enc2utf8(value_str), reserved = TRUE),
        fixed = TRUE
      )
      encoded_pieces <- c(encoded_pieces, paste0(param_name, "=", encoded_value))
    }
    query_url <- paste0(BASE_URL, "?", paste(encoded_pieces, collapse = "&"))
  } else {
    query_url <- BASE_URL
  }

  has_id <- !is.null(association_id)
  has_query_filters <- length(params) > 0

  result_by_id <- NULL
  result_by_query <- NULL

  if (has_id) {
    id_url <- paste0(
      BASE_URL, "/",
      URLencode(enc2utf8(as.character(association_id)), reserved = TRUE)
    )
    result_by_id <- Association(
      data = get_association(id_url, cached = cached),
      mode = mode
    )
  }

  if (!has_id || has_query_filters) {
    result_by_query <- Association(
      data = get_association(query_url, cached = cached),
      mode = mode
    )
  }

  if (has_id && has_query_filters) {
    return(result_by_id & result_by_query)
  }
  if (has_id) {
    return(result_by_id)
  }
  return(result_by_query)
}
