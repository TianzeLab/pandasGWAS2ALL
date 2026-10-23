#' Get Trait data from the server
#'
#' This function retrieves trait data from the GWAS API based on various search
#' criteria.
#'
#' @param efo_trait The trait name or label
#' @param uri The trait URI or unique identifier
#' @param efo_id The trait URI shortform
#' @param pubmed_id pubmedId of the publication
#' @param mapped_gene Gene(s) overlapping the variant. If a variant is
#'   intergenic, the closest 5' and 3' genes are listed
#' @param extended_geneset Show extended matching genes in addition to the mapped
#'   genes
#' @param cached Whether or not to try to get data from the cache
#' @param mode Fat or Thin. Specifies the mode of the returned object
#'
#' @return A Trait object. Attributes of type data.frame have hierarchical
#'   dependencies
#'
#' @examples
#' \dontrun{
#' ch <- get_traits(efo_id = "EFO_0001360")
#' }
#'
#' @importFrom utils URLencode
#'
#' @export
#'
get_traits <- function(
    efo_trait = NULL,
    uri = NULL,
    efo_id = NULL,
    pubmed_id = NULL,
    mapped_gene = NULL,
    extended_geneset = NULL,
    cached = TRUE,
    mode = "Fat"
) {
  BASE_URL <- "https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits"

  params <- list(
    efo_trait = efo_trait,
    uri = uri,
    efo_id = efo_id,
    pubmed_id = pubmed_id,
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

  return(
    Trait(
      data = get_trait(query_url, cached = cached),
      mode = mode
    )
  )
}
