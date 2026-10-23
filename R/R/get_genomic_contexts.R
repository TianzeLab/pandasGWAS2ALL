#' Get GenomicContext data from the server
#'
#' This function retrieves genomic context data from the GWAS API for a given
#' SNP rs identifier and optional genomic context identifier.
#'
#' @param rs_id SNP rsID; if a haplotype it may include more than one rs number
#'   (multiple SNPs comprising the haplotype)
#' @param genomicContext_id primary identifier of genomicContext table
#' @param cached Whether or not to try to get data from the cache
#' @param mode Fat or Thin. Specifies the mode of the returned object
#'
#' @return A GenomicContext object. Attributes of type data.frame have
#'   hierarchical dependencies
#'
#' @examples
#' \dontrun{
#' ch <- get_genomic_contexts(rs_id = "rs3093017")
#' }
#'
#' @importFrom utils URLencode
#'
#' @export
#'
get_genomic_contexts <- function(
    rs_id = NULL,
    genomicContext_id = NULL,
    cached = TRUE,
    mode = "Fat"
) {
  BASE_URL <- "https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphisms"

  if (is.null(rs_id) ||
      trimws(as.character(rs_id)) == "") {
    stop("rs_id cannot be empty")
  }

  query_url <- paste0(
    BASE_URL, "/",
    URLencode(enc2utf8(as.character(rs_id)), reserved = TRUE),
    "/genomic-contexts"
  )
  if (!is.null(genomicContext_id) &&
      trimws(as.character(genomicContext_id)) != "") {
    query_url <- paste0(
      query_url, "/",
      URLencode(enc2utf8(as.character(genomicContext_id)), reserved = TRUE)
    )
  }

  return(
    GenomicContext(
      data = get_genomicContext(query_url, cached = cached),
      mode = mode
    )
  )
}
