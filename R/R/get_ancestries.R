#' Get Ancestry data from the server
#'
#' This function retrieves ancestry data from the GWAS API for a given study
#' accession identifier and optional ancestry identifier.
#'
#' @param accession_id The study's GWAS Catalog accession ID
#' @param ancestry_id primary identifier of ancestry table
#' @param cached Whether or not to try to get data from the cache
#' @param mode Fat or Thin. Specifies the mode of the returned object
#'
#' @return An Ancestry object. Attributes of type data.frame have hierarchical
#'   dependencies
#'
#' @examples
#' \dontrun{
#' ch <- get_ancestries(accession_id = "GCST000854")
#' }
#'
#' @importFrom utils URLencode
#'
#' @export
#'
get_ancestries <- function(
    accession_id = NULL,
    ancestry_id = NULL,
    cached = TRUE,
    mode = "Fat"
) {
  BASE_URL <- "https://www.ebi.ac.uk/gwas/rest/api/v2/studies"

  if (is.null(accession_id) ||
      trimws(as.character(accession_id)) == "") {
    stop("accession_id cannot be empty")
  }

  query_url <- paste0(
    BASE_URL, "/",
    URLencode(enc2utf8(as.character(accession_id)), reserved = TRUE),
    "/ancestries"
  )
  if (!is.null(ancestry_id) &&
      trimws(as.character(ancestry_id)) != "") {
    query_url <- paste0(
      query_url, "/",
      URLencode(enc2utf8(as.character(ancestry_id)), reserved = TRUE)
    )
  }

  return(
    Ancestry(
      data = get_ancestry(query_url, cached = cached),
      mode = mode
    )
  )
}
