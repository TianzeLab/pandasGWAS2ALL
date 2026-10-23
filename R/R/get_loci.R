#' Get Locus data from the server
#'
#' This function retrieves locus data from the GWAS API for a given association
#' identifier and optional locus identifier.
#'
#' @param association_id primary identifier of Association table
#' @param locus_id primary identifier of Locus table
#' @param cached Whether or not to try to get data from the cache
#' @param mode Fat or Thin. Specifies the mode of the returned object
#'
#' @return A Locus object. Attributes of type data.frame have hierarchical
#'   dependencies
#'
#' @examples
#' \dontrun{
#' ch <- get_loci(association_id = "169106883")
#' }
#'
#' @importFrom utils URLencode
#'
#' @export
#'
get_loci <- function(
    association_id = NULL,
    locus_id = NULL,
    cached = TRUE,
    mode = "Fat"
) {
  BASE_URL <- "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"

  if (is.null(association_id) ||
      trimws(as.character(association_id)) == "") {
    stop("association_id cannot be empty")
  }

  query_url <- paste0(
    BASE_URL, "/",
    URLencode(enc2utf8(as.character(association_id)), reserved = TRUE),
    "/loci"
  )
  if (!is.null(locus_id) &&
      trimws(as.character(locus_id)) != "") {
    query_url <- paste0(
      query_url, "/",
      URLencode(enc2utf8(as.character(locus_id)), reserved = TRUE)
    )
  }

  return(
    Locus(
      data = get_locus(query_url, cached = cached),
      mode = mode
    )
  )
}
