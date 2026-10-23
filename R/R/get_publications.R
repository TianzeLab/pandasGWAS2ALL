#' Get Publication data from the server
#'
#' This function retrieves publication data from the GWAS API based on various
#' search criteria.
#'
#' @param pubmed_id pubmedId of the publication
#' @param title Title of manuscript the GWAS is included in
#' @param first_author Last name and initials of first author
#' @param cached Whether or not to try to get data from the cache
#' @param mode Fat or Thin. Specifies the mode of the returned object
#'
#' @return A Publication object. Attributes of type data.frame have hierarchical
#'   dependencies
#'
#' @examples
#' \dontrun{
#' ch <- get_publications(pubmed_id = "35241825")
#' }
#'
#' @importFrom utils URLencode
#'
#' @export
#'
get_publications <- function(
    pubmed_id = NULL,
    title = NULL,
    first_author = NULL,
    cached = TRUE,
    mode = "Fat"
) {
  BASE_URL <- "https://www.ebi.ac.uk/gwas/rest/api/v2/publications"

  params <- list(
    pubmed_id = pubmed_id,
    title = title,
    first_author = first_author
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
    Publication(
      data = get_publication(query_url, cached = cached),
      mode = mode
    )
  )
}
