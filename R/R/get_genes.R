#' Get Gene data from the server
#'
#' This function retrieves gene data from the GWAS API, optionally filtered by
#' gene name.
#'
#' @param gene_name Gene Name
#' @param cached Whether or not to try to get data from the cache
#' @param mode Fat or Thin. Specifies the mode of the returned object
#'
#' @return A Gene object. Attributes of type data.frame have hierarchical
#'   dependencies
#'
#' @examples
#' \dontrun{
#' ch <- get_genes(gene_name = "HBS1L")
#' }
#'
#' @export
#'
get_genes <- function(
    gene_name = NULL,
    cached = TRUE,
    mode = "Fat"
) {
  BASE_URL <- "https://www.ebi.ac.uk/gwas/rest/api/v2/genes"

  if (is.null(gene_name)) {
    query_url <- BASE_URL
  } else {
    query_url <- paste0(BASE_URL, "/", as.character(gene_name))
  }

  return(
    Gene(
      data = get_gene(query_url, cached = cached),
      mode = mode
    )
  )
}
