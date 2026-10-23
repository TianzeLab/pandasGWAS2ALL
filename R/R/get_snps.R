#' Get Snp data from the server
#'
#' This function retrieves SNP data from the GWAS API based on various search
#' criteria.
#'
#' @param rs_id SNP rsID; if a haplotype it may include more than one rs number
#'   (multiple SNPs comprising the haplotype)
#' @param bp_location search for SNPs via their base pair location
#' @param bp_start Search for SNPs on a certain chromosome and within a given
#'   base pair location range
#' @param bp_end Search for SNPs on a certain chromosome and within a given base
#'   pair location range
#' @param pubmed_id pubmedId of the publication
#' @param chromosome Search for SNPs on a certain chromosome
#' @param mapped_gene Gene(s) overlapping the variant. If a variant is
#'   intergenic, the closest 5' and 3' genes are listed
#' @param extended_geneset Show extended matching genes in addition to the mapped
#'   genes
#' @param cached Whether or not to try to get data from the cache
#' @param mode Fat or Thin. Specifies the mode of the returned object
#'
#' @return A Snp object. Attributes of type data.frame have hierarchical
#'   dependencies
#'
#' @examples
#' \dontrun{
#' ch <- get_snps(rs_id = "rs3093017")
#' }
#'
#' @importFrom utils URLencode
#'
#' @export
#'
get_snps <- function(
    rs_id = NULL,
    bp_location = NULL,
    bp_start = NULL,
    bp_end = NULL,
    pubmed_id = NULL,
    chromosome = NULL,
    mapped_gene = NULL,
    extended_geneset = NULL,
    cached = TRUE,
    mode = "Fat"
) {
  BASE_URL <- "https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphisms"

  params <- list(
    rs_id = rs_id,
    bp_location = bp_location,
    bp_start = bp_start,
    bp_end = bp_end,
    pubmed_id = pubmed_id,
    chromosome = chromosome,
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
    Snp(
      data = get_snp(query_url, cached = cached),
      mode = mode
    )
  )
}
