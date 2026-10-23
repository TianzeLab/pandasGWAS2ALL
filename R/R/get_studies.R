#' Get Study data from the server
#'
#' This function retrieves study data from the GWAS API based on various search
#' criteria.
#'
#' @param accession_id The study's GWAS Catalog accession ID
#' @param pubmed_id pubmedId of the publication
#' @param disease_trait Free text description of the trait investigated in this
#'   study
#' @param full_pvalue_set Whether full summary statistics are available for this
#'   study
#' @param efo_id The trait URI shortform
#' @param efo_trait The trait name or label
#' @param cohort Discovery stage cohorts used in this study
#' @param gxe Whether the study investigates a gene-environment interaction
#' @param ancestral_group Ancestry category group label, assigned to reduce
#'   complexity within the data sets and place samples in context
#' @param no_of_individuals Number of individuals in this sample
#' @param show_child_trait display entities for descendants of a parent Efo
#'   Trait Term
#' @param mapped_gene Gene(s) overlapping the variant. If a variant is
#'   intergenic, the closest 5' and 3' genes are listed
#' @param extended_geneset Show extended matching genes in addition to the mapped
#'   genes
#' @param cached Whether or not to try to get data from the cache
#' @param mode Fat or Thin. Specifies the mode of the returned object
#'
#' @return A Study object. Attributes of type data.frame have hierarchical
#'   dependencies
#'
#' @examples
#' \dontrun{
#' ch <- get_studies(accession_id = "GCST000854")
#' }
#'
#' @importFrom utils URLencode
#'
#' @export
#'
get_studies <- function(
    accession_id = NULL,
    pubmed_id = NULL,
    disease_trait = NULL,
    full_pvalue_set = NULL,
    efo_id = NULL,
    efo_trait = NULL,
    cohort = NULL,
    gxe = NULL,
    ancestral_group = NULL,
    no_of_individuals = NULL,
    show_child_trait = NULL,
    mapped_gene = NULL,
    extended_geneset = NULL,
    cached = TRUE,
    mode = "Fat"
) {
  BASE_URL <- "https://www.ebi.ac.uk/gwas/rest/api/v2/studies"

  params <- list(
    accession_id = accession_id,
    pubmed_id = pubmed_id,
    disease_trait = disease_trait,
    full_pvalue_set = full_pvalue_set,
    efo_id = efo_id,
    efo_trait = efo_trait,
    cohort = cohort,
    gxe = gxe,
    ancestral_group = ancestral_group,
    no_of_individuals = no_of_individuals,
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

  return(
    Study(
      data = get_study(query_url, cached = cached),
      mode = mode
    )
  )
}
