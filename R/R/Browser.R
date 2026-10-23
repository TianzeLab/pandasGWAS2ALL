#' Open a PubMed citation in the browser
#'
#' This function launches the web browser and opens a tab for a PubMed citation.
#'
#' @param pubmed_id A PubMed identifier
#'
#' @return None (invisible NULL). A browser tab is opened as a side effect.
#'
#' @examples
#' \dontrun{
#' open_in_pubmed("26301688")
#' }
#'
#' @export
#'
open_in_pubmed <- function(pubmed_id) {
  utils::browseURL(
    paste0("https://pubmed.ncbi.nlm.nih.gov/", pubmed_id)
  )
}


#' Open a dbSNP page in the browser
#'
#' This function launches the web browser at dbSNP and opens a tab for a SNP
#' identifier.
#'
#' @param variant_id A variant (Single Nucleotide Polymorphism) identifier
#'
#' @return None (invisible NULL). A browser tab is opened as a side effect.
#'
#' @examples
#' \dontrun{
#' open_in_dbsnp("rs56261590")
#' }
#'
#' @export
#'
open_in_dbsnp <- function(variant_id) {
  utils::browseURL(
    paste0("https://www.ncbi.nlm.nih.gov/snp/", variant_id)
  )
}


#' Open a GTEx Portal page in the browser
#'
#' This function launches the web browser at the GTEx Portal and opens a tab for
#' a SNP identifier.
#'
#' @param variant_id A variant (Single Nucleotide Polymorphism) identifier
#'
#' @return None (invisible NULL). A browser tab is opened as a side effect.
#'
#' @examples
#' \dontrun{
#' open_in_gtex("rs56261590")
#' }
#'
#' @export
#'
open_in_gtex <- function(variant_id) {
  utils::browseURL(
    paste0("https://gtexportal.org/home/snp/", variant_id)
  )
}


#' Open a Study page in the GWAS Catalog browser
#'
#' Browse GWAS Catalog entity Study from the GWAS Web Graphical User Interface.
#'
#' @param study_id A Study identifier
#'
#' @return None (invisible NULL). A browser tab is opened as a side effect.
#'
#' @examples
#' \dontrun{
#' open_study_in_gwas_catalog("GCST000016")
#' }
#'
#' @export
#'
open_study_in_gwas_catalog <- function(study_id) {
  utils::browseURL(
    paste0("https://www.ebi.ac.uk/gwas/studies/", study_id)
  )
}


#' Open a Variant page in the GWAS Catalog browser
#'
#' Browse GWAS Catalog entity Single Nucleotide Polymorphism from the GWAS Web
#' Graphical User Interface.
#'
#' @param variant_id A variant (Single Nucleotide Polymorphism) identifier
#'
#' @return None (invisible NULL). A browser tab is opened as a side effect.
#'
#' @examples
#' \dontrun{
#' open_variant_in_gwas_catalog("rs146992477")
#' }
#'
#' @export
#'
open_variant_in_gwas_catalog <- function(variant_id) {
  utils::browseURL(
    paste0("https://www.ebi.ac.uk/gwas/variants/", variant_id)
  )
}


#' Open a Trait page in the GWAS Catalog browser
#'
#' Browse GWAS Catalog entity Trait from the GWAS Web Graphical User Interface.
#'
#' @param efo_id An EFO Trait identifier
#'
#' @return None (invisible NULL). A browser tab is opened as a side effect.
#'
#' @examples
#' \dontrun{
#' open_trait_in_gwas_catalog("EFO_0004884")
#' }
#'
#' @export
#'
open_trait_in_gwas_catalog <- function(efo_id) {
  utils::browseURL(
    paste0("https://www.ebi.ac.uk/gwas/efotraits/", efo_id)
  )
}


#' Open a Gene page in the GWAS Catalog browser
#'
#' Browse GWAS Catalog entity Gene from the GWAS Web Graphical User Interface.
#'
#' @param gene_name Gene name
#'
#' @return None (invisible NULL). A browser tab is opened as a side effect.
#'
#' @examples
#' \dontrun{
#' open_gene_in_gwas_catalog("DPP6")
#' }
#'
#' @export
#'
open_gene_in_gwas_catalog <- function(gene_name) {
  utils::browseURL(
    paste0("https://www.ebi.ac.uk/gwas/genes/", gene_name)
  )
}


#' Open a Region page in the GWAS Catalog browser
#'
#' Browse GWAS Catalog entity Region from the GWAS Web Graphical User Interface.
#'
#' @param region_name_or_location Region name or chromosome and base pair location
#'   on the reference genome
#'
#' @return None (invisible NULL). A browser tab is opened as a side effect.
#'
#' @examples
#' \dontrun{
#' # region name
#' open_region_in_gwas_catalog("2q37.1")
#' # location
#' open_region_in_gwas_catalog("chr2:230100001-234700000")
#' }
#'
#' @export
#'
open_region_in_gwas_catalog <- function(region_name_or_location) {
  utils::browseURL(
    paste0("https://www.ebi.ac.uk/gwas/regions/", region_name_or_location)
  )
}


#' Open a Publication page in the GWAS Catalog browser
#'
#' Browse GWAS Catalog entity Publication from the GWAS Web Graphical User
#' Interface.
#'
#' @param pubmed_id A PubMed identifier
#'
#' @return None (invisible NULL). A browser tab is opened as a side effect.
#'
#' @examples
#' \dontrun{
#' open_publication_in_gwas_catalog("25533513")
#' }
#'
#' @export
#'
open_publication_in_gwas_catalog <- function(pubmed_id) {
  utils::browseURL(
    paste0("https://www.ebi.ac.uk/gwas/publications/", pubmed_id)
  )
}
