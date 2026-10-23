#' Bind together GWAS objects of the same type
#'
#' Binds together GWAS objects of the same type. Note that bind() preserves
#' duplicates whereas union() does not.
#'
#' @param a An object of a pandasGWASr S4 class
#' @param b An object of the same type as a
#'
#' @return An object of the same type as a
#'
#' @examples
#' \dontrun{
#' a <- get_snps(rs_id = "rs3093017")
#' b <- get_snps(rs_id = "rs7671745")
#' c <- bind(a, b)
#' }
#'
#' @export
#'
bind <- function(a, b) {
  return(a + b)
}


#' Return the data common to both A and B
#'
#' Returns the data common to both A and B, with no repetitions.
#'
#' @param a An object of a pandasGWASr S4 class
#' @param b An object of the same type as a
#'
#' @return An object of the same type as a
#'
#' @examples
#' \dontrun{
#' a <- get_snps(rs_id = "rs3093017")
#' b <- get_snps(rs_id = "rs7671745")
#' c <- intersect(a, b)
#' }
#'
#' @export
#'
intersect <- function(a, b) {
  return(a & b)
}


#' Return the data in A that is not in B
#'
#' Returns the data in A that is not in B, with no repetitions.
#'
#' @param a An object of a pandasGWASr S4 class
#' @param b An object of the same type as a
#'
#' @return An object of the same type as a
#'
#' @examples
#' \dontrun{
#' a <- get_snps(rs_id = "rs3093017")
#' b <- get_snps(rs_id = "rs7671745")
#' c <- set_diff(a, b)
#' }
#'
#' @export
#'
set_diff <- function(a, b) {
  return(a - b)
}


#' Return the symmetric difference of A and B
#'
#' Returns the data of A and B that are not in their intersection (the symmetric
#' difference), with no repetitions.
#'
#' @param a An object of a pandasGWASr S4 class
#' @param b An object of the same type as a
#'
#' @return An object of the same type as a
#'
#' @examples
#' \dontrun{
#' a <- get_snps(rs_id = "rs3093017")
#' b <- get_snps(rs_id = "rs7671745")
#' c <- set_xor(a, b)
#' }
#'
#' @export
#'
set_xor <- function(a, b) {
  return(a ^ b)
}


#' Return the combined data from A and B
#'
#' Returns the combined data from A and B with no repetitions.
#'
#' @param a An object of a pandasGWASr S4 class
#' @param b An object of the same type as a
#'
#' @return An object of the same type as a
#'
#' @examples
#' \dontrun{
#' a <- get_snps(rs_id = "rs3093017")
#' b <- get_snps(rs_id = "rs7671745")
#' c <- union(a, b)
#' }
#'
#' @export
#'
union <- function(a, b) {
  return(a | b)
}


#' Check if the raw data of a and b are equal
#'
#' Check if the raw data of a and b are equal.
#'
#' @param a An object of a pandasGWASr S4 class
#' @param b An object of the same type as a
#'
#' @return TRUE, if a is equal to b. FALSE, if a is not equal to b.
#'
#' @examples
#' \dontrun{
#' a <- get_snps(rs_id = "rs3093017")
#' b <- get_snps(rs_id = "rs7671745")
#' c <- set_equal(a, b)
#' }
#'
#' @export
#'
set_equal <- function(a, b) {
  return(a == b)
}
