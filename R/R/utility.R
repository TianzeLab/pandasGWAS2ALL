#' Test if the API is available
#'
#' Test if the API is available. If the return value is not TRUE, it means that
#' the API is not available.
#'
#' @return TRUE or FALSE
#'
#' @examples
#' \dontrun{
#' is_available <- is_API_available()
#' }
#'
#' @importFrom httr2 request req_headers req_perform resp_status req_retry
#' @importFrom httr2 req_error req_cookie_preserve
#'
#' @export
#'
is_API_available <- function() {
  tryCatch({
    cookie_path <- tempfile()
    r <- httr2::req_perform(
      httr2::req_error(
        httr2::req_retry(
          httr2::req_headers(
            httr2::req_cookie_preserve(
              httr2::request("https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits/"),
              cookie_path
            ),
            Accept = "application/json"
          ),
          max_tries = 6
        ),
        is_error = ~FALSE
      )
    )
    if (httr2::resp_status(r) == 200) {
      return(TRUE)
    }
    return(FALSE)
  }, error = function(e) {
    return(FALSE)
  })
}
