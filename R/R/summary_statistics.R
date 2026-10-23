#' Summary statistics submodule
#'
#' A set of functions for easy retrieval of summary statistics data based on FTP
#' data. This environment mimics a Python submodule and contains four
#' functions: \code{search}, \code{browser}, \code{download}, and \code{parse}.
#' Call them via \code{summary_statistics$search(...)},
#' \code{summary_statistics$browser(...)},
#' \code{summary_statistics$download(...)}, and
#' \code{summary_statistics$parse(...)}.
#'
#' @details
#' \describe{
#'   \item{search}{Search for where data is stored based on indexing and query
#'     criteria. Returns a data.frame that records where the data is stored.}
#'   \item{browser}{See where the data is stored in the browser.}
#'   \item{download}{Download FTP data to directory ~/pandasgwas_home.}
#'   \item{parse}{Resolve the specified data from ~/pandasgwas_home to a
#'     data.frame recording summary statistics.}
#' }
#'
#' @examples
#' \dontrun{
#' search_DF <- summary_statistics$search(
#'   PubMed_id = "27918534",
#'   study_accession_id = "GCST003966"
#' )
#' summary_statistics$browser(search_DF)
#' summary_statistics$download(search_DF)
#' df <- summary_statistics$parse(search_DF)
#' }
#'
#' @importFrom utils download.file
#' @importFrom httr2 request req_retry req_method req_perform resp_headers
#'
#' @export
#'
summary_statistics <- new.env(parent = baseenv())


summary_statistics$search <- function(PubMed_id = NULL,
                                      study_accession_id = NULL,
                                      EFO_trait_id = NULL,
                                      online_index = FALSE) {
  host <- "ftp.ebi.ac.uk"
  home_path <- file.path(path.expand("~"), "pandasgwas_home")

  if (is.null(PubMed_id) && is.null(study_accession_id) && is.null(EFO_trait_id)) {
    warning("Since there are no input conditions, all index values will be returned.")
  }

  if (online_index) {
    dir.create(home_path, recursive = TRUE, showWarnings = FALSE)
    download.file(
      paste0("https://", host, "/pub/databases/gwas/summary_statistics/harmonised_list.txt"),
      file.path(home_path, "harmonised_list.txt"),
      mode = "wb",
      quiet = TRUE
    )
    lines <- readLines(file.path(home_path, "harmonised_list.txt"))
  } else {
    if (is.null(.pandasGWASr$harmonised_df)) {
      lines <- readLines(system.file("extdata", "harmonised_list.txt", package = "pandasGWASr"))
    } else {
      lines <- NULL
    }
  }

  if (!is.null(lines)) {
    parsed <- vapply(lines, function(line) {
      match1 <- regmatches(line, regexec("\\.(.*/)((.*)-(.*)-(.*).h.tsv.gz)", line))[[1]]
      if (length(match1) == 6) {
        return(c(
          paste0("/pub/databases/gwas/summary_statistics", match1[2]),
          match1[3],
          match1[4],
          match1[5],
          match1[6]
        ))
      }
      match2 <- regmatches(line, regexec("\\.(.*/)((GCST.*).h.tsv.gz)", line))[[1]]
      return(c(
        paste0("/pub/databases/gwas/summary_statistics", match2[2]),
        match2[3],
        NA_character_,
        match2[4],
        NA_character_
      ))
    }, character(5))
    parsed <- t(parsed)
    filter_df <- data.frame(
      raw = lines,
      dir = parsed[, 1],
      file_name = parsed[, 2],
      PubMed_id = parsed[, 3],
      study_accession_id = parsed[, 4],
      EFO_trait_id = parsed[, 5],
      stringsAsFactors = FALSE
    )
    if (!online_index) {
      .pandasGWASr$harmonised_df <- filter_df
    }
  } else {
    filter_df <- .pandasGWASr$harmonised_df
  }

  if (!is.null(PubMed_id)) {
    filter_df <- filter_df[!is.na(filter_df$PubMed_id) & filter_df$PubMed_id == PubMed_id, , drop = FALSE]
  }
  if (!is.null(study_accession_id)) {
    filter_df <- filter_df[!is.na(filter_df$study_accession_id) & filter_df$study_accession_id == study_accession_id, , drop = FALSE]
  }
  if (!is.null(EFO_trait_id)) {
    filter_df <- filter_df[!is.na(filter_df$EFO_trait_id) & filter_df$EFO_trait_id == EFO_trait_id, , drop = FALSE]
  }

  rownames(filter_df) <- NULL
  return(filter_df)
}


summary_statistics$browser <- function(search_DF, interactive = TRUE) {
  host <- "ftp.ebi.ac.uk"
  if (nrow(search_DF) > 5 && interactive) {
    answer <- ask_yes_no_question(
      paste0("You will have more than ", nrow(search_DF), " pages opened.\r\nDo you still want to proceed? (Yes or No)")
    )
    if (answer == "NO") {
      return(invisible(NULL))
    }
  }
  for (i in seq_len(nrow(search_DF))) {
    utils::browseURL(paste0("https://", host, search_DF$dir[i]))
  }
  return(invisible(NULL))
}


summary_statistics$download <- function(search_DF) {
  host <- "ftp.ebi.ac.uk"
  home_path <- file.path(path.expand("~"), "pandasgwas_home")
  dir.create(home_path, recursive = TRUE, showWarnings = FALSE)

  for (i in seq_len(nrow(search_DF))) {
    dir_path <- search_DF$dir[i]
    file_name <- search_DF$file_name[i]
    url <- paste0("https://", host, dir_path, file_name)
    local_path <- file.path(home_path, file_name)

    skip <- FALSE
    if (file.exists(local_path)) {
      local_size <- file.size(local_path)
      if (local_size > 0) {
        r <- tryCatch({
          httr2::req_perform(
            httr2::req_method(
              httr2::req_retry(
                httr2::request(url),
                max_tries = 6
              ),
              "HEAD"
            )
          )
        }, error = function(e) NULL)
        if (!is.null(r)) {
          headers <- httr2::resp_headers(r)
          online_size <- as.numeric(headers[["content-length"]])
          if (!is.na(online_size) && online_size == local_size) {
            cat(sprintf("[SKIP]: %s has been downloaded in %s\n", file_name, home_path))
            skip <- TRUE
          }
        }
      }
    }

    if (!skip) {
      download.file(url, local_path, mode = "wb", quiet = FALSE)
      cat(sprintf("%s (%.2f MB) has been downloaded in %s\n",
                  file_name, file.size(local_path) / 1024 / 1024, home_path))
    }
  }
  return(invisible(NULL))
}


summary_statistics$parse <- function(search_DF, interactive = TRUE) {
  home_path <- file.path(path.expand("~"), "pandasgwas_home")

  total <- sum(file.size(file.path(home_path, search_DF$file_name)), na.rm = TRUE) / 1024 / 1024 / 1024
  if (total > 1 && interactive) {
    answer <- ask_yes_no_question(
      paste0("You will load more than ", sprintf("%.2f", total), " GB data.\r\nDo you still want to proceed? (Yes or No)")
    )
    if (answer == "NO") {
      return(invisible(NULL))
    }
  }

  needed_cols <- c("variant_id", "p_value", "chromosome", "base_pair_location",
                   "odds_ratio", "ci_lower", "ci_upper", "beta", "standard_error",
                   "effect_allele", "other_allele", "effect_allele_frequency",
                   "hm_variant_id", "hm_odds_ratio", "hm_ci_lower", "hm_ci_upper",
                   "hm_beta", "hm_effect_allele", "hm_other_allele",
                   "hm_effect_allele_frequency", "hm_code")

  dfs <- vector("list", nrow(search_DF))
  for (i in seq_len(nrow(search_DF))) {
    path <- file.path(home_path, search_DF$file_name[i])
    con <- gzfile(path)
    header <- strsplit(readLines(con, n = 1), "\t")[[1]]
    close(con)
    col_classes <- ifelse(header %in% needed_cols, NA, "NULL")
    df <- read.table(gzfile(path), header = TRUE, sep = "\t", quote = "",
                     colClasses = col_classes, stringsAsFactors = FALSE,
                     comment.char = "", check.names = FALSE)
    df$PubMed_id <- search_DF$PubMed_id[i]
    df$study_accession_id <- search_DF$study_accession_id[i]
    df$EFO_trait_id <- search_DF$EFO_trait_id[i]
    dfs[[i]] <- df
  }

  result <- do.call(rbind, dfs)
  row.names(result) <- NULL
  return(result)
}
