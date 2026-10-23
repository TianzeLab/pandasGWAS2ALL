



# 缓存变量（在zzz.R中定义）

ask_yes_no_question <- function(question) {
  yes_no_answer <- ""
  while (toupper(yes_no_answer) != "YES" &&
         toupper(yes_no_answer) != "NO") {
    yes_no_answer <- readline(question)
  }
  return(toupper(yes_no_answer))
}


#' Get data from GWAS API
#'
#' This function retrieves data from the GWAS API, handling caching, pagination,
#' composite key extraction from self links and removal of link fields.
#'
#' @param url The URL to request data from
#' @param cache_impl A cachem cache object used for caching responses
#' @param cached Logical value indicating whether to use cached responses
#' @param type The type of data to retrieve (e.g., "associations", "studies")
#'
#' @return A list of data objects retrieved from the API
#'
#' @importFrom httr2 request req_headers req_perform resp_status
#' @importFrom httr2 resp_body_string req_cookie_preserve req_retry req_error
#' @importFrom jsonlite fromJSON
#' @importFrom progress progress_bar
#' @importFrom cachem cache_mem
#' @importFrom rlang hash
#'
#' @noRd
#'
get_data <- function(url,
                     cache_impl = NULL,
                     cached = TRUE,
                     type = NULL) {
  if (is.null(cache_impl)) {
    cache_impl <- cachem::cache_mem()
  }

  cookie_path <- tempfile()

  # Check cache. Python keys on the exact URL; cachem only permits lowercase
  # letters/numbers/_/- as keys, so use the URL's deterministic hash
  if (cached && cache_impl$exists(hash(url))) {
    r <- cache_impl$get(hash(url))
  } else {
    r <- httr2::req_perform(
      httr2::req_error(
        httr2::req_retry(
          httr2::req_headers(
            httr2::req_cookie_preserve(
              httr2::request(url),
              cookie_path
            ),
            Accept = "application/json"
          ),
          max_tries = 6
        ),
        is_error = ~FALSE
      )
    )
  }

  if (httr2::resp_status(r) == 200) {
    if (is.null(r$body) || length(r$body) == 0) {
      return(list())
    }
    content <- httr2::resp_body_string(r)

    parsed_data <- jsonlite::fromJSON(content, simplifyVector = FALSE)
    if (length(parsed_data) == 0) {
      return(list())
    }
    cache_impl$set(hash(url), r)

    if (!is.null(parsed_data[["_embedded"]])) {
      results_list <- parsed_data[["_embedded"]][[type]]
      if (is.null(results_list)) {
        results_list <- list()
      }

      # Pagination requires both _links and page information
      if (!is.null(parsed_data[["_links"]]) &&
          !is.null(parsed_data[["page"]])) {
        total_elements <- parsed_data[["page"]][["totalElements"]]

        pb <- NULL
        if (total_elements > 0) {
          pb <- progress::progress_bar$new(
            format = "  downloading [:bar] :current/:total eta: :eta",
            total = total_elements,
            clear = FALSE,
            width = 60
          )
          pb$update(length(results_list) / total_elements)
        }

        next_link <- parsed_data[["_links"]][["next"]]

        while (!is.null(next_link)) {
          next_url <- next_link[["href"]]

          if (cached && cache_impl$exists(hash(next_url))) {
            r <- cache_impl$get(hash(next_url))
          } else {
            r <- httr2::req_perform(
              httr2::req_error(
                httr2::req_retry(
                  httr2::req_headers(
                    httr2::req_cookie_preserve(
                      httr2::request(next_url),
                      cookie_path
                    ),
                    Accept = "application/json"
                  ),
                  max_tries = 6
                ),
                is_error = ~FALSE
              )
            )
          }

          if (httr2::resp_status(r) == 200) {
            parsed_next <- jsonlite::fromJSON(
              httr2::resp_body_string(r),
              simplifyVector = FALSE
            )
            page_results <- parsed_next[["_embedded"]][[type]]
            if (is.null(page_results)) {
              page_results <- list()
            }
            results_list <- c(results_list, page_results)
            cache_impl$set(hash(next_url), r)
            if (!is.null(pb)) {
              pb$update(length(results_list) / total_elements)
            }
            next_link <- parsed_next[["_links"]][["next"]]
          } else if (httr2::resp_status(r) == 404) {
            return(list())
          } else {
            stop(paste(
              "The request for", next_url,
              "failed: response code was", httr2::resp_status(r),
              "content was\n", httr2::resp_body_string(r)
            ))
          }
        }

        if (!is.null(pb)) {
          pb$terminate()
        }
      }

      # Clean every dict item in place:
      # extract composite keys from self href, drop links and _links
      for (item_idx in seq_along(results_list)) {
        item <- results_list[[item_idx]]
        if (is.list(item) &&
            (length(item) == 0 || !is.null(names(item)))) {
          href <- item[["_links"]][["self"]][["href"]]

          if (type == "ancestries" && !is.null(href)) {
            match_parts <- regmatches(
              href,
              regexec("/studies/([^/]+)/ancestries/([^/?#]+)", href)
            )[[1]]
            if (length(match_parts) == 3) {
              item[["accession_id"]] <- match_parts[2]
              item[["ancestry_id"]] <- match_parts[3]
            }
          } else if (type == "genomic_contexts" && !is.null(href)) {
            match_parts <- regmatches(
              href,
              regexec(
                "/single-nucleotide-polymorphisms/([^/]+)/genomic-contexts/([^/?#]+)",
                href
              )
            )[[1]]
            if (length(match_parts) == 3) {
              item[["rs_id"]] <- match_parts[2]
              item[["genomicContext_id"]] <- match_parts[3]
            }
          } else if (type == "loci" && !is.null(href)) {
            match_parts <- regmatches(
              href,
              regexec("/associations/([^/]+)/loci/([^/?#]+)", href)
            )[[1]]
            if (length(match_parts) == 3) {
              item[["association_id"]] <- match_parts[2]
              item[["locus_id"]] <- match_parts[3]
            }
          }

          item[["links"]] <- NULL
          item[["_links"]] <- NULL
          results_list[[item_idx]] <- item
        }
      }

      return(results_list)
    } else {
      # A top-level object carrying a page field but no _embedded yields nothing
      if (!is.null(parsed_data[["page"]])) {
        return(list())
      }

      # Clean the single object in place
      href <- parsed_data[["_links"]][["self"]][["href"]]

      if (type == "ancestries" && !is.null(href)) {
        match_parts <- regmatches(
          href,
          regexec("/studies/([^/]+)/ancestries/([^/?#]+)", href)
        )[[1]]
        if (length(match_parts) == 3) {
          parsed_data[["accession_id"]] <- match_parts[2]
          parsed_data[["ancestry_id"]] <- match_parts[3]
        }
      } else if (type == "genomic_contexts" && !is.null(href)) {
        match_parts <- regmatches(
          href,
          regexec(
            "/single-nucleotide-polymorphisms/([^/]+)/genomic-contexts/([^/?#]+)",
            href
          )
        )[[1]]
        if (length(match_parts) == 3) {
          parsed_data[["rs_id"]] <- match_parts[2]
          parsed_data[["genomicContext_id"]] <- match_parts[3]
        }
      } else if (type == "loci" && !is.null(href)) {
        match_parts <- regmatches(
          href,
          regexec("/associations/([^/]+)/loci/([^/?#]+)", href)
        )[[1]]
        if (length(match_parts) == 3) {
          parsed_data[["association_id"]] <- match_parts[2]
          parsed_data[["locus_id"]] <- match_parts[3]
        }
      }

      parsed_data[["links"]] <- NULL
      parsed_data[["_links"]] <- NULL
      return(list(parsed_data))
    }
  } else if (httr2::resp_status(r) == 404) {
    return(list())
  } else {
    stop(paste(
      "The request for", url,
      "failed: response code was", httr2::resp_status(r),
      "content was\n", httr2::resp_body_string(r)
    ))
  }
}

get_association <- function(url, cached = TRUE) {
  return(
    get_data(
      url,
      cache_impl = .pandasGWASr$association_cache,
      cached = cached,
      type = "associations"
    )
  )
}

get_study <- function(url, cached = TRUE) {
  return(
    get_data(
      url,
      cache_impl = .pandasGWASr$study_cache,
      cached = cached,
      type = "studies"
    )
  )
}

get_trait <- function(url, cached = TRUE) {
  return(
    get_data(
      url,
      cache_impl = .pandasGWASr$trait_cache,
      cached = cached,
      type = "efo_traits"
    )
  )
}

get_snp <- function(url, cached = TRUE) {
  return(
    get_data(
      url,
      cache_impl = .pandasGWASr$snp_cache,
      cached = cached,
      type = "snps"
    )
  )
}

get_ancestry <- function(url, cached = TRUE) {
  return(
    get_data(
      url,
      cache_impl = .pandasGWASr$ancestry_cache,
      cached = cached,
      type = "ancestries"
    )
  )
}

get_publication <- function(url, cached = TRUE) {
  return(
    get_data(
      url,
      cache_impl = .pandasGWASr$publication_cache,
      cached = cached,
      type = "publications"
    )
  )
}

get_genomicContext <- function(url, cached = TRUE) {
  return(
    get_data(
      url,
      cache_impl = .pandasGWASr$genomicContext_cache,
      cached = cached,
      type = "genomic_contexts"
    )
  )
}

get_gene <- function(url, cached = TRUE) {
  return(
    get_data(
      url,
      cache_impl = .pandasGWASr$gene_cache,
      cached = cached,
      type = "genes"
    )
  )
}

get_locus <- function(url, cached = TRUE) {
  return(
    get_data(
      url,
      cache_impl = .pandasGWASr$locus_cache,
      cached = cached,
      type = "loci"
    )
  )
}

#' Reinitialize cache
#'
#' This function reinitializes the cache for the specified data type or all caches.
#'
#' @param field The field to reinitialize the cache for, or "All" for all caches
#' @param maxsize The maximum size of the cache
#' @param ttl The time to live for the cache in seconds
#'
#' @export
#'
reinit_cache <- function(field = "All", maxsize = 2 * 1024, ttl = 60 * 60 * 24) {
  fields <- c("Association", "Study", "Trait", "Snp", "Ancestry", "Publication", "GenomicContext", "Gene", "Locus")
  fields_and_all <- c("All", "Association", "Study", "Trait", "Snp", "Ancestry", "Publication", "GenomicContext", "Gene", "Locus")
  
  if (!field %in% fields_and_all) {
    stop(paste("The field must be one of", paste(fields_and_all, collapse = ", ")))
  }
  
  if (field == "Association") {
    .pandasGWASr$association_cache <- cache_mem(max_age = ttl, max_size = maxsize)
  } else if (field == "Study") {
    .pandasGWASr$study_cache <- cache_mem(max_age = ttl, max_size = maxsize)
  } else if (field == "Trait") {
    .pandasGWASr$trait_cache <- cache_mem(max_age = ttl, max_size = maxsize)
  } else if (field == "Snp") {
    .pandasGWASr$snp_cache <- cache_mem(max_age = ttl, max_size = maxsize)
  } else if (field == "Ancestry") {
    .pandasGWASr$ancestry_cache <- cache_mem(max_age = ttl, max_size = maxsize)
  } else if (field == "Publication") {
    .pandasGWASr$publication_cache <- cache_mem(max_age = ttl, max_size = maxsize)
  } else if (field == "GenomicContext") {
    .pandasGWASr$genomicContext_cache <- cache_mem(max_age = ttl, max_size = maxsize)
  } else if (field == "Gene") {
    .pandasGWASr$gene_cache <- cache_mem(max_age = ttl, max_size = maxsize)
  } else if (field == "Locus") {
    .pandasGWASr$locus_cache <- cache_mem(max_age = ttl, max_size = maxsize)
  } else { # All
    .pandasGWASr$association_cache <- cache_mem(max_age = ttl, max_size = maxsize)
    .pandasGWASr$study_cache <- cache_mem(max_age = ttl, max_size = maxsize)
    .pandasGWASr$trait_cache <- cache_mem(max_age = ttl, max_size = maxsize)
    .pandasGWASr$snp_cache <- cache_mem(max_age = ttl, max_size = maxsize)
    .pandasGWASr$ancestry_cache <- cache_mem(max_age = ttl, max_size = maxsize)
    .pandasGWASr$publication_cache <- cache_mem(max_age = ttl, max_size = maxsize)
    .pandasGWASr$genomicContext_cache <- cache_mem(max_age = ttl, max_size = maxsize)
    .pandasGWASr$gene_cache <- cache_mem(max_age = ttl, max_size = maxsize)
    .pandasGWASr$locus_cache <- cache_mem(max_age = ttl, max_size = maxsize)
  }
}
