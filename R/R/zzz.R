# 包的初始化文件

# 定义包全局环境
.pandasGWASr <- new.env(parent = emptyenv())

#' @importFrom cachem cache_mem
#' @noRd
.onLoad <- function(libname, pkgname) {
  # 设置缓存过期时间为24小时（单位：秒）
  ttl <- 60 * 60 * 24
  
  # 创建包全局内存TTL缓存
  .pandasGWASr$association_cache <- cache_mem(max_age = ttl)
  .pandasGWASr$study_cache <- cache_mem(max_age = ttl)
  .pandasGWASr$trait_cache <- cache_mem(max_age = ttl)
  .pandasGWASr$snp_cache <- cache_mem(max_age = ttl)
  .pandasGWASr$ancestry_cache <- cache_mem(max_age = ttl)
  .pandasGWASr$publication_cache <- cache_mem(max_age = ttl)
  .pandasGWASr$genomicContext_cache <- cache_mem(max_age = ttl)
  .pandasGWASr$gene_cache <- cache_mem(max_age = ttl)
  .pandasGWASr$locus_cache <- cache_mem(max_age = ttl)
}

#' Clear some or all of the cache pools
#'
#' This function clears some or all of the cache pools.
#'
#' @param field It can be one of the following: 'All', 'Association', 'Study', 'Trait', 'Snp', 'Ancestry', 'Publication', 'GenomicContext', 'Gene', 'Locus'
#'
#' @export
#'
clear_cache <- function(field = "All") {
  fields_and_all <- c("All", "Association", "Study", "Trait", "Snp", "Ancestry", "Publication", "GenomicContext", "Gene", "Locus")

  if (!field %in% fields_and_all) {
    stop(paste("The field must be one of", paste(fields_and_all, collapse = ", ")))
  }

  if (field == "Association") {
    .pandasGWASr$association_cache$reset()
  } else if (field == "Study") {
    .pandasGWASr$study_cache$reset()
  } else if (field == "Trait") {
    .pandasGWASr$trait_cache$reset()
  } else if (field == "Snp") {
    .pandasGWASr$snp_cache$reset()
  } else if (field == "Ancestry") {
    .pandasGWASr$ancestry_cache$reset()
  } else if (field == "Publication") {
    .pandasGWASr$publication_cache$reset()
  } else if (field == "GenomicContext") {
    .pandasGWASr$genomicContext_cache$reset()
  } else if (field == "Gene") {
    .pandasGWASr$gene_cache$reset()
  } else if (field == "Locus") {
    .pandasGWASr$locus_cache$reset()
  } else { # All
    .pandasGWASr$association_cache$reset()
    .pandasGWASr$study_cache$reset()
    .pandasGWASr$trait_cache$reset()
    .pandasGWASr$snp_cache$reset()
    .pandasGWASr$ancestry_cache$reset()
    .pandasGWASr$publication_cache$reset()
    .pandasGWASr$genomicContext_cache$reset()
    .pandasGWASr$gene_cache$reset()
    .pandasGWASr$locus_cache$reset()
  }
}
