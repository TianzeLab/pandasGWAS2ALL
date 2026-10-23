import requests
import json
import progressbar
from typing import List, Dict
import re
from requests.adapters import HTTPAdapter
from cachetools import TTLCache

fields = ['Association', 'Study', 'Trait', 'Snp', 'Ancestry', 'Publication', 'GenomicContext', 'Gene', 'Locus']
fields_and_all = ['All', 'Association', 'Study', 'Trait', 'Snp', 'Ancestry', 'Publication', 'GenomicContext', 'Gene',
                  'Locus']
association_cache = TTLCache(maxsize=2 * 1024, ttl=60 * 60 * 24)
study_cache = TTLCache(maxsize=2 * 1024, ttl=60 * 60 * 24)
trait_cache = TTLCache(maxsize=2 * 1024, ttl=60 * 60 * 24)
snp_cache = TTLCache(maxsize=2 * 1024, ttl=60 * 60 * 24)
ancestry_cache = TTLCache(maxsize=2 * 1024, ttl=60 * 60 * 24)
publication_cache = TTLCache(maxsize=2 * 1024, ttl=60 * 60 * 24)
genomiccontext_cache = TTLCache(maxsize=2 * 1024, ttl=60 * 60 * 24)
gene_cache = TTLCache(maxsize=2 * 1024, ttl=60 * 60 * 24)
locus_cache = TTLCache(maxsize=2 * 1024, ttl=60 * 60 * 24)


def get_association(url: str, cached=True) -> List[Dict]:
    return get_data(url, cache_impl=association_cache, cached=cached, type='associations')


def get_study(url: str, cached=True) -> List[Dict]:
    return get_data(url, cache_impl=study_cache, cached=cached, type='studies')


def get_trait(url: str, cached=True) -> List[Dict]:
    return get_data(url, cache_impl=trait_cache, cached=cached, type='efo_traits')


def get_snp(url: str, cached=True) -> List[Dict]:
    return get_data(url, cache_impl=snp_cache, cached=cached, type='snps')


def get_ancestry(url: str, cached=True) -> List[Dict]:
    return get_data(url, cache_impl=ancestry_cache, cached=cached, type='ancestries')


def get_publication(url: str, cached=True) -> List[Dict]:
    return get_data(url, cache_impl=publication_cache, cached=cached, type='publications')


def get_genomicContext(url: str, cached=True) -> List[Dict]:
    return get_data(url, cache_impl=genomiccontext_cache, cached=cached, type='genomic_contexts')


def get_gene(url: str, cached=True) -> List[Dict]:
    return get_data(url, cache_impl=gene_cache, cached=cached, type='genes')


def get_locus(url: str, cached=True) -> List[Dict]:
    return get_data(url, cache_impl=locus_cache, cached=cached, type='loci')


def reinit_cache(field: str = 'All', maxsize: int = 2 * 1024, ttl: int = 60 * 60 * 24) -> None:
    """
    Reinitialize some or all of the cache pools.

    Args:
        field: It can be one of the following: 'All', 'Association', 'Study', 'Trait', 'Snp', 'Ancestry', 'Publication', 'GenomicContext', 'Gene','Locus'
        maxsize: The number of HTTP requests that can be stored in the cache pool.
        ttl: The unit is seconds. Expiration time.

    Returns:
        None

    ```Python
    from pandasgwas import get_publications, reinit_cache

    # Reinitialize all the cache pools. Each cache pool can only store 10 HTTP requests, and the expiration time is 10 seconds.
    reinit_cache(field = 'All', maxsize = 10, ttl = 10)
    pub = get_publications()
    # Reinitialize the cache used by get_publications()
    reinit_cache(field = 'Publication', maxsize = 10, ttl = 10)
    pub = get_publications()
    ```
    """
    global association_cache
    global study_cache
    global trait_cache
    global snp_cache
    global ancestry_cache
    global publication_cache
    global genomiccontext_cache
    global gene_cache
    global locus_cache
    if field not in fields_and_all:
        raise Exception('The field must one of %s' % str(fields_and_all))
    if field == "Association":
        association_cache = TTLCache(maxsize=maxsize, ttl=ttl)
    elif field == "Study":
        study_cache = TTLCache(maxsize=maxsize, ttl=ttl)
    elif field == "Trait":
        trait_cache = TTLCache(maxsize=maxsize, ttl=ttl)
    elif field == "Snp":
        snp_cache = TTLCache(maxsize=maxsize, ttl=ttl)
    elif field == "Ancestry":
        ancestry_cache = TTLCache(maxsize=maxsize, ttl=ttl)
    elif field == "Publication":
        publication_cache = TTLCache(maxsize=maxsize, ttl=ttl)
    elif field == "GenomicContext":
        genomiccontext_cache = TTLCache(maxsize=maxsize, ttl=ttl)
    elif field == "Gene":
        gene_cache = TTLCache(maxsize=maxsize, ttl=ttl)
    elif field == "Locus":
        locus_cache = TTLCache(maxsize=maxsize, ttl=ttl)
    else:
        association_cache = TTLCache(maxsize=maxsize, ttl=ttl)
        study_cache = TTLCache(maxsize=maxsize, ttl=ttl)
        trait_cache = TTLCache(maxsize=maxsize, ttl=ttl)
        snp_cache = TTLCache(maxsize=maxsize, ttl=ttl)
        ancestry_cache = TTLCache(maxsize=maxsize, ttl=ttl)
        publication_cache = TTLCache(maxsize=maxsize, ttl=ttl)
        genomiccontext_cache = TTLCache(maxsize=maxsize, ttl=ttl)
        gene_cache = TTLCache(maxsize=maxsize, ttl=ttl)
        locus_cache = TTLCache(maxsize=maxsize, ttl=ttl)


def clear_cache(field: str = 'All') -> None:
    """
    Clear some or all of the cache pools.

    Args:
        field: It can be one of the following: 'All', 'Association', 'Study', 'Trait', 'Snp', 'Ancestry', 'Publication', 'GenomicContext', 'Gene', 'Locus'

    Returns:
        None

    ```Python
    from pandasgwas.get_publication import get_publications
    from pandasgwas.client import clear_cache

    # Clear all caches.
    clear_cache('All')
    pub = get_publications()
    # Clear the cache used by get_publications()
    clear_cache('Publication')
    pub = get_publications()
    ```

    """
    if field not in fields_and_all:
        raise Exception('The field must one of %s' % str(fields_and_all))
    if field == 'All':
        for s in fields:
            eval(s.lower() + '_cache.clear()')
    else:
        eval(field.lower() + '_cache.clear()')


def get_data(url: str, cache_impl=None, cached=True, type=None) -> List[Dict]:
    if cache_impl is None:
        cache_impl = {}

    def _clean_item(item: dict, type: str) -> dict:
        href = item.get("_links", {}).get("self", {}).get("href")

        if type == "ancestries" and href:
            match = re.search(r"/studies/([^/]+)/ancestries/([^/?#]+)", href)
            if match:
                item["accession_id"] = match.group(1)
                item["ancestry_id"] = match.group(2)

        elif type == "genomic_contexts" and href:
            match = re.search(
                r"/single-nucleotide-polymorphisms/([^/]+)/genomic-contexts/([^/?#]+)",
                href
            )
            if match:
                item["rs_id"] = match.group(1)
                item["genomicContext_id"] = match.group(2)

        elif type == "loci" and href:
            match = re.search(r"/associations/([^/]+)/loci/([^/?#]+)", href)
            if match:
                item["association_id"] = match.group(1)
                item["locus_id"] = match.group(2)

        item.pop("links", None)
        item.pop("_links", None)

        return item

    with requests.Session() as s:
        s.mount('https://', HTTPAdapter(max_retries=5))
        if url in cache_impl and cached:
            r = cache_impl[url]
        else:
            r = s.get(url)

        if r.status_code == 200:
            parsed_data = json.loads(r.text)
            if len(parsed_data) == 0:
                return []
            cache_impl[url] = r

            if parsed_data.get('_embedded') is not None:
                results_list = parsed_data.get('_embedded').get(type, [])
                if parsed_data.get('_links') is not None and parsed_data.get('page') is not None:
                    bar = progressbar.ProgressBar(
                        max_value=parsed_data.get('page').get('totalElements')
                    ).start()
                    progress = len(results_list)
                    bar.update(progress)
                    next_link = parsed_data.get('_links').get('next')
                    while next_link is not None:
                        next_url = next_link.get('href')
                        if (next_url in cache_impl) and cached:
                            r = cache_impl[next_url]
                        else:
                            r = s.get(next_url)
                            if r.status_code == 200:
                                cache_impl[next_url] = r
                            elif r.status_code == 404:
                                return []
                            else:
                                raise Exception(
                                    'The request for %s failed: response code was %d, content was \n%s'
                                    % (next_url, r.status_code, r.text)
                                )
                        parsed_data = json.loads(r.text)
                        current_items = parsed_data.get('_embedded', {}).get(type, [])
                        results_list.extend(current_items)
                        progress += len(current_items)
                        bar.update(progress)
                        next_link = parsed_data.get('_links', {}).get('next')
                    bar.finish()
                for item in results_list:
                    if isinstance(item, dict):
                        _clean_item(item, type)
                return results_list
            else:
                if parsed_data.get('page') is not None:
                    return []
                _clean_item(parsed_data, type)
                parsed_data.pop('_links', None)
                return [parsed_data]
        elif r.status_code == 404:
            return []
        else:
            raise Exception(
                'The request for %s failed: response code was %d, content was \n%s'
                % (url, r.status_code, r.text)
            )

def ask_yes_no_question(question: str) -> str:
    yes_no_answer = ""
    while yes_no_answer != "YES" and yes_no_answer != "NO":
        yes_no_answer = input(question).upper()
    return yes_no_answer
