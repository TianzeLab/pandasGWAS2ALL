from urllib.parse import urlencode

from pandasgwas.client import get_snp
from pandasgwas.Snp import Snp


BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphisms"


def get_snps(
    rs_id: str = None,
    bp_location: int = None,
    bp_start: int = None,
    bp_end: int = None,
    pubmed_id: str = None,
    chromosome: str = None,
    mapped_gene: str = None,
    extended_geneset: bool = None,
    cached: bool = True,
    mode: str = 'Fat'
) -> Snp:
    """
    Get Snp data from the server.

    Args:
        rs_id: SNP rsId; if a haplotype it may include more than one rs number (multiple SNPs comprising the haplotype).
        bp_location: search for SNPs via their base pair location.
        bp_start: Search for SNPs on a certain chromosome and within a given base pair location range.
        bp_end: Search for SNPs on a certain chromosome and within a given base pair location range.
        pubmed_id: pubmedId of the publication.
        chromosome: Search for SNPs on a certain chromosome.
        mapped_gene: Gene(s) overlapping the variant. If a variant is intergenic, the closest 5' and 3' genes are listed.
        extended_geneset: Show extended matching genes in addition to the mapped genes.
        cached: Whether or not to try to get data from the cache.
        mode: Fat or Thin. Specifies the mode of the returned object.


    Returns:
        A Snp object. Attributes of type DataFrame have hierarchical dependencies.

    ```Python
    from pandasGWAS2.get_snps import get_snps

    ch = get_snps(rs_id='rs3093017')
    ```
    """
    params = {
        "rs_id": rs_id,
        "bp_location": bp_location,
        "bp_start": bp_start,
        "bp_end": bp_end,
        "pubmed_id": pubmed_id,
        "chromosome": chromosome,
        "mapped_gene": mapped_gene,
        "extended_geneset": extended_geneset,
    }

    params = {k: v for k, v in params.items() if v is not None}
    query_url = BASE_URL if not params else f"{BASE_URL}?{urlencode(params)}"

    return Snp(
        data=get_snp(query_url, cached=cached),
        mode=mode
    )
