from urllib.parse import quote

from pandasgwas.client import get_genomicContext
from pandasgwas.GenomicContext import GenomicContext


BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/single-nucleotide-polymorphisms"


def get_genomic_contexts(
    rs_id: str = None,
    genomicContext_id: str = None,
    cached: bool = True,
    mode: str = 'Fat'
) -> GenomicContext:
    """
    Get GenomicContext data from the server.

    Args:
        rs_id: SNP rsID; if a haplotype it may include more than one rs number (multiple SNPs comprising the haplotype).
        genomicContext_id: primary identifier of genomicContext table.
        cached: Whether or not to try to get data from the cache.
        mode: Fat or Thin. Specifies the mode of the returned object.


    Returns:
        A GenomicContext object. Attributes of type DataFrame have hierarchical dependencies.

    ```Python
    from pandasGWAS2.get_genomic_contexts import get_genomic_contexts

    ch = get_genomic_contexts(rs_id='rs3093017')
    ```
    """
    if rs_id is None or str(rs_id).strip() == "":
        raise Exception("rs_id cannot be empty")

    query_url = f"{BASE_URL}/{quote(str(rs_id), safe='')}/genomic-contexts"
    if genomicContext_id is not None and str(genomicContext_id).strip() != "":
        query_url = f"{query_url}/{quote(str(genomicContext_id), safe='')}"

    return GenomicContext(
        data=get_genomicContext(query_url, cached=cached),
        mode=mode
    )

