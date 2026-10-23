from urllib.parse import urlencode

from pandasgwas.client import get_trait
from pandasgwas.Trait import Trait


BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits"


def get_traits(
    efo_trait: str = None,
    uri: str = None,
    efo_id: str = None,
    pubmed_id: str = None,
    mapped_gene: str = None,
    extended_geneset: bool = None,
    cached: bool = True,
    mode: str = 'Fat'
) -> Trait:
    """
    Get Trait data from the server.

    Args:
        efo_trait: The trait name or label.
        uri: The trait URI or unique identifier.
        efo_id: The trait URI shortform.
        pubmed_id: pubmedId of the publication.
        mapped_gene: Gene(s) overlapping the variant. If a variant is intergenic, the closest 5' and 3' genes are listed.
        extended_geneset: Show extended matching genes in addition to the mapped genes.
        cached: Whether or not to try to get data from the cache.
        mode: Fat or Thin. Specifies the mode of the returned object.


    Returns:
        A Trait object. Attributes of type DataFrame have hierarchical dependencies.

    ```Python
    from pandasGWAS2.get_traits import get_traits

    ch = get_traits(efo_id='EFO_0001360')
    ```
    """
    params = {
        "efo_trait": efo_trait,
        "uri": uri,
        "efo_id": efo_id,
        "pubmed_id": pubmed_id,
        "mapped_gene": mapped_gene,
        "extended_geneset": extended_geneset,
    }

    params = {k: v for k, v in params.items() if v is not None}
    query_url = BASE_URL if not params else f"{BASE_URL}?{urlencode(params)}"

    return Trait(
        data=get_trait(query_url, cached=cached),
        mode=mode
    )
