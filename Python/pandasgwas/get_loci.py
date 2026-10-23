from urllib.parse import quote

from pandasgwas.client import get_locus
from pandasgwas.Locus import Locus


BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"


def get_loci(
    association_id: str = None,
    locus_id: str = None,
    cached: bool = True,
    mode: str = 'Fat'
) -> Locus:
    """
    Get Locus data from the server.

    Args:
        association_id: primary identifier of Association table.
        locus_id: primary identifier of Locus table.
        cached: Whether or not to try to get data from the cache.
        mode: Fat or Thin. Specifies the mode of the returned object.


    Returns:
        A Locus object. Attributes of type DataFrame have hierarchical dependencies.

    ```Python
    from pandasGWAS2.get_loci import get_loci

    ch = get_loci(association_id='169106883')
    ```
    """
    if association_id is None or str(association_id).strip() == "":
        raise Exception("association_id cannot be empty")

    query_url = f"{BASE_URL}/{quote(str(association_id), safe='')}/loci"
    if locus_id is not None and str(locus_id).strip() != "":
        query_url = f"{query_url}/{quote(str(locus_id), safe='')}"

    return Locus(
        data=get_locus(query_url, cached=cached),
        mode=mode
    )
