from urllib.parse import quote

from pandasgwas.client import get_ancestry
from pandasgwas.Ancestry import Ancestry


BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/studies"


def get_ancestries(
    accession_id: str = None,
    ancestry_id: str = None,
    cached: bool = True,
    mode: str = 'Fat'
) -> Ancestry:
    """
    Get Ancestry data from the server.

    Args:
        accession_id: The study’s GWAS Catalog accession ID
        ancestry_id: primary identifier of ancestry table
        cached: Whether or not to try to get data from the cache.
        mode: Fat or Thin. Specifies the mode of the returned object.


    Returns:
        A Ancestry object. Attributes of type DataFrame have hierarchical dependencies.

    ```Python
    from pandasGWAS2.get_ancestries import get_ancestries

    ch = get_ancestries(accession_id='GCST000854')
    ```
    """
    if accession_id is None or str(accession_id).strip() == "":
        raise Exception("accession_id cannot be empty")

    query_url = f"{BASE_URL}/{quote(str(accession_id), safe='')}/ancestries"
    if ancestry_id is not None and str(ancestry_id).strip() != "":
        query_url = f"{query_url}/{quote(str(ancestry_id), safe='')}"

    return Ancestry(
        data=get_ancestry(query_url, cached=cached),
        mode=mode
    )
