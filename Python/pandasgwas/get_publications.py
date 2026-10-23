from urllib.parse import urlencode

from pandasgwas.client import get_publication
from pandasgwas.Publication import Publication


BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/publications"


def get_publications(
    pubmed_id: str = None,
    title: str = None,
    first_author: str = None,
    cached: bool = True,
    mode: str = 'Fat'
) -> Publication:
    """
    Get Publication data from the server.

    Args:
        pubmed_id: pubmedId of the publication.
        title: Title of manuscript the GWAS is included in.
        first_author: Last name and initials of first author.
        cached: Whether or not to try to get data from the cache.
        mode: Fat or Thin. Specifies the mode of the returned object.


    Returns:
        A Publication object. Attributes of type DataFrame have hierarchical dependencies.

    ```Python
    from pandasGWAS2.get_publications import get_publications

    ch = get_publications(pubmed_id='35241825')
    ```
    """
    params = {
        "pubmed_id": pubmed_id,
        "title": title,
        "first_author": first_author,
    }

    params = {k: v for k, v in params.items() if v is not None}
    query_url = BASE_URL if not params else f"{BASE_URL}?{urlencode(params)}"

    return Publication(
        data=get_publication(query_url, cached=cached),
        mode=mode
    )
