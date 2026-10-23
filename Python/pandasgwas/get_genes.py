from pandasgwas.client import get_gene
from pandasgwas.Gene import Gene


BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/genes"


def get_genes(
    gene_name: str = None,
    cached: bool = True,
    mode: str = 'Fat'
) -> Gene:
    """
    Get Gene data from the server.

    Args:
        gene_name: Gene Name.
        cached: Whether or not to try to get data from the cache.
        mode: Fat or Thin. Specifies the mode of the returned object.


    Returns:
        A Gene object. Attributes of type DataFrame have hierarchical dependencies.

    ```Python
    from pandasGWAS2.get_genes import get_genes

    ch = get_genes(gene_name='HBS1L')
    ```
    """
    query_url = BASE_URL if gene_name is None else f"{BASE_URL}/{gene_name}"

    return Gene(
        data=get_gene(query_url, cached=cached),
        mode=mode
    )
