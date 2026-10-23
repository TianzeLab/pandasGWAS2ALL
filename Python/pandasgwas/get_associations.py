from urllib.parse import urlencode, quote

from pandasgwas.client import get_association
from pandasgwas.Association import Association


BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/associations"


def get_associations(
    association_id: str = None,
    pubmed_id: str = None,
    rs_id: str = None,
    full_pvalue_set: bool = None,
    accession_id: str = None,
    efo_trait: str = None,
    efo_id: str = None,
    show_child_trait: bool = None,
    mapped_gene: str = None,
    extended_geneset: bool = None,
    cached: bool = True,
    mode: str = 'Fat'
) -> Association:
    """
    Get Association data from the server.

    Args:
        association_id: Primary identifier of association table.
        pubmed_id: pubmedId of the publication.
        rs_id: SNP rsID; if a haplotype it may include more than one rs number (multiple SNPs comprising the haplotype).
        full_pvalue_set: Whether full summary statistics are available for this study.
        accession_id: The study's GWAS Catalog accession ID.
        efo_trait: The trait name or label.
        efo_id: The trait URI shortform.
        show_child_trait: display entities for descendants of a parent Efo Trait Term.
        mapped_gene: Gene(s) overlapping the variant. If a variant is intergenic, the closest 5' and 3' genes are listed.
        extended_geneset: Show extended matching genes in addition to the mapped genes.
        cached: Whether or not to try to get data from the cache.
        mode: Fat or Thin. Specifies the mode of the returned object.


    Returns:
        An Association object. Attributes of type DataFrame have hierarchical dependencies.

    ```Python
    from pandasGWAS2.get_associations import get_associations

    ch = get_associations(association_id='169106883')
    ```
    """
    params = {
        "pubmed_id": pubmed_id,
        "rs_id": rs_id,
        "full_pvalue_set": full_pvalue_set,
        "accession_id": accession_id,
        "efo_trait": efo_trait,
        "efo_id": efo_id,
        "show_child_trait": show_child_trait,
        "mapped_gene": mapped_gene,
        "extended_geneset": extended_geneset,
    }

    params = {k: v for k, v in params.items() if v is not None}
    query_url = BASE_URL if not params else f"{BASE_URL}?{urlencode(params)}"

    has_id = association_id is not None
    has_query_filters = bool(params)

    result_by_id = None
    result_by_query = None

    if has_id:
        id_url = f"{BASE_URL}/{quote(association_id)}"
        result_by_id = Association(
            data=get_association(id_url, cached=cached),
            mode=mode
        )

    if not has_id or has_query_filters:
        result_by_query = Association(
            data=get_association(query_url, cached=cached),
            mode=mode
        )

    if has_id and has_query_filters:
        return result_by_id & result_by_query
    if has_id:
        return result_by_id
    return result_by_query
