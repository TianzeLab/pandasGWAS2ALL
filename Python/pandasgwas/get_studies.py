from urllib.parse import urlencode

from pandasgwas.client import get_study
from pandasgwas.Study import Study


BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2/studies"


def get_studies(
    accession_id: str = None,
    pubmed_id: str = None,
    disease_trait: str = None,
    full_pvalue_set: bool = None,
    efo_id: str = None,
    efo_trait: str = None,
    cohort: str = None,
    gxe: bool = None,
    ancestral_group: str = None,
    no_of_individuals: int = None,
    show_child_trait: bool = None,
    mapped_gene: str = None,
    extended_geneset: bool = None,
    cached: bool = True,
    mode: str = 'Fat'
) -> Study:
    """
    Get Study data from the server.

    Args:
        accession_id: The study’s GWAS Catalog accession ID.
        pubmed_id: pubmedId of the publication.
        disease_trait: Free text description of the trait investigated in this study.
        full_pvalue_set: Whether full summary statistics are available for this study.
        efo_id: The trait URI shortform.
        efo_trait: The trait name or label.
        cohort: Discovery stage cohorts used in this study.
        gxe: Whether the study investigates a gene-environment interaction.
        ancestral_group: Ancestry category group label, assigned to reduce complexity within the data sets and place samples in context.
        no_of_individuals: Number of individuals in this sample.
        show_child_trait: display entities for descendants of a parent Efo Trait Term.
        mapped_gene: Gene(s) overlapping the variant. If a variant is intergenic, the closest 5' and 3' genes are listed.
        extended_geneset: Show extended matching genes in addition to the mapped genes.
        cached: Whether or not to try to get data from the cache.
        mode: Fat or Thin. Specifies the mode of the returned object.


    Returns:
        A Study object. Attributes of type DataFrame have hierarchical dependencies.

    ```Python
    from pandasGWAS2.get_studies import get_studies

    ch = get_studies(accession_id='GCST000854')
    ```
    """
    params = {
        "accession_id": accession_id,
        "pubmed_id": pubmed_id,
        "disease_trait": disease_trait,
        "full_pvalue_set": full_pvalue_set,
        "efo_id": efo_id,
        "efo_trait": efo_trait,
        "cohort": cohort,
        "gxe": gxe,
        "ancestral_group": ancestral_group,
        "no_of_individuals": no_of_individuals,
        "show_child_trait": show_child_trait,
        "mapped_gene": mapped_gene,
        "extended_geneset": extended_geneset,
    }

    params = {k: v for k, v in params.items() if v is not None}
    query_url = BASE_URL if not params else f"{BASE_URL}?{urlencode(params)}"

    return Study(
        data=get_study(query_url, cached=cached),
        mode=mode
    )
