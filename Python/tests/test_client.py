from pandasgwas.client import get_trait

def test_get_trait():
    a_one = get_trait("https://www.ebi.ac.uk/gwas/rest/api/v2/efo-traits/EFO_0009713")
    assert len(a_one) == 1