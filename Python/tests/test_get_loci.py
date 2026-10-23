from pandasgwas import get_loci


def test_get_loci():
    locus_a = get_loci(association_id='5860')
    locus_a = get_loci(association_id='169106883', mode='Thin')
    locus_a = get_loci(association_id='169106883')

    assert len(locus_a) >= 1
    assert locus_a.loci.shape[0] >= 1
    assert locus_a.loci.iloc[0]['association_id'] == '169106883'

    locus_b = get_loci(association_id='169106883', locus_id='169106882')

    assert len(locus_b) == 1
    assert locus_b.loci.shape[0] == 1
    assert locus_b.loci.iloc[0]['association_id'] == '169106883'
    assert locus_b.loci.iloc[0]['locus_id'] == '169106882'
    assert locus_b.strongest_risk_alleles.shape[0] >= 1

    assert locus_b[0] == locus_b
    assert locus_b['169106883|169106882'] == locus_b
