from pandasgwas import get_ancestries


def test_get_ancestries():
    ancestry_a = get_ancestries(accession_id='GCST000854', mode='Thin')
    ancestry_a = get_ancestries(accession_id='GCST000854')

    assert len(ancestry_a) >= 1
    assert ancestry_a.ancestries.shape[0] >= 1
    assert ancestry_a.ancestries.iloc[0]['accession_id'] == 'GCST000854'

    ancestry_b = get_ancestries(accession_id='GCST000853', ancestry_id='7529')

    assert len(ancestry_b) == 1
    assert ancestry_b.ancestries.shape[0] == 1
    assert ancestry_b.ancestries.iloc[0]['accession_id'] == 'GCST000853'
    assert ancestry_b.ancestries.iloc[0]['ancestry_id'] == '7529'

    assert ancestry_b[0] == ancestry_b
    assert ancestry_b['GCST000853|7529'] == ancestry_b
