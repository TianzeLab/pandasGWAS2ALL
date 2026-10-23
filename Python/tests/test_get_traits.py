from pandasgwas import get_traits


def test_get_traits():
    trait_a = get_traits(efo_id='EFO_0009713', mode='Thin')
    trait_a = get_traits(efo_id='EFO_0009713')

    assert len(trait_a) == 1
    assert trait_a.traits.shape[0] == 1
    assert trait_a.traits.iloc[0]['efo_id'] == 'EFO_0009713'

    trait_b = get_traits(pubmed_id='35241825')
    assert len(trait_b) >= 1

    trait_c = get_traits(efo_trait='celiac disease')
    assert len(trait_c) >= 1

    trait_d = get_traits(uri='http://www.ebi.ac.uk/efo/EFO_0009713')
    assert len(trait_d) >= 1

    trait_e = get_traits(mapped_gene='CCR6')
    assert len(trait_e) >= 1

    trait_f = get_traits(mapped_gene='CCR6', extended_geneset=True)
    assert len(trait_f) >= 1

    trait_a = get_traits(efo_id='EFO_0009713')
    trait_b = get_traits(pubmed_id='35241825')

    assert trait_a[0] == trait_a
    assert trait_a['EFO_0009713'] == trait_a

    assert len(trait_e[0:2]) == 2
    assert len(trait_e[range(2)]) == 2
    assert trait_e[0] == trait_e[trait_e.raw_data[0]['efo_id']]

    assert len(trait_a + trait_b[0]) >= 1
    assert len(trait_b - trait_b[0]) == len(trait_b) - 1
    assert len(trait_b & trait_b[0]) == 1
    assert len(trait_b | trait_b[0]) == len(trait_b)
    assert len(trait_b ^ trait_b[0]) == len(trait_b) - 1

    trait_c = get_traits(efo_id='EFO_0009713')
    assert trait_c == trait_a

