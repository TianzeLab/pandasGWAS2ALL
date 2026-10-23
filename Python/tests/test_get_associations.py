from pandasgwas import get_associations


def test_get_associations():
    association_a = get_associations(association_id='169106883', mode='Thin')
    association_a = get_associations(association_id='169106883')

    assert len(association_a) == 1
    assert association_a.associations.shape[0] == 1
    assert association_a.associations.iloc[0]['association_id'] == '169106883'

    association_b = get_associations(pubmed_id='35241825')
    assert len(association_b) >= 1

    association_c = get_associations(rs_id='rs3093017')
    assert len(association_c) >= 1

    association_d = get_associations(accession_id='GCST000854')
    assert len(association_d) >= 1

    association_e = get_associations(efo_trait='celiac disease')
    assert len(association_e) >= 1

    association_f = get_associations(efo_id='MONDO_0005130')
    assert len(association_f) >= 1


    association_g = get_associations(mapped_gene='CCR6')
    assert len(association_g) >= 1

    association_h = get_associations(accession_id='GCST000854', full_pvalue_set=False)
    assert len(association_h) >= 1

    association_i = get_associations(efo_trait="celiac disease", show_child_trait=True)
    assert len(association_i) >= 1
    association_j = get_associations(mapped_gene='CCR6', extended_geneset=True)
    assert len(association_j) >= 1

    association_a = get_associations(association_id='169106883')
    association_b = get_associations(accession_id='GCST000854')

    assert association_a[0] == association_a
    assert association_a['169106883'] == association_a

    assert len(association_b[0:2]) == 2
    assert len(association_b[range(2)]) == 2
    assert association_b[0] == association_b[association_b.raw_data[0]['association_id']]

    assert len(association_a + association_b[0]) >= 1
    assert len(association_b - association_b[0]) == len(association_b) - 1
    assert len(association_b & association_b[0]) == 1
    assert len(association_b | association_b[0]) == len(association_b)
    assert len(association_b ^ association_b[0]) == len(association_b) - 1

    association_c = get_associations(association_id='169106883')
    assert association_c == association_a
