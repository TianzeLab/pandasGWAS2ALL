from pandasgwas import get_studies


def test_get_studies():
    study_a = get_studies(accession_id='GCST000854', mode='Thin')
    study_a = get_studies(accession_id='GCST000854')

    assert len(study_a) == 1
    assert study_a.studies.shape[0] == 1
    assert study_a.studies.iloc[0]['accession_id'] == 'GCST000854'

    study_b = get_studies(pubmed_id='35241825')
    assert len(study_b) >= 1

    study_c = get_studies(disease_trait="Early-onset Parkinson's disease")
    assert len(study_c) >= 1

    study_d = get_studies(accession_id='GCST000854', full_pvalue_set=False)
    assert len(study_d) >= 1

    study_e = get_studies(efo_id='MONDO_0005130')
    assert len(study_e) >= 1


    study_f = get_studies(efo_trait='celiac disease')
    assert len(study_f) >= 1

    study_g = get_studies(cohort='BioImage')
    assert len(study_g) >= 1

    study_k = get_studies(efo_id='MONDO_0005130', show_child_trait=True)
    assert len(study_k) >= 1

    study_l = get_studies(mapped_gene='CCR6')
    assert len(study_l) >= 1

    study_m = get_studies(mapped_gene='CCR6', extended_geneset=True)
    assert len(study_m) >= 1

    study_a = get_studies(accession_id='GCST000854')
    study_b = get_studies(pubmed_id='35241825')

    assert study_a[0] == study_a
    assert study_a['GCST000854'] == study_a

    assert len(study_b[0:2]) == 2
    assert len(study_b[range(2)]) == 2
    assert study_b[0] == study_b[study_b.raw_data[0]['accession_id']]

    assert len(study_a + study_b[0]) >= 1
    assert len(study_b - study_b[0]) == len(study_b) - 1
    assert len(study_b & study_b[0]) == 1
    assert len(study_b | study_b[0]) == len(study_b)
    assert len(study_b ^ study_b[0]) == len(study_b) - 1

    study_c = get_studies(accession_id='GCST000854')
    assert study_c == study_a
