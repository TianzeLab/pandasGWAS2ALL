from pandasgwas import get_publications


def test_get_publications():
    publication_a = get_publications(pubmed_id='35241825', mode='Thin')
    publication_a = get_publications(pubmed_id='35241825')

    assert len(publication_a) == 1
    assert publication_a.publications.shape[0] == 1
    assert publication_a.publications.iloc[0]['pubmed_id'] == '35241825'

    publication_b = get_publications(title='Genome-wide association study')
    assert len(publication_b) >= 1

    publication_a = get_publications(pubmed_id='35241825')
    publication_b = get_publications(title='Genome-wide association study')

    assert publication_a[0] == publication_a
    assert publication_a['35241825'] == publication_a

    assert len(publication_b[0:2]) == 2
    assert len(publication_b[range(2)]) == 2
    assert publication_b[0] == publication_b[publication_b.raw_data[0]['pubmed_id']]

    assert len(publication_a + publication_b[0]) >= 1
    assert len(publication_b - publication_b[0]) == len(publication_b) - 1
    assert len(publication_b & publication_b[0]) == 1
    assert len(publication_b | publication_b[0]) == len(publication_b)
    assert len(publication_b ^ publication_b[0]) == len(publication_b) - 1

    publication_d = get_publications(pubmed_id='35241825')
    assert publication_d == publication_a

