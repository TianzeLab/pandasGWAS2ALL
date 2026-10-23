from pandasgwas import get_genomic_contexts


def test_get_genomic_contexts():
    genomic_context_a = get_genomic_contexts(rs_id='rs3093017', mode='Thin')
    genomic_context_a = get_genomic_contexts(rs_id='rs3093017')

    assert len(genomic_context_a) >= 1
    assert genomic_context_a.genomic_contexts.shape[0] >= 1
    assert genomic_context_a.genomic_contexts.iloc[0]['rs_id'] == 'rs3093017'

    genomic_context_b = get_genomic_contexts(rs_id='rs3093017', genomicContext_id='206217817')

    assert len(genomic_context_b) == 1
    assert genomic_context_b.genomic_contexts.shape[0] == 1
    assert genomic_context_b.genomic_contexts.iloc[0]['rs_id'] == 'rs3093017'
    assert genomic_context_b.genomic_contexts.iloc[0]['genomicContext_id'] == '206217817'
    assert genomic_context_b.genes.shape[0] >= 1
    assert genomic_context_b.locations.shape[0] >= 1

    assert genomic_context_b[0] == genomic_context_b
    assert genomic_context_b['rs3093017|206217817'] == genomic_context_b
