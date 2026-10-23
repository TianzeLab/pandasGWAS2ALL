from pandasgwas import get_snps


def test_get_snps():
    snp_a = get_snps(rs_id='rs3093017', mode='Thin')
    snp_a = get_snps(rs_id='rs3093017')

    assert len(snp_a) == 1
    assert snp_a.snps.shape[0] == 1
    assert snp_a.snps.iloc[0]['rs_id'] == 'rs3093017'

    snp_d = get_snps(pubmed_id='35241825')
    assert len(snp_d) >= 1

    snp_e = get_snps(mapped_gene='CCR6')
    assert len(snp_e) >= 1

    snp_f = get_snps(mapped_gene='CCR6', extended_geneset=True)
    assert len(snp_f) >= 1

    snp_a = get_snps(rs_id='rs3093017')
    snp_b = get_snps(pubmed_id='35241825')

    assert snp_a[0] == snp_a
    assert snp_a['rs3093017'] == snp_a

    assert len(snp_b[0:2]) == 2
    assert len(snp_b[range(2)]) == 2
    assert snp_b[0] == snp_b[snp_b.raw_data[0]['rs_id']]

    assert len(snp_a + snp_b[0]) >= 1
    assert len(snp_b - snp_b[0]) == len(snp_b) - 1
    assert len(snp_b & snp_b[0]) == 1
    assert len(snp_b | snp_b[0]) == len(snp_b)
    assert len(snp_b ^ snp_b[0]) == len(snp_b) - 1

    snp_c = get_snps(rs_id='rs3093017')
    assert snp_c == snp_a

