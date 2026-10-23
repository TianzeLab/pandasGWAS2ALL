from pandasgwas import get_genes


def test_get_genes():
    gene_a = get_genes(gene_name='HBS1L', mode='Thin')
    gene_a = get_genes(gene_name='HBS1L')

    assert len(gene_a) == 1
    assert gene_a.genes.shape[0] == 1
    assert gene_a.genes.iloc[0]['gene_name'] == 'HBS1L'

    gene_a = get_genes(gene_name='HBS1L')
    gene_b = get_genes(gene_name='SAMD5')

    assert gene_a[0] == gene_a
    assert gene_a['HBS1L'] == gene_a

    assert gene_b[0] == gene_b[gene_b.raw_data[0]['gene_name']]

    assert len(gene_a + gene_b[0]) >= 1
    assert len(gene_b - gene_b[0]) == len(gene_b) - 1
    assert len(gene_b & gene_b[0]) == 1
    assert len(gene_b | gene_b[0]) == len(gene_b)
    assert len(gene_b ^ gene_b[0]) == len(gene_b) - 1

    gene_c = get_genes(gene_name='HBS1L')
    assert gene_c == gene_a



