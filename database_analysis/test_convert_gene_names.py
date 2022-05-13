from convert_gene_names import get_gene


def test_convert_all_genes():
    genes = ['XM_00001']

    real_name = get_gene(genes[0])
    assert real_name is None
