from convert_gene_names import get_gene


def test_convert_all_genes():
    genes = ['XM_00001', "ENST00000434970.2", "HBB"]

    real_name = get_gene(genes[0])
    orn  =get_gene(genes[2])
    assert real_name is None

