from logger import logger
from convert_gene_names import get_gene
from ncbi import eutilities
from database_analysis import sql_operations as sql


def test_convert_all_genes():
    genes = ['XM_00001']

    real_name = get_gene(genes[0])
    assert real_name is None
