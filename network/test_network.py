from networkx import is_bipartite

from network.networkX import create_graph, draw_graph, load_graph, calculate_centralities
import pytest

genes = ['Cd320', 'Ndrg3', 'Aldoa', 'Bckdk', 'SLC7A1', 'ADAM17', 'NUMBL', 'FOXJ3', 'XPO6', 'AP3M2']
mirnas = ['mmu-miR-122-5p', 'mmu-miR-122-5p', 'mmu-miR-122-5p', 'mmu-miR-122-5p', 'hsa-miR-122-5p',
          'hsa-miR-122-5p', 'hsa-miR-122-5p', 'hsa-miR-122-5p', 'hsa-miR-122-5p', 'hsa-miR-122-5p']
relationship = [('Cd320', 'mmu-miR-122-5p'), ('Ndrg3', 'mmu-miR-122-5p'), ('Aldoa', 'mmu-miR-122-5p'),
                ('Bckdk', 'mmu-miR-122-5p'), ('SLC7A1', 'hsa-miR-122-5p'), ('ADAM17', 'hsa-miR-122-5p'),
                ('NUMBL', 'hsa-miR-122-5p'), ('FOXJ3', 'hsa-miR-122-5p'), ('XPO6', 'hsa-miR-122-5p'),
                ('AP3M2', 'hsa-miR-122-5p')]


def test_create_network():
    graph = create_graph(genes=genes, mirnas=mirnas, relationsip=relationship)
    assert len(graph.nodes()) == 12
    assert len(graph.edges()) == 10
    draw_graph(graph)
def test_draw():
    graph = create_graph(genes=genes, mirnas=mirnas, relationsip=relationship)
    draw_graph(graph)
    pass

def test_is_bipartite():
    graph = load_graph("small_graph.pkl")
    x = is_bipartite(graph)

def test_calculate_centralities():
    graph = load_graph("big_graph.pkl")
    calculate_centralities(graph)