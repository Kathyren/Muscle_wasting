import networkx as nx
import node_evaluation as ne
from network.networkX import create_graph, draw_graph, load_graph


def test_distance_to_target():
    graph = load_graph("test_w_mirnas_organs_n_systems.pkl")
    sources = ['muscle', 'lymphatic', 'hsa-miR-122-5p', 'hsa-miR-122-5p', 'CLCN1', 'SCN4A']
    target = 'muscle'
    paths = []
    for source in sources:
        route = ne.distance_to_target(graph, source, target)
        paths.append(route)
    assert len(paths) == len(sources)
    assert len(paths[0]) == 1, f"The source and target are the same, this should be only the target"
    test_path = paths[2]
    for i, v in enumerate(test_path[0:-1]):
        source = test_path[i]
        end = test_path[i + 1]
        edge = graph[source][end]
        assert 'weight' in edge, f"The edge only has the values {edge.keys()} but not weight."
