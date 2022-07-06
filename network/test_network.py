from networkx import is_bipartite

from network.networkX import create_graph, draw_graph, load_graph, calculate_centralities, \
    create_graph_from_dictionaries
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


def test_create_graph_from_dictioanries():
    mirnas_dictList = [{
        "data": {
            "id": "166",
            "enhancedLabel_Passthrough": "label: attribute=\"name\" labelsize=12 labelAlignment=center outline=true outlineColor=black outlineTransparency=130 outlineWidth=5 background=false color=white dropShadow=false",
            "IntAct_AC": "EBI-10206780",
            "IntAct_Mutation": False,
            "IntAct_Taxon_Id": "9606",
            "IntAct_Preferred_Id": "P35523",
            "IntAct_Preferred_Id_Database": "uniprotkb",
            "IntAct_Species": "Homo sapiens",
            "IntAct_Type_MI_identifier": "MI:0326",
            "shared_name": "CLCN1",
            "IntAct_Type": "protein",
            "IntAct_Preferred_Id_Database_MI_identifier": "MI:0486",
            "IntAct_Identifier_Identifiers": ["EBI-10206783", "EBI-10206784", "EBI-30414706", "EBI-10206780",
                                              "EBI-10206782"],
            "IntAct_Description": "Chloride channel protein 1",
            "name": "CLCN1",
            "SUID": 166,
            "selected": False,
            "IntAct_Feature_Features": ["EBI-25146010", "EBI-25065318", "EBI-24991613", "EBI-23932890", "EBI-23509792",
                                        "EBI-23122996"]
        },
        "position": {
            "x": -48.68993826842569,
            "y": 35.43257138655283
        },
        "selected": False
    }, {
        "data": {
            "id": "164",
            "enhancedLabel_Passthrough": "label: attribute=\"name\" labelsize=12 labelAlignment=center outline=true outlineColor=black outlineTransparency=130 outlineWidth=5 background=false color=white dropShadow=false",
            "IntAct_AC": "EBI-7062247",
            "IntAct_Mutation": False,
            "IntAct_Taxon_Id": "9606",
            "IntAct_Preferred_Id": "Q9UHD4",
            "IntAct_Preferred_Id_Database": "uniprotkb",
            "IntAct_Species": "Homo sapiens",
            "IntAct_Type_MI_identifier": "MI:0326",
            "shared_name": "CIDEB",
            "IntAct_Type": "protein",
            "IntAct_Preferred_Id_Database_MI_identifier": "MI:0486",
            "IntAct_Identifier_Identifiers": ["EBI-30215058", "EBI-30215059", "EBI-30215060", "EBI-30215061",
                                              "EBI-30215062", "EBI-30215063", "EBI-7062251", "EBI-7062252",
                                              "EBI-7062253", "EBI-7062247", "EBI-7062250"],
            "IntAct_Description": "Cell death activator CIDE-B",
            "name": "CIDEB",
            "SUID": 164,
            "selected": True,
            "IntAct_Feature_Features": ["EBI-25146017", "EBI-23932897"]
        },
        "position": {
            "x": -54.39876059290216,
            "y": 145.18372932861945
        },
        "selected": True},
        {
        "data": {
            "id": "162",
            "enhancedLabel_Passthrough": "label: attribute=\"name\" labelsize=12 labelAlignment=center outline=true outlineColor=black outlineTransparency=130 outlineWidth=5 background=false color=white dropShadow=false",
            "IntAct_AC": "EBI-12859340",
            "IntAct_Mutation": False,
            "IntAct_Taxon_Id": "9606",
            "IntAct_Preferred_Id": "Q9NQX1-2",
            "IntAct_Preferred_Id_Database": "uniprotkb",
            "IntAct_Species": "Homo sapiens",
            "IntAct_Type_MI_identifier": "MI:0326",
            "shared_name": "PRDM5",
            "IntAct_Type": "protein",
            "IntAct_Preferred_Id_Database_MI_identifier": "MI:0486",
            "IntAct_Identifier_Identifiers": ["EBI-30096885", "EBI-12859340", "EBI-12859341"],
            "IntAct_Description": "PR domain zinc finger protein 5",
            "name": "PRDM5",
            "SUID": 162,
            "selected": False,
            "IntAct_Feature_Features": ["EBI-25065325", "EBI-23509799"]
        },
        "position": {
            "x": -201.89180646623885,
            "y": 33.63696071079798
        },
        "selected": False
        }]
    graph = create_graph_from_dictionaries(genes=genes, mirnas=mirnas_dictList, relationsip=relationship)
    print (list(graph.nodes))
    print(list(graph.edges))


test_create_graph_from_dictioanries()