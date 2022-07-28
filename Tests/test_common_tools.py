import pytest
import requests
import numpy as np
from bs4 import BeautifulSoup

import common_tools as ct
from common_tools.common_tools import get_list_from_file

short_csv = [{'database': 'Antagomirbase', 'year': '2011', 'Organism': 'Unspecified',
              'Website': 'http://bioinfopresidencycollegekolkata.edu.in/antagomirs.html', 'Download': 'No',
              'Available': 'No',
              'PubmedID': '21904438'},
             {'database': 'ARN', 'year': '2016', 'Organism': 'Unspecified', 'Website': 'http://210.27.80.93/arn/',
              'Download': 'No', 'Available': 'Yes', 'PubmedID': '27503118'}]


def test_get_list_from_file():
    databases = get_list_from_file(file_name="../../resources/data/human_only_databases_list.txt")
    assert isinstance(databases, list)


def test_get_paper_data(monkeypatch):
    """
    This test will make sure to have all the 8 fields of the datasets;
    'database', 'year', 'Organism', 'Website', 'Download', 'Available', 'PubmedID', and 'cite_number'
    :param monkeypatch:
    :return:
    """
    # TODO Change the categories and change the name from db to something more general
    # helper.category_list_address = "../../paper_helper/resources/data/db_categories.yml"

    # monkeypatch.setattr(common_tools, "get_csv_into_dictionary", lambda *args, **kwargs: short_csv)
    # dics = helper.get_paper_data_from_file(file_name="../test_resources/papers_data.csv")
    # assert isinstance(dics, list)
    # dic = dics[0]
    # assert len(dic) == 8
    pass


def test_get_yml_from_file():
    # TODO same here
    # dic = yml_to_dict(file_name="../../resources/data/db_categories.yml")
    # assert isinstance(dic, dict)
    # assert dic['Regulation network']
    pass


def test_write_list_of_dict():
    ct.write_list_of_dict(short_csv)


def test_write_dictionary_to_tsv():
    ct.write_dictionary_to_tsv(file_name="x.tsv", keys=short_csv[0].keys(), list_dic=short_csv)


def test_merge_article_files_by():
    ct.merge_article_files_by(files_2_merge=["../../resources/data/pareto_fronts_databases/papers_data.csv",
                                             "../../resources/data/pareto_fronts_databases/pubmed_results_databses.csv",
                                             "../../resources/data/pareto_fronts_databases/pubmed_results_mirna.csv"])
    pass
