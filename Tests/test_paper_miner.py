import pytest
import paper_mining.paper_miner as miner
import os.path
from os import path

def test_script_execution():
    """
    This test will just make sure that the script is executable
    :return:
    """
    pass


def test_efetch_manula():
    """
    This test is to check the manual made efetch that is a copy with a sligthly modification to the
    one in eutilz.
    :return:
    """
    #miner.
    pass

def test_paper_miner():
    miner.paper_miner(regular_exp="(((mirna[Title]) AND (network[Title])))",
                      min_papers=1,
                      max_papers=1,
                      output="test.tsv",
                      print_diagram=False)
    assert path.exists("pre_process_test.tsv")
