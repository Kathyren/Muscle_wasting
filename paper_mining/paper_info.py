# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import getopt
import sys
from pprint import pprint

import eutilities as eutilities



def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


def get_papers_from_NCBI(search, retmax=100000):
    """
    This function will get term to see for and give all the results as find them
    :param retmax: How many results are you expecting, the max possible by eutils api is 100,000
    :param search:
    :param get_first:
    :return:
    """
    connection = eutilities.EutilsConnection(eutilities.NCBIDatabases.Pubmed)
    id = connection.fetch_queries_ids(term=search, retmax=retmax, db="pubmed")
    result = connection.get_papers_information(db_id=id, db="pubmed")
    return result
    pass

