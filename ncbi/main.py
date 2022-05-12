# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import getopt
import sys
from pprint import pprint

import ncbi.eutilities as eutilities



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
    result = connection.get_ids_information(db_id=id, db="pubmed")
    return result
    pass

def main(argv):
    regular_exp = ''
    output = ''
    retmax = None
    help = f'The format for ncbi main searcher is: \n ' \
           f'main.py -s <search> -o <output> \n' \
           f'\t search:\tThe search that is going to be executed in pubmed as pubmed format is. ' \
           'ADD DOUBLE QUOTES ""\n' \
           '\t output:\tThe name of the file (and address if it is different from current)' \
           ' where the output will be saved\n' \

    try:
        opts, args = getopt.getopt(argv, "hs:o", ["search=", "output="])
    except getopt.GetoptError:
        print(help)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(help)
            sys.exit()
        elif opt in ("-s", "--search"):
            regular_exp = arg
        elif opt in ("-o", "--output"):
            output = arg
        elif opt in ("-n"):
            retmax = arg
    papers = get_papers_from_NCBI(search=regular_exp, retmax=retmax)
    print(f"The search {regular_exp} got {len(papers)} results")
    pprint(papers)
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main(sys.argv[1:])



# See PyCharm help at https://www.jetbrains.com/help/pycharm/
