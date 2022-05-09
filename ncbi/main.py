# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import ncbi.eutilities as eutilities



def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


def get_papers_from_NCBI(search, get_first=False):
    """
    This function will get term to see for and give all the results as find them
    :param search:
    :param get_first:
    :return:
    """
    connection = eutilities.EutilsConnection(eutilities.NCBIDatabases.Pubmed)
    id = connection.fetch_queries_ids(term=search, get_first=get_first, db="pubmed")
    result = connection.get_ids_information(db_id=id, db="pubmed")
    return result
    pass


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    connection = eutilities.EutilsConnection(eutilities.NCBIDatabases.Nucleotides)
    id = connection.fetch_queries_ids(term='XM_537211', get_first=True)[0]
    result = connection.get_id_information(db_id=id)
    gene = result.gene
    print(gene)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
