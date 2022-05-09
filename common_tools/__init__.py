import csv
def write_dictionary_to_tsv(file_name="x.csv", keys=[], list_dic=None):
    if list_dic is None:
        dic = [{}]
    with open(file_name, mode='w', newline='') as outfile:
        dict_writer = csv.DictWriter(outfile,fieldnames=keys, delimiter='\t')
        dict_writer.writeheader()
        dict_writer.writerows(list_dic)
