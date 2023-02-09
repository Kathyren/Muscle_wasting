import csv
import bios
import requests
from bs4 import BeautifulSoup


def write_dictionary_to_tsv(file_name="x.csv", keys=[], list_dic=None):
    if list_dic is None:
        dic = [{}]
    with open(file_name, mode='w', newline='', encoding="utf8") as outfile:
        dict_writer = csv.DictWriter(outfile, fieldnames=keys, delimiter='\t')
        dict_writer.writeheader()
        dict_writer.writerows(list_dic)


def get_list_from_file(file_name="resources/data/human_only_databases_list.txt"):
    with open(file_name, 'r') as file:
        file = file.read()
        file_list = file.split("\n")
    return file_list


def get_tsv_into_dictionary(file_name="resources/data/human_only_databases_list.txt"):
    with open(file_name, 'r') as f:
        dic = [{k: str(v) for k, v in row.items()}
               for row in csv.DictReader(f, skipinitialspace=True, delimiter='\t')]
    return dic


def get_csv_into_dictionary(file_name="resources/data/human_only_databases_list.txt"):
    with open(file_name, 'r', encoding="utf8") as f:
        dic = [{k: str(v) for k, v in row.items()}
               for row in csv.DictReader(f, skipinitialspace=True)]
    return dic


def write_dictionary_to_csv(file_name="x.csv", dic=None):
    if dic is None:
        dic = {}
    with open(file_name, mode='w') as outfile:
        writer = csv.writer(outfile)
        mydict = {rows[0]: rows[1] for rows in dic}


def yml_to_dict(
        file_name=r"C:\Users\crtuser\Documents\PhD\Project\repos\miRNA_small_tools\paper_helper\resources\data\db_categories.yml"):
    my_dict = bios.read(file_name)
    return my_dict


def get_soup_from_html(url):
    req = requests.get(url)
    soup = BeautifulSoup(req.content, 'html.parser')
    return soup


def write_list_of_dict(list_dict, file_name="datasets_pareto_front.csv"):
    with open(file_name, 'w') as f:
        header = []
        for key in list_dict[0]:
            header.append(key)
            string_value = str(header).replace("[", "").replace("]", "").replace("'", "\"") + "\n"
        f.write(string_value)
        for element in list_dict:
            line = []
            for key, value in element.items():
                if isinstance(value, str):
                    line.append(value.replace(",", "-"))
                else:
                    line.append(str(value))
            string_value = str(line).replace("[", "").replace("]", "").replace("'", "\"") + "\n"
            f.write(string_value)
