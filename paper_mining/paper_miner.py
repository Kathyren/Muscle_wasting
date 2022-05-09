#### Main funtion for paper_miner


import argparse
import sys, getopt
import common_tools as ct
from ncbi.main import get_papers_from_NCBI

default_output = "paper_miner_output.csv"


def paper_miner(regular_exp, min_papers, max_papers, output, print_diagram):
    """

    :param regular_exp: STRING The regular expression that is going to be look for in pubmed
    :param min_papers: INT The minimum amount of papers expected on the output
    :param max_papers: INT The max amount of papers expected from the output
    :param output: STRING Where to save the output
    :param print_diagram: BOOL If you want to see the pareto diagram
    :return: 0 if not papers where found and 1 if papers where saved
    """
    papers_dictionary = get_papers_from_NCBI(search=regular_exp)
    if len(papers_dictionary) > 0:
        save_papers("pre_process_" + output, papers_dictionary)
        return 1
    else:
        return 0


def save_papers(output_file, papers_info):
    ct.write_dictionary_to_tsv(output_file, keys=papers_info[0].keys(), list_dic=papers_info)


def main(argv):
    regular_exp = ''
    min_papers = 1
    max_papers = 0
    output = ''
    print_diagram = False
    try:
        opts, args = getopt.getopt(argv, "hs:n:x:o:p", ["search=", "min=", "max=", "output=", "print="])
    except getopt.GetoptError:
        print(f'The format for Paper miner is: ')
        print('paper_miner.py -s <search> -n <min> -x <max> -o <output> \n'
              '\t search:\tThe search that is going to be executed in pubmed as pubmed format is. '
              'ADD DOUBLE QUOTES ""\n'
              '\t min   :\tThe minimum amount of papers that are desired form the mining\n'
              '\t max   :\tThe maximum amount of papers that are desired from the mining\n'
              '\t output:\tThe name of the file (and address if it is different from current)'
              ' where the output will be saved\n'
              '\t print :\t 0 if you do not want the Pareto front graph, 1 if you want it on the'
              ' same location than the output')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('test.py -i <inputfile> -o <outputfile>')
            sys.exit()
        elif opt in ("-s", "--search"):
            regular_exp = arg
        elif opt in ("-n", "--min"):
            min_papers = arg
        elif opt in ("-x", "--max"):
            max_papers = arg
        elif opt in ("-o", "--output"):
            output = arg
        elif opt in ("-p", "--print"):
            print_diagram = arg
    if regular_exp == '':
        print("The regular expression is essential to find the papers in pubmed and it cannot be empty")
        return 0
    if output == '':
        print(f"There was not input for the output, will use the default {default_output}")
        output = default_output
    # print(regular_exp, min_papers, max_papers, output, print_diagram)
    return paper_miner(regular_exp=regular_exp, min_papers=min_papers, max_papers=max_papers,
                       output=output, print_diagram=print_diagram)


if __name__ == "__main__":
    main(sys.argv[1:])
