# ### Main function for paper_miner


import sys, getopt

import common_tools as ct
from paper_info import get_papers_from_NCBI
from paper_optimization import EvaluatePapers

default_output = "paper_miner_output.csv"


def paper_miner(regular_exp, min_papers, max_papers=None, output="output.txt", print_diagram=False, input=None):
    """

    :param regular_exp: STRING The regular expression that is going to be look for in pubmed
    :param min_papers: INT The minimum amount of papers expected on the output
    :param max_papers: INT The max amount of papers from PubMed
    :param output: STRING Where to save the output
    :param print_diagram: BOOL If you want to see the pareto diagram
    :return: 0 if not papers where found and 1 if papers where saved
    """
    print("Getting the papers c:")
    if input is None:
        papers_dictionary = get_papers_from_NCBI(search=regular_exp, retmax=max_papers)
        if len(papers_dictionary) == 0:
            print(f"No articles were found with the query {regular_exp}. Please use the PubMed sintax.")
            return 0
        print(f"{len(papers_dictionary)} results found!, saving output")
        save_papers("pre_process_" + output, papers_dictionary)
    else:
        papers_dictionary = ct.get_tsv_into_dictionary(input)
        if len(papers_dictionary) == 0:
            print(f"No articles were found on {input}. Nothing else to do.")
            return 0
    if len(papers_dictionary) <= min_papers:
        print(f"No sufficient papers were found to satisfy the min "
              f"articles desired or they are the exact number. Nothing else to do.")
        return 0
    print("Evaluating best articles")
    #ep = EvaluatePapers(papers_info=papers_dictionary)
    #pareto = ep.get_pareto_cites_year(plot=True, min_papers=min_papers)
    #print(f"{len(pareto)} articles were selected!, saving output in {output}")
    #ct.write_list_of_dict(pareto, file_name=output)
    print(f"Successfully saved!!")


def save_papers(output_file, papers_info):
    ct.write_dictionary_to_tsv(output_file, keys=papers_info[0].keys(), list_dic=papers_info)


def main(argv):
    regular_exp = ''
    min_papers = 1
    max_papers = 100000
    output = ''
    print_diagram = False
    help = f'The format for Paper miner is: \n ' \
           f'paper_miner.py -s <search> -n <min> -x <max> -o <output> \n' \
           f'\t search:\tThe search that is going to be executed in pubmed as pubmed format is. ' \
           'ADD DOUBLE QUOTES ""\n' \
           '\t min   :\tThe minimum amount of papers that are desired form the mining\n' \
           '\t max   :\tThe maximum amount of papers that are desired from the mining\n' \
           '\t output:\tThe name of the file (and address if it is different from current)' \
           ' where the output will be saved\n' \
           '\t print :\t 0 if you do not want the Pareto front graph, 1 if you want it on the' \
           ' same location than the output'
    try:
        opts, args = getopt.getopt(argv, "hs:n:x:o:p", ["search=", "min=", "max=", "output=", "print="])
    except getopt.GetoptError:
        print(help)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(help)
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
    return paper_miner(regular_exp=regular_exp, min_papers=int(min_papers), max_papers=int(max_papers),
                       output=output, print_diagram=print_diagram)


def fake_main():
    regular_exp = '(mirna) AND (sarcopenia[Title])'
    regular_exp = 'HMB'

    min_papers = 1
    max_papers = 9999
    output = 'myeloma.tsv'
    print_diagram = False
    input = None # "pre_process_test_commandline.tsv"

    for y in range(1960,2000, 10):
        regular_exp = f'(Multiple Myeloma) AND (("{y}"[Date - Publication] : "{y+10}"[Date - Publication]))'
        paper_miner(regular_exp=regular_exp, min_papers=min_papers, max_papers=max_papers,
                       output=output, print_diagram=print_diagram, input=input)


if __name__ == "__main__":
    #main(sys.argv[1:])
    # print("Holi")
    fake_main()
