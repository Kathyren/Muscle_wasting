### This will run the pathway analysis. It will ask as parameter the path to the yml
# As optional, a string with the path of the pathways,
# and the pvalue=0.05, threshold=None, interest='stat', pathway_pvalue=0.05.
import sys
import os.path

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
os.chdir(BASE_DIR)
pwd = os.getcwd()
import argparse
import Pathway


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run the Pathways analysis")
    parser.add_argument('--dds_files', type=str, help='YAML file with the names of the comparisons and the path to the DDS file')
    parser.add_argument('--pathways_dbs', type=str, help='The file path with the name of the pathway dabases (collection) to use.'\
                         'For example "go_molecular_function"', default=None)
    parser.add_argument('--pvalue', type=float, help="pvalue to determine when a pathway is significantly enrriched.", default=0.05)
    parser.add_argument('--threshold', type=float, help="pvalue to determine when a pathwway is enriched.", default=None)
    parser.add_argument('--pathway_pvalue', type=float, help="pvalue to determine when a pathway is significantly enrriched.", default=0.05)
    parser.add_argument('--interest', type=str, help="What is going to be the base of desicion. log2FoldChange or stat.", default='stat')
    parser.add_argument('--save_name', type=str, help='Path to save the results')
    parser.add_argument('--save_individual', type=str, help='Path to save the results of each comparison individually', default=None)
    parser.add_argument('--top', type=int, help='number of top elements to show', default=None)

    args = parser.parse_args()
    if args.pvalue<=0 or args.pvalue>1:
            raise ValueError('Pvalue must be >0 and <=1')
    if args.threshold and args.threshold<0:
            raise ValueError('Threshold must be >=0')
    if args.interest not in ['stat', 'log2FoldChange']:
            raise ValueError('Interest must be in ["stat", "log2FoldChange"]')
    if args.pathway_pvalue<=0 or args.pathway_pvalue>1:
            raise ValueError('Pvalue must be >0 and <=1')
    pathways = Pathway.Pathways(dds_dict=args.dds_files, pvalue=args.pvalue,
                                     threshold=args.threshold, interest=args.interest, 
                                     pathway_pvalue= args.pathway_pvalue)
    if args.pathways_dbs:
           pathways.load_sel_df_from_file(args.pathways_dbs)

    pathways.save_enrriched_pathways_with_genes_ORA(args.save_name)

    if args.save_individual is not None:
        if args.top is None:
            top = 10
        else:
             top = args.top
        pathways.save_individual_pathway_results(args.save_individual, n_top=top, plots=True)
    print('Done')
