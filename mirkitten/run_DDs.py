import argparse
import sys
import os

import sys
import os.path

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
os.chdir(BASE_DIR)
pwd = os.getcwd()
import DDS

## This is to run and get the DDS. The arguments are the yml file with the name of the comparisons and the path where the DDS file is

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run DEA')
    parser.add_argument('--dds_files', type=str, help='YAML file with the names of the comparisons and the path to the DDS file')
    parser.add_argument('--save_name', type=str, help='Path to save the results')
    parser.add_argument('--interest', type=str, default='stat', help='Interest to combine the DDS files (default: stat)')
    parser.add_argument('--threshold', type=float, default=0, help='Threshold for the interest (default: 0)')
    parser.add_argument('--pvalue', type=float, default=0.05, help='P-value threshold for DE genes (default: 0.05)')
    args = parser.parse_args()
    # dds_files_path= dds_files, pvalue=0.05, threshold=0, interest='stat'
    dds = DDS.DDS(dds_files_path=args.dds_files, 
                  pvalue=args.pvalue, 
                  threshold=args.threshold, 
                  interest=args.interest)
    dds.combine_dds_interest_genes()
    dds.save_dds_df(path=args.save_name)
    print('Done')
