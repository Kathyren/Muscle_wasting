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
    args = parser.parse_args()
    dds = DDS.DDS(dds_files_path=args.dds_files)
    dds.combine_dds_interest_genes()
    dds.save_dds_df(path=args.save_name)
    print('Done')
