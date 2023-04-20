import argparse
import logging
import os 
from pathlib import Path


import pandas as pd
from tqdm import tqdm

from scripts.utils_v2 import *
from scripts.features import *

from scripts.supplementary import *


# TODO not important: improve create_targetscan_db() to include to_csv parts.
# TODO: improve this function. add excel stuff

def create_supplementary_files(): # improved
    """
    Creates supplementary files for the project.

    The function creates a directory for the supplementary files if it doesn't exist
    already. Then, it creates a TargetScan database, and writes it to a tab-separated
    file in the supplementary files directory. Similarly, it creates a miRBase database,
    and writes it to a tab-separated file in the same directory. Finally, it calls the
    parse_clash_data() function to parse CLASH data.

    Raises:
        FileNotFoundError: If the required file doesn't exist.
        PermissionError: If you don't have permission to access the file or directory.
        pd.errors.EmptyDataError: If the data file appears to be empty.
        Exception: If an unexpected error occurs.

    Args:
        None

    Returns:
        None
    """
    
    SUPPLEMENTARY_DIR = Path("data/supplementary_files")
    SUPPLEMENTARY_DIR.mkdir(parents=True, exist_ok=True)

    try:
        targetscan_df = create_targetscan_db()
        targetscan_path = SUPPLEMENTARY_DIR / "targetscan.tsv"
        with targetscan_path.open("w", newline="") as f:
            targetscan_df.to_csv(f, sep="\t", index=False)

        mirbase_df = create_mirbase_db()
        mirbase_path = SUPPLEMENTARY_DIR / "mirbase.tsv"
        with mirbase_path.open("w", newline="") as f:
            mirbase_df.to_csv(f, sep="\t", index=False)

        # Uncomment the following code if you want to create a ta_sps.tsv file
        # ta_sps_df = read_ta_sps_data("data/raw_supplementary_files/bartel_2011_ta_sps_data.xlsx")
        # ta_sps_path = SUPPLEMENTARY_DIR / "ta_sps.tsv"
        # with ta_sps_path.open("w", newline="") as f:
        # ta_sps_df.to_csv(f, sep="\t", index=False)

        parse_clash_data()

    except FileNotFoundError as e:
        print(f"Error: {e}. Make sure that the required file exists.")
    except PermissionError as e:
        print(f"Error: {e}. Make sure that you have permission to access the file or directory.")
    except pd.errors.EmptyDataError as e:
        print(f"Error: {e}. The data file appears to be empty.")
    except Exception as e:
        print(f"Error: {e}. An unexpected error occurred.")



def analyze_data(sequence, targetscan_df, output_file, find_non_CLASH_types):
    
    with tqdm(total=13, desc='Processing data') as pbar:

        # Step 1: Find matches
        pbar.set_description('Finding matches')
        pbar.update(1)     
        df = find_matches(sequence, targetscan_df)

        # Step 2: Write results to disk
        pbar.set_description('Writing preliminary results to disk')
        pbar.update(1)
        
        try:
            df.to_csv(output_file, index=False)
        except IOError as e:
            logging.error("An error occurred while writing preliminary results to the output file: %s", str(e))
            logging.error("Please check that the output file path is correct and that you have write permissions.")
            return

        # Step 3: Add flag column
        df["flag_column"] = None

# TODO: sort them by their priorities (CLASH II first)

        # Step 4: Find sites
        pbar.set_description('Finding 8mer sites')
        pbar.update(1)
        df = find_8mer_sites(df)

        
        pbar.set_description('Finding 7mer-a1 sites')
        pbar.update(1)
        df = find_7mer_a1_sites(df)

        
        pbar.set_description('Finding 7mer-m8 sites')
        pbar.update(1)
        df = find_7mer_m8_sites(df)

        
        pbar.set_description('Finding CLASH-II sites')
        pbar.update(1)
        df = find_CLASH_II_sites(df)

        
        pbar.set_description('Finding CLASH-III sites')
        pbar.update(1)
        df = find_CLASH_III_sites(df)

        
        pbar.set_description('Finding CLASH-IV sites')
        pbar.update(1)
        df = find_CLASH_IV_sites(df)

        
        pbar.set_description('Finding CLASH-V sites')
        pbar.update(1)
        df = find_CLASH_V_sites(df)


        # Step 5: Find optional sites
        if find_non_CLASH_types:
            pbar.set_description('Finding compensatory sites')
            pbar.update(1)
            df = find_compensatory_sites(df)

            pbar.set_description('Finding seed matches with 1 mismatch')
            pbar.update(1)
            df = find_seed_with_1_mismatch(df)

            pbar.set_description('Finding centered sites')
            pbar.update(1)
            df = find_centered_site(df)

        # keeping only the matches & dropping temporary flag column
        df = df[df['flag_column'] == 1].drop('flag_column', axis=1)
        
        # ensuring alignment_string column dtype is string
        df['alignment_string'] = df['alignment_string'].astype(str)
        
        
        try:
            logging.info('Writing the results to the output file: %s', output_file)
            df.to_csv(output_file, sep="\t", index=False)
            
            
        except IOError as e:
            logging.error("An error occurred while writing the results to the output file: %s", str(e))
            logging.error("Please check that the output file path is correct and that you have write permissions.")
            return
  



def main():

    parser = argparse.ArgumentParser(description='Process some data.')
    parser.add_argument('-i',
                        '--input',
                        type=str,
                        help='Path to input file',
                        default="bag2.fa")

    parser.add_argument('-o',
                        '--output',
                        type=str,
                        help='Path to output file',
                        default="results/results.tsv")

    parser.add_argument('--find-non-CLASH-types',
                        action='store_true',
                        help='Whether to find non-CLASH types of sites',
                        default=False)
    
    parser.add_argument('--dataset',
                        type=str,
                        help='Dataset to use (TargetScan (ts) or miRBase (mb))',
                        choices=['ts', 'mb'],
                        default='ts')
    
    parser.add_argument('--preprocess',
                        action='store_true',
                        help='Whether to preprocess the data',
                        default=False)

    args = parser.parse_args()


### calling miRNA preprocessing if necessary
    if args.preprocess:
        logging.info("Running preprocessing step...")
        create_supplementary_files()
        logging.info("Preprocessing step completed successfully. Please run the script again without the --preprocess argument.")
        return

### pick and load miRNA dataset of choice

    if args.dataset == 'ts':
        try:
            targetscan_df = pd.read_csv(
                "data/supplementary_files/targetscan.tsv", sep="\t")
        except FileNotFoundError as e:
            logging.error("TargetScan file not found, please run with --preprocess argument first: %s", str(e))
            return
    elif args.dataset == 'mb':
        try:
            mirbase_df = pd.read_csv(
                "data/supplementary_files/mirbase.tsv", sep="\t")
        except FileNotFoundError as e:
            logging.error("miRBase file not found, please run with --preprocess argument first: %s", str(e))
            return
    else:
        logging.error("Invalid dataset: %s, please specify either 'ts' for TargetScan or 'mb' for miRBase.", args.dataset)
        return


### loading the sequence 
    logging.info("Loading input sequence...")
    sequence = parse_fasta_file(args.input)
    logging.info("Input sequence loaded successfully.")


### checking if output folder exists
    if not os.path.exists('results'):
        os.makedirs('results')
        
### run
    analyze_data(sequence, targetscan_df, args.output, args.find_non_CLASH_types)




if __name__ == '__main__':
    main()
