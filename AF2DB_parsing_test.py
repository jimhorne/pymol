# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 11:17:12 2024

@author: james

Takes around 0.5 - 0.9s per entry
"""

import time
start_time = time.time()

import requests
from urllib.parse import urlsplit
import os
import pandas as pd

def fetch_info(UniprotID):
    """ Accesses AF2 predictions from EBI AF2 DB via a GET request to 
    AF2 DB API
    """
    base_url = 'https://alphafold.ebi.ac.uk/api/prediction/'
    url = f'{base_url}{UniprotID}'
    
    print(f'Processing {UniprotID}')
    
    try:
        response = requests.get(url)
        if response.status_code == 200:
            return response
        elif response.status_code == 404:
            print('Error 404: Resource not found')
            return None
        else:
            print(f'Error unexpected response: {response.status_code}')
            return None
    except requests.exceptions.RequestException as e:
        # This will catch network-related errors
        print(f'Error occurred: {e}')
        return None


def download_file(url):
    """ Download a file from a URL 
    """
    # extract filename from base url - this will save the file into the PWD of the script unless absolute path appended
    filename = os.path.basename(urlsplit(url).path)
    filename = os.path.join(r'PAP_pdbs', filename)
    
    # Check if a file already exists first before trying to download again
    if os.path.exists(filename):
        print(f'File {filename} already exists. Skipping download.')
        return True
    
    try:
        # send GET request to url
        response = requests.get(url, stream=True)
        if response.status_code == 200:
            with open(filename, 'wb') as file:
                for chunk in response.iter_content(chunk_size=8192): # file written in chunks specified in bytes to be more memory efficient
                    file.write(chunk)
            print(f'File saved as {filename}')
            return True
        else:
            print(f'Error: failed to download file. Status code {response.status_code}')
            return False
    
    except requests.exceptions.RequestException as e:
        print(f'Error occurred: {e}')
        return False
                

def load_spreadsheet(file_path, xlsx_file=False):
    '''Load data in a spreadsheet format (either TSV or XLSX) into memory as
    a pandas dataframe. Expects a header row
    '''
    try:
        if xlsx_file:
            df = pd.read_excel(file_path, header=0)
        else:
            df = pd.read_csv(file_path, sep='\t', header=0)
        
        print(f'Check dataframe headers loaded in are correct: {df.columns.tolist()}')
        print(f'Check dataframe size is correct (excl. header row): {df.shape}')
        return df
    
    except Exception as e:
        print(f"Failed to process file {file_path}: {e}")


def retrieve_AF2db_PDB(df):
    """ Take a dataframe containing uniprot accessions in first column (0th col) and 
    retrieve AF2DB URL for predicted structure in PDB format. Download to local dir.
    Assumes first column contains the Uniprot accessions
    Args:
        df (pandas dataframe): dataframe where the first column contains your uniprot
        accession numbers
    """
    failed_entry_collector = [] # log missing entries (most commonly get 404 errors where not present in AF2DB)
    metadata_collector = []
    count = 0
    
    for i, row in df.iterrows():
        temp_dict = {}
        acc_value = row.iloc[0] # assumes accessions are stored in first column
        af2db_response = fetch_info(acc_value)
        if af2db_response is not None: # will only continue when status code was 200 and returned data
            data = af2db_response.json()[0]
        else:
            failed_entry_collector.append(acc_value)
            continue
        url = data['pdbUrl']
        success = download_file(url)
        count += 1
        
        # Create metadata dictionary for entry
        # Define keys to copy
        keys_to_copy = [
            'entryId', 'gene', 'sequenceVersionDate', 'uniprotAccession', 'uniprotId', 
            'uniprotDescription', 'taxId', 'organismScientificName', 'uniprotSequence', 
            'modelCreatedDate', 'paeImageUrl', 'paeDocUrl'
            ]
        
        # Copy values from data to temp_dict
        temp_dict = {key: data[key] for key in keys_to_copy} # dict comprehension
        temp_dict['downloaded'] = 'Y' if success else 'N'
        
        # Append dictionary to master metadata list
        metadata_collector.append(temp_dict)
        
    # Save metadata to output file
    export_dict_list_to_spreadsheet(metadata_collector, 'AF2db_PAP_metadata') # supply filename as string
    print(f'Retrieved URLs for {count}/{len(df)} entries.')
    print(f'The following entries were not retrieved: {failed_entry_collector}. See output print log for details.')
    

def export_dict_list_to_spreadsheet(dict_list, filename):
    """
    Converts a list of dictionaries into a spreadsheet with the keys as headers and each
    dictionary as a row.
    Args:
        dict_list (list): a list object containing dictionaries
        filename (str): the filename to save as (excluding extension)
    """
    df = pd.DataFrame(dict_list)
    df.to_csv(f'{filename}.tsv', sep='\t', index=False)


database_dir = r'C:\Users\james\OneDrive\Cambridge Postdoc Project Docs\Projects\YbjP (MacAB-TolC)\Bioinformatics\InterPro\PAPs\CD-Hit_PAP_all\50pc_reclust'
database_file = 'output.tsv'
full_path = os.path.join(database_dir, database_file)
#full_path = 'tolC_PF02321_uni90_nofrag_EcTolC_55pc-above-cluster_onlyacc_trunc.csv'

df = load_spreadsheet(full_path)

retrieve_AF2db_PDB(df)


end_time = time.time()
execution_time = end_time - start_time
print('\nExecution time:', execution_time, "seconds")