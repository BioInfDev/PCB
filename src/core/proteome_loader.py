import os, requests

def ref_proteom(taxonomy_id='9606', extension='fasta') -> None:
    '''
    Download reviewed canonical proteins from SwissProt database.
    
    Downloads only reviewed canonical proteins for the specified taxonomy ID
    from UniProt SwissProt database and saves them in FASTA format.
    
    Args:
        taxonomy_id: NCBI taxonomy identifier for the organism (default is '9606' for human)
        extension: file extension for the proteom file (default is 'fasta')
    
    Returns:
        None
    
    Side effects:
        Saves the downloaded proteome to ./data/proteoms/ref_proteom_{taxonomy_id}.{extension}
    
    '''
    url = f'https://rest.uniprot.org/uniprotkb/stream'
    params = {
        'format': 'fasta', 
        'query': f'reviewed:true AND taxonomy_id:{taxonomy_id}'
    }
    
    response = requests.get(url=url, params=params)
    response.encoding = 'utf-8'
    
    with open(f'./data/ref_proteom_{taxonomy_id}.{extension}', 'wt') as rp:
        rp.write(response.text)

def ref_nucleolom() -> None:
    '''
    Download nucleolar proteins from Human Protein Atlas.
    
    Retrieves nucleolar protein data from Human Protein Atlas, extracts UniProt
    accessions, and downloads corresponding protein sequences from UniProt.
    Includes proteins localized to Nucleoli, Nucleoli fibrillar center, and Nucleoli rim.
    
    Returns:
        None
    
    Side effects:
        Saves the downloaded nucleolar proteome to ./data/proteoms/ref_nucleolom_9606.fasta
    
    Notes:
        - Last update: 1437 proteins (with uniprot_id)
        - Uses chunked downloads (500 proteins per request) due to UniProt API limitations
    '''
    hpa_url = 'https://www.proteinatlas.org/api/search_download.php'
    hpa_params = {
        'search': 'subcell_location:Nucleoli,Nucleoli fibrillar center,Nucleoli rim',
        'columns': 'up',
        'compress': 'no',
        'format': 'json'
    }

    hpa_response = requests.get(url=hpa_url, params=hpa_params)
    hpa_response.encoding = 'utf-8'
    
    accessions = [descriptor['Uniprot'][0] for descriptor in hpa_response.json() if descriptor['Uniprot']]
    chunk_size = 500 # Uniprot API does not allow downloading in chunks larger than 500 fasta
    with open('./data/ref_nucleolom_9606.fasta', 'wt') as rn:
        chunk = ''
        for i, acc in enumerate(accessions):
            chunk += acc + ',' 
            if i % chunk_size == 0 or i == len(accessions) - 1:
                uniprot_url = f'https://rest.uniprot.org/uniprotkb/accessions?accessions={chunk}&format=fasta'
                uniprot_response = requests.get(url=uniprot_url)
                rn.writelines(uniprot_response.text)
                chunk = ''

def fusion_proteom(rp_path):
    '''
    Fuse multiple protein sequences from FASTA file into a single continuous sequence.
    
    Args:
        rp_path: path to the input FASTA file containing multiple protein sequences
    
    Returns:
        None
    
    Side effects:
        Saves the fused sequence to a text file with '_fused.txt' suffix
        in the same directory as the input file
    '''
    with open(rp_path, 'rt') as rp:
        rp = rp.read().split('>')[1:]
    with open(f'{os.path.splitext(rp_path)[0]}_fused.txt', 'wt') as fp:
        for protein in rp:
            seq = protein.split('\n')[1:]
            fp.writelines(seq)