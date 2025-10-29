# Automated dataset retrieval script
import pandas as pd
import requests
from pyteomics import mzml, mgf
import numpy as np

# Priority datasets for TCS analysis
datasets = {
    'PXD009225': 'Phosphoproteome_Fnovicida_U112',
    'PXD035145': 'Dual_proteomics_macrophage_infection',
    'PXD013074': 'OMV_stress_response',
    'PXD016669': 'Oxidative_stress_LdcF',
    'PXD001584': 'Arginine_limitation',
    'PXD022406': 'Cross_species_membrane',
    'PXD006759': 'Host_phosphoproteomics_early',
    'PXD005747': 'Dendritic_cell_SILAC',
    'PXD006908': 'Metabolic_regulator_fba',
    'PXD019739': 'GapA_interactome'
}

# Download raw files via PRIDE API
def download_proteomexchange(accession):
    api_url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/files/byProject?accession={accession}"
    response = requests.get(api_url)
    files = response.json()
    
    # Filter for raw MS files and result files
    raw_files = [f for f in files if f['fileType'] == 'RAW']
    result_files = [f for f in files if f['fileType'] in ['RESULT', 'SEARCH']]
    
    return raw_files, result_files
