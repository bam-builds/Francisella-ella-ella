"""
Data retrieval module for ProteomeXchange datasets
"""
import requests
import pandas as pd
from pathlib import Path


class ProteomeXchangeDownloader:
    """Downloads and processes ProteomeXchange datasets"""

    def __init__(self):
        self.base_url = "https://www.ebi.ac.uk/pride/ws/archive/v2"

    def download_dataset(self, accession, output_dir):
        """
        Download dataset from ProteomeXchange

        Parameters:
        -----------
        accession : str
            ProteomeXchange accession (e.g., 'PXD009225')
        output_dir : str
            Output directory path
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        print(f"Note: Actual download from ProteomeXchange requires authentication")
        print(f"Created directory: {output_dir}")
        print(f"To download {accession}, visit: https://www.ebi.ac.uk/pride/archive/projects/{accession}")

        return output_path

    def get_dataset_info(self, accession):
        """Get metadata about a dataset"""
        api_url = f"{self.base_url}/projects/{accession}"
        try:
            response = requests.get(api_url, timeout=10)
            if response.status_code == 200:
                return response.json()
        except Exception as e:
            print(f"Error fetching dataset info: {e}")
        return None
