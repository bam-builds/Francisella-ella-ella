#!/usr/bin/env python3
"""
Download MaxQuant proteinGroups.txt files from PRIDE for Francisella TCS analysis

This script requires network access to PRIDE repositories. Run this in an environment
with unrestricted network access (not in a sandboxed environment).
"""
import requests
import json
import os
import time
from pathlib import Path
import urllib.request
import ftplib
import sys

# Define datasets of interest
DATASETS = {
    'PXD035145': {
        'date': '2024/04',
        'description': 'Rytter 2024 - macrophage infection',
        'ftp_url': 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/04/PXD035145/'
    },
    'PXD005747': {
        'date': '2017/10',
        'description': 'Link 2018 - phosphoproteomics',
        'ftp_url': 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/10/PXD005747/'
    },
    'PXD013074': {
        'date': '2019/10',
        'description': 'Klimentova 2019 - stress conditions',
        'ftp_url': 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/10/PXD013074/'
    },
    'PXD025439': {
        'date': '2023/01',
        'description': 'Valikangas 2022 - temperature series',
        'ftp_url': 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/01/PXD025439/'
    },
    'PXD016669': {
        'date': '2020/01',
        'description': 'Felix 2021 - ldcF mutant',
        'ftp_url': 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD016669/'
    },
    'PXD001584': {
        'date': '2015/01',
        'description': 'Ramond 2015 - argP mutant',
        'ftp_url': 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/01/PXD001584/'
    },
    'PXD019739': {
        'date': '2020/09',
        'description': 'Kopeckova 2020 - GapA interactome',
        'ftp_url': 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD019739/'
    },
    'PXD022406': {
        'date': '2021/03',
        'description': 'Klimentova 2021 - cross-species',
        'ftp_url': 'https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/PXD022406/'
    }
}

def test_network_access():
    """Test if we have network access to PRIDE"""
    test_url = "https://ftp.pride.ebi.ac.uk/pride/data/archive/"
    try:
        response = requests.head(test_url, timeout=10)
        return response.status_code != 403
    except Exception as e:
        print(f"Network test failed: {e}")
        return False

def list_ftp_files_http(ftp_url):
    """Try to list files via HTTP directory listing"""
    try:
        response = requests.get(ftp_url, timeout=30)
        response.raise_for_status()

        # Parse HTML directory listing
        content = response.text
        files = []

        # Simple HTML parsing for file links
        for line in content.split('\n'):
            if 'href=' in line.lower():
                start = line.lower().find('href="') + 6
                end = line.find('"', start)
                if start > 5 and end > start:
                    filename = line[start:end]
                    if not filename.startswith('?') and not filename.startswith('/') and filename != '../':
                        files.append(filename)

        return files
    except Exception as e:
        print(f"  HTTP listing error: {e}")
        return []

def download_file_http(url, dest_path):
    """Download a file via HTTP with progress tracking"""
    try:
        response = requests.get(url, stream=True, timeout=60)
        response.raise_for_status()

        total_size = int(response.headers.get('content-length', 0))
        block_size = 8192

        print(f"  File size: {total_size / 1024 / 1024:.1f} MB")

        with open(dest_path, 'wb') as f:
            downloaded = 0
            for chunk in response.iter_content(block_size):
                if chunk:
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        print(f"\r  Downloading: {percent:.1f}%", end='', flush=True)

        print()  # New line
        return True
    except Exception as e:
        print(f"\n  Download error: {e}")
        if dest_path.exists():
            dest_path.unlink()  # Remove partial file
        return False

def generate_download_script():
    """Generate a bash script for manual downloads"""
    script_path = Path('scripts/download_pride_manual.sh')

    with open(script_path, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write('# Manual download script for PRIDE datasets\n')
        f.write('# Run this script in an environment with network access\n\n')
        f.write('mkdir -p data/pride\n\n')

        for accession, info in DATASETS.items():
            f.write(f'# {accession}: {info["description"]}\n')
            f.write(f'echo "Downloading {accession}..."\n')
            f.write(f'mkdir -p data/pride/{accession}\n')
            f.write(f'cd data/pride/{accession}\n')
            f.write(f'wget -r -np -nH --cut-dirs=5 -R "index.html*" -e robots=off \\\n')
            f.write(f'  "{info["ftp_url"]}" \\\n')
            f.write(f'  -A "*proteinGroups.txt,*peptides.txt,*evidence.txt,*summary.txt,mqpar.xml"\n')
            f.write(f'cd ../../..\n\n')

    script_path.chmod(0o755)
    print(f"\nGenerated manual download script: {script_path}")
    return script_path

def main():
    print("PRIDE MaxQuant Data Downloader")
    print("="*60)

    # Test network access
    print("\nTesting network access to PRIDE...")
    if not test_network_access():
        print("WARNING: Network access to PRIDE appears to be blocked.")
        print("This environment may have restricted network access.\n")
        print("Alternative approaches:")
        print("1. Run this script in a local environment with network access")
        print("2. Use the generated bash script for wget downloads")
        print("3. Manually download from PRIDE website\n")

        generate_download_script()

        print("\nDataset URLs for manual download:")
        for accession, info in DATASETS.items():
            print(f"\n{accession}: {info['description']}")
            print(f"  URL: {info['ftp_url']}")

        return

    # Use current directory structure
    base_dir = Path('data/pride')
    base_dir.mkdir(parents=True, exist_ok=True)

    # Summary of available files
    available_data = {}
    downloaded_files = []

    for accession, info in DATASETS.items():
        print(f"\n{'='*60}")
        print(f"{accession}: {info['description']}")
        print(f"{'='*60}")

        dataset_dir = base_dir / accession
        dataset_dir.mkdir(parents=True, exist_ok=True)

        # List files via HTTP
        files = list_ftp_files_http(info['ftp_url'])

        if not files:
            print(f"  Could not list files. Try manual download from:")
            print(f"  {info['ftp_url']}")
            continue

        print(f"  Found {len(files)} total files")

        # Look for MaxQuant output files
        maxquant_files = [f for f in files if any(x in f.lower() for x in
                         ['proteingroups.txt', 'peptides.txt', 'evidence.txt',
                          'summary.txt', 'mqpar.xml'])]

        if maxquant_files:
            print(f"  Found {len(maxquant_files)} MaxQuant output files:")
            for mq_file in maxquant_files:
                print(f"    - {mq_file}")

            available_data[accession] = maxquant_files

            # Download proteinGroups.txt if available
            for filename in maxquant_files:
                if 'proteingroups.txt' in filename.lower():
                    local_file = dataset_dir / filename

                    if local_file.exists():
                        print(f"\n  {filename} already downloaded")
                        downloaded_files.append((accession, filename, local_file))
                    else:
                        print(f"\n  Downloading {filename}...")
                        file_url = info['ftp_url'] + filename
                        if download_file_http(file_url, local_file):
                            print(f"  Successfully downloaded to {local_file}")
                            downloaded_files.append((accession, filename, local_file))
                        else:
                            print(f"  Failed to download {filename}")
        else:
            print(f"  No MaxQuant output files found")

    # Summary
    print("\n" + "="*60)
    print("SUMMARY: MaxQuant data availability and downloads")
    print("="*60)

    if available_data:
        print("\nDatasets with MaxQuant files:")
        for accession, files in available_data.items():
            print(f"  {accession}: {len(files)} MaxQuant files")
    else:
        print("\nNo MaxQuant files found in any dataset")

    if downloaded_files:
        print("\nSuccessfully downloaded files:")
        for accession, filename, path in downloaded_files:
            size = path.stat().st_size / 1024 / 1024
            print(f"  {accession}/{filename} ({size:.1f} MB)")
            print(f"    -> {path}")
    else:
        print("\nNo files downloaded")

if __name__ == "__main__":
    main()
