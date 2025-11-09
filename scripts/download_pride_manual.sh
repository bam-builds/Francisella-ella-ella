#!/bin/bash
# Manual download script for PRIDE datasets
# Run this script in an environment with network access

mkdir -p data/pride

# PXD035145: Rytter 2024 - macrophage infection
echo "Downloading PXD035145..."
mkdir -p data/pride/PXD035145
cd data/pride/PXD035145
wget -r -np -nH --cut-dirs=5 -R "index.html*" -e robots=off \
  "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/04/PXD035145/" \
  -A "*proteinGroups.txt,*peptides.txt,*evidence.txt,*summary.txt,mqpar.xml"
cd ../../..

# PXD005747: Link 2018 - phosphoproteomics
echo "Downloading PXD005747..."
mkdir -p data/pride/PXD005747
cd data/pride/PXD005747
wget -r -np -nH --cut-dirs=5 -R "index.html*" -e robots=off \
  "https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/10/PXD005747/" \
  -A "*proteinGroups.txt,*peptides.txt,*evidence.txt,*summary.txt,mqpar.xml"
cd ../../..

# PXD013074: Klimentova 2019 - stress conditions
echo "Downloading PXD013074..."
mkdir -p data/pride/PXD013074
cd data/pride/PXD013074
wget -r -np -nH --cut-dirs=5 -R "index.html*" -e robots=off \
  "https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/10/PXD013074/" \
  -A "*proteinGroups.txt,*peptides.txt,*evidence.txt,*summary.txt,mqpar.xml"
cd ../../..

# PXD025439: Valikangas 2022 - temperature series
echo "Downloading PXD025439..."
mkdir -p data/pride/PXD025439
cd data/pride/PXD025439
wget -r -np -nH --cut-dirs=5 -R "index.html*" -e robots=off \
  "https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/01/PXD025439/" \
  -A "*proteinGroups.txt,*peptides.txt,*evidence.txt,*summary.txt,mqpar.xml"
cd ../../..

# PXD016669: Felix 2021 - ldcF mutant
echo "Downloading PXD016669..."
mkdir -p data/pride/PXD016669
cd data/pride/PXD016669
wget -r -np -nH --cut-dirs=5 -R "index.html*" -e robots=off \
  "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD016669/" \
  -A "*proteinGroups.txt,*peptides.txt,*evidence.txt,*summary.txt,mqpar.xml"
cd ../../..

# PXD001584: Ramond 2015 - argP mutant
echo "Downloading PXD001584..."
mkdir -p data/pride/PXD001584
cd data/pride/PXD001584
wget -r -np -nH --cut-dirs=5 -R "index.html*" -e robots=off \
  "https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/01/PXD001584/" \
  -A "*proteinGroups.txt,*peptides.txt,*evidence.txt,*summary.txt,mqpar.xml"
cd ../../..

# PXD019739: Kopeckova 2020 - GapA interactome
echo "Downloading PXD019739..."
mkdir -p data/pride/PXD019739
cd data/pride/PXD019739
wget -r -np -nH --cut-dirs=5 -R "index.html*" -e robots=off \
  "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD019739/" \
  -A "*proteinGroups.txt,*peptides.txt,*evidence.txt,*summary.txt,mqpar.xml"
cd ../../..

# PXD022406: Klimentova 2021 - cross-species
echo "Downloading PXD022406..."
mkdir -p data/pride/PXD022406
cd data/pride/PXD022406
wget -r -np -nH --cut-dirs=5 -R "index.html*" -e robots=off \
  "https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/PXD022406/" \
  -A "*proteinGroups.txt,*peptides.txt,*evidence.txt,*summary.txt,mqpar.xml"
cd ../../..

