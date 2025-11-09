# PRIDE Data Download Instructions

## Overview

This directory contains scripts to download MaxQuant proteinGroups.txt files from PRIDE (PRoteomics IDEntifications) archive for Francisella TCS analysis.

## Datasets

The following 8 PRIDE datasets have been identified for analysis:

| Dataset ID | Description | Date |
|------------|-------------|------|
| PXD035145 | Rytter 2024 - macrophage infection | 2024/04 |
| PXD005747 | Link 2018 - phosphoproteomics | 2017/10 |
| PXD013074 | Klimentova 2019 - stress conditions | 2019/10 |
| PXD025439 | Valikangas 2022 - temperature series | 2023/01 |
| PXD016669 | Felix 2021 - ldcF mutant | 2020/01 |
| PXD001584 | Ramond 2015 - argP mutant | 2015/01 |
| PXD019739 | Kopeckova 2020 - GapA interactome | 2020/09 |
| PXD022406 | Klimentova 2021 - cross-species | 2021/03 |

## Download Methods

### Method 1: Python Script (Recommended)

The Python script `download_pride_data.py` attempts to automatically download files from PRIDE:

```bash
python3 scripts/download_pride_data.py
```

**Features:**
- Tests network connectivity to PRIDE
- Lists available MaxQuant files
- Downloads proteinGroups.txt files with progress tracking
- Handles errors gracefully
- Falls back to manual download instructions if network is restricted

**Requirements:**
- Python 3.6+
- requests library: `pip install requests`
- Unrestricted network access to ftp.pride.ebi.ac.uk

### Method 2: Bash Script (Alternative)

If the Python script cannot access PRIDE (e.g., in restricted environments), use the generated wget script:

```bash
chmod +x scripts/download_pride_manual.sh
./scripts/download_pride_manual.sh
```

This will recursively download MaxQuant output files using wget.

**Requirements:**
- wget installed
- Unrestricted network access

### Method 3: Manual Download

Visit each dataset URL directly in a web browser:

1. **PXD035145**: https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/04/PXD035145/
2. **PXD005747**: https://ftp.pride.ebi.ac.uk/pride/data/archive/2017/10/PXD005747/
3. **PXD013074**: https://ftp.pride.ebi.ac.uk/pride/data/archive/2019/10/PXD013074/
4. **PXD025439**: https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/01/PXD025439/
5. **PXD016669**: https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/01/PXD016669/
6. **PXD001584**: https://ftp.pride.ebi.ac.uk/pride/data/archive/2015/01/PXD001584/
7. **PXD019739**: https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/09/PXD019739/
8. **PXD022406**: https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/PXD022406/

Download the following files from each dataset:
- `*proteinGroups.txt` (primary file needed)
- `*peptides.txt` (optional)
- `*evidence.txt` (optional)
- `*summary.txt` (optional)
- `mqpar.xml` (MaxQuant parameters, optional)

Save files to: `data/pride/<DATASET_ID>/`

## File Organization

Downloaded files will be organized as:

```
data/pride/
├── PXD035145/
│   └── proteinGroups.txt
├── PXD005747/
│   └── proteinGroups.txt
├── PXD013074/
│   └── proteinGroups.txt
...
```

## Troubleshooting

### Network Access Issues

If you encounter network errors (403 Forbidden, DNS resolution failures):

1. This indicates restricted network access (common in sandboxed environments)
2. Run the scripts on your local machine with full network access
3. Use a VPN if accessing from restricted networks
4. Contact your network administrator if persistent

### API Access Blocked

The PRIDE Archive API (v2) may block automated requests. The scripts have been designed to:
- Use FTP/HTTP access instead of API
- Implement proper error handling
- Provide alternative download methods

### Large File Downloads

MaxQuant proteinGroups.txt files can be large (10-100+ MB):
- Ensure sufficient disk space (~1-2 GB total)
- Use stable network connection
- Downloads show progress indicators
- Partial downloads are automatically cleaned up on failure

## Next Steps

After downloading the data:

1. Verify files downloaded correctly:
   ```bash
   ls -lh data/pride/*/proteinGroups.txt
   ```

2. Proceed with TCS protein analysis using the DEqMS workflow
3. Refer to the main analysis scripts for data processing

## References

- PRIDE Archive: https://www.ebi.ac.uk/pride/
- ProteomeXchange: http://www.proteomexchange.org/
- MaxQuant: https://www.maxquant.org/

## Support

For issues with:
- **PRIDE data access**: Contact PRIDE support at pride-support@ebi.ac.uk
- **Script errors**: Check script logs and network connectivity
- **Analysis questions**: Refer to project documentation
