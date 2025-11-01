#!/usr/bin/env python3
"""
Francisella Two-Component Systems Analysis
Simplified runnable version for demonstration
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style('whitegrid')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12

# Import custom modules
import sys
sys.path.insert(0, str(Path(__file__).parent / 'src' / 'python'))
from data_retrieval import ProteomeXchangeDownloader
from secretion_prediction import SecretionPredictor
from ptm_analysis import PhosphoproteomicsAnalyzer

print('=' * 60)
print('Francisella TCS Analysis Pipeline')
print('=' * 60)
print('Libraries loaded successfully!\n')


def generate_sample_data():
    """Generate sample data for demonstration"""
    print("Generating sample data...")

    # Create sample protein groups data
    np.random.seed(42)
    n_proteins = 200

    protein_data = {
        'Protein IDs': [f'P{i:05d}' for i in range(n_proteins)],
        'Gene names': [
            'QseB', 'PmrA', 'BfpR', 'KdpE', 'QseC', 'KdpD', 'FTN_1453',
            'IglA', 'IglB', 'IglC', 'IglD', 'PdpA', 'PdpB',
            'ChiA', 'ChiB', 'CbpA', 'PepO'
        ] + [f'Gene{i}' for i in range(n_proteins - 17)],
        'Peptides': np.random.randint(2, 20, n_proteins),
        'Intensity': np.random.lognormal(20, 2, n_proteins)
    }
    protein_groups = pd.DataFrame(protein_data)

    # Create sample phosphorylation data
    n_phospho = 50
    phospho_data = {
        'Proteins': np.random.choice(protein_data['Gene names'][:20], n_phospho),
        'Position': np.random.randint(10, 300, n_phospho),
        'Amino acid': np.random.choice(['S', 'T', 'Y'], n_phospho, p=[0.6, 0.3, 0.1]),
        'Localization prob': np.random.beta(8, 2, n_phospho),
        'Intensity': np.random.lognormal(18, 2, n_phospho)
    }
    phospho_sites = pd.DataFrame(phospho_data)

    # Create sample differential expression data
    de_data = {
        'Protein IDs': protein_data['Protein IDs'],
        'Gene names': protein_data['Gene names'],
        'log2FC': np.random.normal(0, 2, n_proteins),
        'adj.P.Val': np.random.beta(0.5, 5, n_proteins)
    }
    de_results = pd.DataFrame(de_data)

    # Create sample protein sequences
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    sequences = []
    for i, gene in enumerate(protein_data['Gene names'][:20]):
        # Create a realistic-ish protein sequence
        length = np.random.randint(100, 500)
        # Add signal peptide to some proteins
        if np.random.random() > 0.7:
            seq = 'MKR' + ''.join(np.random.choice(list('AVILMFW'), 15)) + \
                  ''.join(np.random.choice(list(amino_acids), length - 18))
        else:
            seq = ''.join(np.random.choice(list(amino_acids), length))
        sequences.append({'protein_id': protein_data['Protein IDs'][i],
                         'gene_name': gene, 'sequence': seq})

    sequences_df = pd.DataFrame(sequences)

    # Save data files
    output_dir = Path('data/processed')
    output_dir.mkdir(parents=True, exist_ok=True)

    protein_groups.to_csv(output_dir / 'PXD009225_proteinGroups.txt', sep='\t', index=False)
    phospho_sites.to_csv(output_dir / 'PXD009225_Phospho_STY_Sites.txt', sep='\t', index=False)
    de_results.to_csv('results/tables/differential_expression.csv', index=False)
    sequences_df.to_csv('data/processed/protein_sequences.csv', index=False)

    print(f"  Created {len(protein_groups)} protein entries")
    print(f"  Created {len(phospho_sites)} phosphorylation sites")
    print(f"  Created differential expression data")
    print(f"  Created {len(sequences_df)} protein sequences\n")

    return protein_groups, phospho_sites, de_results, sequences_df


def main():
    """Main analysis pipeline"""

    # Step 1: Generate sample data
    print("Step 1: Data Generation")
    print("-" * 60)
    protein_groups, phospho_sites, de_results, sequences_df = generate_sample_data()

    # Step 2: Load and filter data
    print("\nStep 2: Data Filtering")
    print("-" * 60)
    proteins_filtered = protein_groups[
        (~protein_groups['Protein IDs'].str.contains('REV__|CON__', na=False)) &
        (protein_groups['Peptides'] >= 2)
    ]
    print(f"Filtered to {len(proteins_filtered)} high-confidence proteins")

    phospho_filtered = phospho_sites[phospho_sites['Localization prob'] > 0.75]
    print(f"Filtered to {len(phospho_filtered)} high-confidence phosphosites\n")

    # Step 3: Identify TCS proteins
    print("Step 3: TCS Protein Identification")
    print("-" * 60)
    tcs_proteins = {
        'Response Regulators': ['QseB', 'PmrA', 'BfpR', 'KdpE'],
        'Sensor Kinases': ['QseC', 'KdpD', 'FTN_1453'],
        'FPI Genes': ['IglA', 'IglB', 'IglC', 'IglD', 'PdpA', 'PdpB'],
        'Secreted': ['ChiA', 'ChiB', 'CbpA', 'PepO']
    }

    tcs_found = {}
    for category, proteins in tcs_proteins.items():
        found = []
        for protein in proteins:
            matches = proteins_filtered[
                proteins_filtered['Gene names'].str.contains(protein, na=False, case=False)
            ]
            if not matches.empty:
                found.append(protein)
        tcs_found[category] = found
        print(f"{category}: {len(found)} found - {found}")

    # Step 4: Phosphoproteomics analysis
    print("\n\nStep 4: Phosphoproteomics Analysis")
    print("-" * 60)
    analyzer = PhosphoproteomicsAnalyzer()

    tcs_phospho = phospho_filtered[
        phospho_filtered['Proteins'].str.contains(
            '|'.join(['QseB', 'PmrA', 'QseC', 'KdpD', 'BfpR']),
            na=False, case=False
        )
    ]

    if not tcs_phospho.empty:
        print(f"Found {len(tcs_phospho)} phosphorylation sites on TCS proteins:")
        for idx, row in tcs_phospho.iterrows():
            print(f"  {row['Proteins']}: {row['Amino acid']}{row['Position']} "
                  f"(Prob: {row['Localization prob']:.2f})")
    else:
        print("No TCS phosphosites found in this dataset")

    # Create phosphosite distribution plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))

    aa_counts = phospho_filtered['Amino acid'].value_counts()
    ax1.bar(aa_counts.index, aa_counts.values, color='steelblue')
    ax1.set_title('Phosphorylated Amino Acids')
    ax1.set_xlabel('Amino Acid')
    ax1.set_ylabel('Number of Sites')

    ax2.hist(phospho_filtered['Localization prob'], bins=30, edgecolor='black', color='coral')
    ax2.axvline(x=0.75, color='red', linestyle='--', linewidth=2, label='Quality threshold')
    ax2.set_title('Phosphosite Localization Probability')
    ax2.set_xlabel('Localization Probability')
    ax2.set_ylabel('Number of Sites')
    ax2.legend()

    plt.tight_layout()
    plt.savefig('results/figures/phosphosite_analysis.png', dpi=300, bbox_inches='tight')
    print("\nSaved figure: results/figures/phosphosite_analysis.png")

    # Step 5: Secretion signal prediction
    print("\n\nStep 5: Secretion Signal Prediction")
    print("-" * 60)
    predictor = SecretionPredictor()

    predictions = []
    for idx, row in sequences_df.iterrows():
        result = predictor.predict(row['sequence'])
        predictions.append({
            'protein_id': row['protein_id'],
            'gene_name': row['gene_name'],
            'has_signal': result['has_signal'],
            'probability': result['probability']
        })

    prediction_df = pd.DataFrame(predictions)
    secreted = prediction_df[prediction_df['probability'] > 0.7]

    print(f"Analyzed {len(sequences_df)} protein sequences")
    print(f"Found {len(secreted)} proteins with high-confidence secretion signals")
    print(f"Secretion rate: {len(secreted)/len(predictions)*100:.1f}%")
    print(f"\nTop secreted proteins:")
    for idx, row in secreted.head(5).iterrows():
        print(f"  {row['gene_name']}: {row['probability']:.2f}")

    # Step 6: Differential expression analysis
    print("\n\nStep 6: Differential Expression Analysis")
    print("-" * 60)

    log2fc_threshold = 1.5
    pvalue_threshold = 0.05

    sig_up = de_results[
        (de_results['log2FC'] > log2fc_threshold) &
        (de_results['adj.P.Val'] < pvalue_threshold)
    ]

    sig_down = de_results[
        (de_results['log2FC'] < -log2fc_threshold) &
        (de_results['adj.P.Val'] < pvalue_threshold)
    ]

    print(f"Significantly upregulated: {len(sig_up)} proteins")
    print(f"Significantly downregulated: {len(sig_down)} proteins")

    # Volcano plot
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot all proteins
    ax.scatter(de_results['log2FC'], -np.log10(de_results['adj.P.Val']),
               alpha=0.5, s=20, color='gray', label='Not significant')

    # Highlight significant proteins
    if len(sig_up) > 0:
        ax.scatter(sig_up['log2FC'], -np.log10(sig_up['adj.P.Val']),
                   alpha=0.7, s=40, color='red', label=f'Up ({len(sig_up)})')
    if len(sig_down) > 0:
        ax.scatter(sig_down['log2FC'], -np.log10(sig_down['adj.P.Val']),
                   alpha=0.7, s=40, color='blue', label=f'Down ({len(sig_down)})')

    # Add threshold lines
    ax.axhline(y=-np.log10(pvalue_threshold), color='black', linestyle='--', alpha=0.5)
    ax.axvline(x=log2fc_threshold, color='black', linestyle='--', alpha=0.5)
    ax.axvline(x=-log2fc_threshold, color='black', linestyle='--', alpha=0.5)

    ax.set_xlabel('log2 Fold Change', fontsize=12)
    ax.set_ylabel('-log10(adjusted P-value)', fontsize=12)
    ax.set_title('Volcano Plot: Differential Expression Analysis', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('results/figures/volcano_plot.png', dpi=300, bbox_inches='tight')
    print("\nSaved figure: results/figures/volcano_plot.png")

    # Step 7: Integration
    print("\n\nStep 7: TCS-Regulated Secreted Proteins")
    print("-" * 60)

    integrated = prediction_df.merge(
        de_results,
        left_on='protein_id',
        right_on='Protein IDs',
        how='inner'
    )

    secreted_de = integrated[
        (integrated['probability'] > 0.7) &
        (abs(integrated['log2FC']) > log2fc_threshold) &
        (integrated['adj.P.Val'] < pvalue_threshold)
    ]

    print(f"Found {len(secreted_de)} secreted proteins with differential expression")

    if len(secreted_de) > 0:
        print("\nTop candidates:")
        for idx, row in secreted_de.head(5).iterrows():
            direction = 'UP' if row['log2FC'] > 0 else 'DOWN'
            print(f"  {row['gene_name']}: {direction} (log2FC={row['log2FC']:.2f}, "
                  f"p={row['adj.P.Val']:.3e}, secretion_prob={row['probability']:.2f})")

    # Step 8: Export results
    print("\n\nStep 8: Exporting Results")
    print("-" * 60)

    output_dir = Path('results/tables')
    output_dir.mkdir(parents=True, exist_ok=True)

    # Export to CSV
    if len(secreted_de) > 0:
        secreted_de.to_csv(output_dir / 'tcs_regulated_secreted_proteins.csv', index=False)
        print("Saved: results/tables/tcs_regulated_secreted_proteins.csv")

    if len(tcs_phospho) > 0:
        tcs_phospho.to_csv(output_dir / 'tcs_phosphorylation_sites.csv', index=False)
        print("Saved: results/tables/tcs_phosphorylation_sites.csv")

    # Summary statistics
    print("\n\nSummary Statistics")
    print("=" * 60)
    summary = {
        'Total proteins analyzed': len(proteins_filtered),
        'Phosphorylation sites': len(phospho_filtered),
        'TCS phosphosites': len(tcs_phospho),
        'Predicted secreted proteins': len(secreted),
        'Upregulated proteins': len(sig_up),
        'Downregulated proteins': len(sig_down),
        'TCS-regulated secreted proteins': len(secreted_de)
    }

    for metric, count in summary.items():
        print(f"{metric:.<45} {count}")

    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)


if __name__ == '__main__':
    main()

