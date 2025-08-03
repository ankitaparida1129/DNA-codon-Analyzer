#!/usr/bin/env python3
"""
DNA Codon Usage Analyzer
A bioinformatics tool to analyze codon usage bias in DNA sequences
Author: Your Name
Date: August 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
import requests
from Bio import SeqIO
from Bio.Seq import Seq
#from Bio.SeqUtils import CodonUsage
import io
import warnings
warnings.filterwarnings('ignore')

class CodonAnalyzer:
    def __init__(self):
        """Initialize the CodonAnalyzer with genetic code dictionary"""
        self.genetic_code = {
            'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
            'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
            'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
            'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
            'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
            'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
            'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
            'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
            'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
            'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
            'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
            'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
            'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
        }
        
        self.codon_counts = defaultdict(int)
        self.amino_acid_counts = defaultdict(int)
        self.sequences_analyzed = 0
        
    def clean_sequence(self, sequence):
        """Clean DNA sequence by removing non-nucleotide characters"""
        return ''.join(c.upper() for c in sequence if c.upper() in 'ATCG')
    
    def extract_codons(self, sequence):
        """Extract codons from DNA sequence"""
        sequence = self.clean_sequence(sequence)
        codons = []
        
        # Find the best reading frame (longest without stop codons)
        best_frame = 0
        max_length = 0
        
        for frame in range(3):
            temp_codons = [sequence[i:i+3] for i in range(frame, len(sequence)-2, 3) 
                          if len(sequence[i:i+3]) == 3]
            
            # Count codons before first stop codon
            length = 0
            for codon in temp_codons:
                if codon in self.genetic_code and self.genetic_code[codon] == '*':
                    break
                length += 1
            
            if length > max_length:
                max_length = length
                best_frame = frame
        
        # Extract codons from best frame
        for i in range(best_frame, len(sequence)-2, 3):
            codon = sequence[i:i+3]
            if len(codon) == 3 and codon in self.genetic_code:
                if self.genetic_code[codon] != '*':  # Skip stop codons
                    codons.append(codon)
        
        return codons
    
    def analyze_sequence(self, sequence, sequence_id="Unknown"):
        """Analyze a single DNA sequence for codon usage"""
        codons = self.extract_codons(sequence)
        
        if not codons:
            print(f"Warning: No valid codons found in sequence {sequence_id}")
            return
        
        for codon in codons:
            self.codon_counts[codon] += 1
            amino_acid = self.genetic_code[codon]
            self.amino_acid_counts[amino_acid] += 1
        
        self.sequences_analyzed += 1
        print(f"Analyzed sequence {sequence_id}: {len(codons)} codons found")
    
    def calculate_codon_usage_bias(self):
        """Calculate Relative Synonymous Codon Usage (RSCU)"""
        # Group codons by amino acid
        aa_codons = defaultdict(list)
        for codon, aa in self.genetic_code.items():
            if aa != '*':  # Skip stop codons
                aa_codons[aa].append(codon)
        
        rscu_values = {}
        
        for aa, codons in aa_codons.items():
            if self.amino_acid_counts[aa] == 0:
                continue
                
            # Calculate RSCU for each codon of this amino acid
            for codon in codons:
                observed = self.codon_counts[codon]
                expected = self.amino_acid_counts[aa] / len(codons)
                
                if expected > 0:
                    rscu = observed / expected
                else:
                    rscu = 0
                
                rscu_values[codon] = rscu
        
        return rscu_values
    
    def generate_report(self):
        """Generate comprehensive codon usage report"""
        total_codons = sum(self.codon_counts.values())
        rscu_values = self.calculate_codon_usage_bias()
        
        print(f"\n{'='*60}")
        print(f"CODON USAGE ANALYSIS REPORT")
        print(f"{'='*60}")
        print(f"Sequences analyzed: {self.sequences_analyzed}")
        print(f"Total codons: {total_codons}")
        print(f"Unique codons found: {len(self.codon_counts)}")
        
        # Create detailed DataFrame
        report_data = []
        for codon in sorted(self.codon_counts.keys()):
            aa = self.genetic_code[codon]
            count = self.codon_counts[codon]
            frequency = count / total_codons * 100
            rscu = rscu_values.get(codon, 0)
            
            report_data.append({
                'Codon': codon,
                'Amino_Acid': aa,
                'Count': count,
                'Frequency_%': round(frequency, 2),
                'RSCU': round(rscu, 3)
            })
        
        df = pd.DataFrame(report_data)
        return df
    
    def create_visualizations(self, output_dir="plots"):
        """Create comprehensive visualizations"""
        import os
        os.makedirs(output_dir, exist_ok=True)
        
        rscu_values = self.calculate_codon_usage_bias()
        total_codons = sum(self.codon_counts.values())
        
        # Set style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # 1. Codon Frequency Heatmap
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Create matrix for heatmap
        bases = ['T', 'C', 'A', 'G']
        codon_matrix = np.zeros((4, 16))
        codon_labels = []
        
        for i, first in enumerate(bases):
            row_labels = []
            for j, (second, third) in enumerate([(s, t) for s in bases for t in bases]):
                codon = first + second + third
                row_labels.append(codon)
                if codon in self.codon_counts:
                    codon_matrix[i, j] = self.codon_counts[codon] / total_codons * 100
        
        codon_labels = [f"{bases[i//4]}{bases[(i%4)//1]}{bases[i%4]}" for i in range(16)]
        
        im = ax.imshow(codon_matrix, cmap='YlOrRd', aspect='auto')
        ax.set_xticks(range(16))
        ax.set_xticklabels([f"{s}{t}" for s in bases for t in bases], rotation=45)
        ax.set_yticks(range(4))
        ax.set_yticklabels(bases)
        ax.set_xlabel('Second & Third Position')
        ax.set_ylabel('First Position')
        ax.set_title('Codon Usage Frequency Heatmap (%)')
        
        # Add colorbar
        cbar = plt.colorbar(im)
        cbar.set_label('Frequency (%)')
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/codon_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. RSCU Distribution
        fig, ax = plt.subplots(figsize=(15, 8))
        
        codons = list(rscu_values.keys())
        rscu_vals = list(rscu_values.values())
        
        colors = ['red' if x > 1.5 else 'blue' if x < 0.5 else 'gray' for x in rscu_vals]
        
        bars = ax.bar(range(len(codons)), rscu_vals, color=colors, alpha=0.7)
        ax.axhline(y=1, color='black', linestyle='--', alpha=0.5, label='Expected (RSCU=1)')
        ax.set_xlabel('Codons')
        ax.set_ylabel('RSCU Value')
        ax.set_title('Relative Synonymous Codon Usage (RSCU)')
        ax.set_xticks(range(len(codons)))
        ax.set_xticklabels(codons, rotation=45, ha='right')
        ax.legend()
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/rscu_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 3. Amino Acid Usage
        fig, ax = plt.subplots(figsize=(12, 6))
        
        aa_names = {
            'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
            'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
            'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
            'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
        }
        
        aa_counts = [(aa_names.get(aa, aa), count) for aa, count in self.amino_acid_counts.items()]
        aa_counts.sort(key=lambda x: x[1], reverse=True)
        
        aas, counts = zip(*aa_counts) if aa_counts else ([], [])
        
        ax.bar(aas, counts, color='skyblue', alpha=0.8)
        ax.set_xlabel('Amino Acids')
        ax.set_ylabel('Count')
        ax.set_title('Amino Acid Usage Distribution')
        plt.xticks(rotation=45)
        
        plt.tight_layout()
        plt.savefig(f'{output_dir}/amino_acid_usage.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"\nVisualizations saved in '{output_dir}/' directory:")
        print("- codon_heatmap.png")
        print("- rscu_distribution.png") 
        print("- amino_acid_usage.png")

def download_sample_data():
    """Download sample genomic data for analysis"""
    print("Downloading sample genomic sequences...")
    
    # Sample sequences from different organisms (simplified for demo)
    sample_sequences = {
        "E_coli_sample": "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGTTCAGTGGTTCGTAGGGCTTTGCCCCGCCTGTTTTGACCGCTGTTATGTCTGGTTCCGCGTCCCGTGATGAAACAGGAATCGTCATCCGCGGGACGCAGCAGGGCAAG",
        
        "Human_sample": "ATGGCGGCGCTGAGCGGTGGCGGCGGCGGCGCTGAGCGGTGGCGGCGGCGGCGCTGAGCGGTGGCGGCGGCGGCGCTGAGCGGTGGCGGCGGCGGCGCTGAGCGGTGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGCGGCGGCGGCGGCGCTGAGCGGA",
        
        "Yeast_sample": "ATGTCTTCTGTTCCAAGAAACTCAACACTGGTGAAATTGGTGGTGGTGGTGGTGGTGGTGCTGGTGGTGGTGGTGGTGCTGGTGCTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGTGCTGGTGGTGGTGGT"
    }
    
    # Save sequences to files
    with open('sample_sequences.fasta', 'w') as f:
        for seq_id, sequence in sample_sequences.items():
            f.write(f">{seq_id}\n{sequence}\n")
    
    print("Sample sequences saved to 'sample_sequences.fasta'")
    return sample_sequences

def main():
    """Main analysis pipeline"""
    print("DNA Codon Usage Analyzer")
    print("=" * 40)
    
    # Initialize analyzer
    analyzer = CodonAnalyzer()
    
    # Download sample data
    sequences = download_sample_data()
    
    # Analyze sequences
    print("\nAnalyzing sequences...")
    for seq_id, sequence in sequences.items():
        analyzer.analyze_sequence(sequence, seq_id)
    
    # Generate report
    report_df = analyzer.generate_report()
    
    # Save detailed report
    report_df.to_csv('codon_usage_report.csv', index=False)
    print(f"\nDetailed report saved to 'codon_usage_report.csv'")
    
    # Create visualizations
    analyzer.create_visualizations()
    
    # Display summary statistics
    print(f"\n{'='*60}")
    print("SUMMARY STATISTICS")
    print(f"{'='*60}")
    
    # Most and least used codons
    most_used = report_df.nlargest(5, 'Count')[['Codon', 'Amino_Acid', 'Count', 'Frequency_%']]
    least_used = report_df.nsmallest(5, 'Count')[['Codon', 'Amino_Acid', 'Count', 'Frequency_%']]
    
    print("\nTop 5 Most Used Codons:")
    print(most_used.to_string(index=False))
    
    print("\nTop 5 Least Used Codons:")
    print(least_used.to_string(index=False))
    
    # Codon bias analysis
    biased_codons = report_df[(report_df['RSCU'] > 1.5) | (report_df['RSCU'] < 0.5)]
    print(f"\nCodons with significant bias (RSCU > 1.5 or < 0.5): {len(biased_codons)}")
    
    if len(biased_codons) > 0:
        print(biased_codons[['Codon', 'Amino_Acid', 'RSCU']].to_string(index=False))
    
    print(f"\n{'='*60}")
    print("Analysis completed successfully!")
    print("Check the generated files:")
    print("- codon_usage_report.csv (detailed data)")
    print("- plots/ directory (visualizations)")
    print(f"{'='*60}")

if __name__ == "__main__":
    main()
