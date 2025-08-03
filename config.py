"""
Configuration settings for DNA Codon Usage Analyzer
"""

# Analysis parameters
MIN_SEQUENCE_LENGTH = 300  # Minimum sequence length for analysis
MIN_CODONS_FOR_ANALYSIS = 100  # Minimum codons needed for reliable analysis

# Visualization settings
FIGURE_SIZE = (12, 8)
DPI = 300
COLOR_PALETTE = "husl"

# RSCU thresholds for bias detection
RSCU_HIGH_BIAS = 1.5  # Over-represented threshold
RSCU_LOW_BIAS = 0.5   # Under-represented threshold

# Output settings
OUTPUT_DIR = "results"
PLOTS_DIR = "plots"
REPORTS_DIR = "reports"

# File formats
SEQUENCE_FORMATS = ['fasta', 'fa', 'fna', 'ffn']
REPORT_FORMAT = 'csv'

# Genetic code table (standard)
GENETIC_CODE = {
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

# Amino acid properties
AMINO_ACID_PROPERTIES = {
    'A': {'name': 'Alanine', 'type': 'nonpolar'},
    'R': {'name': 'Arginine', 'type': 'positive'},
    'N': {'name': 'Asparagine', 'type': 'polar'},
    'D': {'name': 'Aspartate', 'type': 'negative'},
    'C': {'name': 'Cysteine', 'type': 'polar'},
    'Q': {'name': 'Glutamine', 'type': 'polar'},
    'E': {'name': 'Glutamate', 'type': 'negative'},
    'G': {'name': 'Glycine', 'type': 'nonpolar'},
    'H': {'name': 'Histidine', 'type': 'positive'},
    'I': {'name': 'Isoleucine', 'type': 'nonpolar'},
    'L': {'name': 'Leucine', 'type': 'nonpolar'},
    'K': {'name': 'Lysine', 'type': 'positive'},
    'M': {'name': 'Methionine', 'type': 'nonpolar'},
    'F': {'name': 'Phenylalanine', 'type': 'nonpolar'},
    'P': {'name': 'Proline', 'type': 'nonpolar'},
    'S': {'name': 'Serine', 'type': 'polar'},
    'T': {'name': 'Threonine', 'type': 'polar'},
    'W': {'name': 'Tryptophan', 'type': 'nonpolar'},
    'Y': {'name': 'Tyrosine', 'type': 'polar'},
    'V': {'name': 'Valine', 'type': 'nonpolar'}
}
