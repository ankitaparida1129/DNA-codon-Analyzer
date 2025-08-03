#!/usr/bin/env python3
"""
Generate test genomic sequences for codon usage analysis
"""

import random
import argparse

def generate_realistic_sequence(length=3000, organism_type="balanced"):
    """Generate realistic DNA sequences with organism-specific codon bias"""
    
    # Define codon preferences for different organism types
    codon_preferences = {
        "balanced": {  # Equal usage
            'TTT': 1, 'TTC': 1, 'TTA': 1, 'TTG': 1,
            'TCT': 1, 'TCC': 1, 'TCA': 1, 'TCG': 1,
            'TAT': 1, 'TAC': 1, 'TGT': 1, 'TGC': 1, 'TGG': 1,
            'CTT': 1, 'CTC': 1, 'CTA': 1, 'CTG': 1,
            'CCT': 1, 'CCC': 1, 'CCA': 1, 'CCG': 1,
            'CAT': 1, 'CAC': 1, 'CAA': 1, 'CAG': 1,
            'CGT': 1, 'CGC': 1, 'CGA': 1, 'CGG': 1,
            'ATT': 1, 'ATC': 1, 'ATA': 1, 'ATG': 1,
            'ACT': 1, 'ACC': 1, 'ACA': 1, 'ACG': 1,
            'AAT': 1, 'AAC': 1, 'AAA': 1, 'AAG': 1,
            'AGT': 1, 'AGC': 1, 'AGA': 1, 'AGG': 1,
            'GTT': 1, 'GTC': 1, 'GTA': 1, 'GTG': 1,
            'GCT': 1, 'GCC': 1, 'GCA': 1, 'GCG': 1,
            'GAT': 1, 'GAC': 1, 'GAA': 1, 'GAG': 1,
            'GGT': 1, 'GGC': 1, 'GGA': 1, 'GGG': 1
        },
        
        "gc_rich": {  # Prefers G/C endings (like some bacteria)
            'TTT': 0.3, 'TTC': 1.5, 'TTA': 0.2, 'TTG': 1.2,
            'TCT': 0.5, 'TCC': 1.8, 'TCA': 0.4, 'TCG': 1.6,
            'TAT': 0.4, 'TAC': 1.6, 'TGT': 0.3, 'TGC': 1.7, 'TGG': 1,
            'CTT': 0.3, 'CTC': 1.7, 'CTA': 0.2, 'CTG': 2.0,
            'CCT': 0.4, 'CCC': 1.8, 'CCA': 0.3, 'CCG': 1.9,
            'CAT': 0.4, 'CAC': 1.6, 'CAA': 0.5, 'CAG': 1.8,
            'CGT': 0.8, 'CGC': 2.0, 'CGA': 0.3, 'CGG': 1.9,
            'ATT': 0.4, 'ATC': 1.6, 'ATA': 0.2, 'ATG': 1,
            'ACT': 0.5, 'ACC': 1.8, 'ACA': 0.4, 'ACG': 1.7,
            'AAT': 0.4, 'AAC': 1.6, 'AAA': 0.6, 'AAG': 1.4,
            'AGT': 0.4, 'AGC': 1.6, 'AGA': 0.3, 'AGG': 0.5,
            'GTT': 0.4, 'GTC': 1.6, 'GTA': 0.3, 'GTG': 1.7,
            'GCT': 0.5, 'GCC': 1.8, 'GCA': 0.4, 'GCG': 1.7,
            'GAT': 0.4, 'GAC': 1.6, 'GAA': 0.6, 'GAG': 1.4,
            'GGT': 0.6, 'GGC': 1.8, 'GGA': 0.5, 'GGG': 1.1
        },
        
        "at_rich": {  # Prefers A/T endings (like some eukaryotes)
            'TTT': 1.5, 'TTC': 0.5, 'TTA': 1.8, 'TTG': 0.8,
            'TCT': 1.6, 'TCC': 0.4, 'TCA': 1.7, 'TCG': 0.3,
            'TAT': 1.6, 'TAC': 0.4, 'TGT': 1.7, 'TGC': 0.3, 'TGG': 1,
            'CTT': 1.7, 'CTC': 0.3, 'CTA': 1.8, 'CTG': 0.2,
            'CCT': 1.6, 'CCC': 0.2, 'CCA': 1.7, 'CCG': 0.1,
            'CAT': 1.6, 'CAC': 0.4, 'CAA': 1.5, 'CAG': 0.2,
            'CGT': 1.2, 'CGC': 0.0, 'CGA': 1.7, 'CGG': 0.1,
            'ATT': 1.6, 'ATC': 0.4, 'ATA': 1.8, 'ATG': 1,
            'ACT': 1.5, 'ACC': 0.2, 'ACA': 1.6, 'ACG': 0.3,
            'AAT': 1.6, 'AAC': 0.4, 'AAA': 1.4, 'AAG': 0.6,
            'AGT': 1.6, 'AGC': 0.4, 'AGA': 1.7, 'AGG': 1.5,
            'GTT': 1.6, 'GTC': 0.4, 'GTA': 1.7, 'GTG': 0.3,
            'GCT': 1.5, 'GCC': 0.2, 'GCA': 1.6, 'GCG': 0.3,
            'GAT': 1.6, 'GAC': 0.4, 'GAA': 1.4, 'GAG': 0.6,
            'GGT': 1.4, 'GGC': 0.2, 'GGA': 1.5, 'GGG': 0.9
        }
    }
    
    preferences = codon_preferences.get(organism_type, codon_preferences["balanced"])
    
    # Create weighted codon list
    weighted_codons = []
    for codon, weight in preferences.items():
        weighted_codons.extend([codon] * int(weight * 10))
    
    # Generate sequence
    sequence = "ATG"  # Start with start codon
    
    # Add codons to reach desired length
    while len(sequence) < length - 3:
        codon = random.choice(weighted_codons)
        sequence += codon
    
    # End with stop codon
    sequence += "TAA"
    
    return sequence

def main():
    parser = argparse.ArgumentParser(description='Generate test genomic sequences')
    parser.add_argument('--organisms', nargs='+', 
                       choices=['balanced', 'gc_rich', 'at_rich'],
                       default=['balanced', 'gc_rich', 'at_rich'],
                       help='Organism types to generate')
    parser.add_argument('--length', type=int, default=3000,
                       help='Sequence length (default: 3000)')
    parser.add_argument('--output', default='test_sequences.fasta',
                       help='Output filename (default: test_sequences.fasta)')
    
    args = parser.parse_args()
    
    print(f"Generating test sequences...")
    
    with open(args.output, 'w') as f:
        for i, org_type in enumerate(args.organisms):
            for replicate in range(3):  # Generate 3 replicates per organism
                seq_id = f"{org_type}_replicate_{replicate + 1}"
                sequence = generate_realistic_sequence(args.length, org_type)
                
                f.write(f">{seq_id} | Generated {org_type} sequence, length={len(sequence)}\n")
                
                # Write sequence in 80-character lines
                for j in range(0, len(sequence), 80):
                    f.write(sequence[j:j+80] + "\n")
    
    print(f"Generated sequences saved to '{args.output}'")
    print(f"Organisms: {', '.join(args.organisms)}")
    print(f"Sequence length: {args.length}")

if __name__ == "__main__":
    main()
