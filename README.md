# DNA Codon Usage Analyzer

A comprehensive bioinformatics tool for analyzing codon usage bias in genomic sequences

Understanding how organisms use their genetic code is fundamental to molecular biology. While the genetic code provides 64 possible codons to encode just 20 amino acids (plus stop signals), organisms don't use all codons equally. This project investigates **codon usage bias** the fascinating phenomenon where cells show distinct preferences for certain codons over their synonymous alternatives.

This analysis reveals critical biological insights including how genes are optimized for expression, evolutionary relationships between species, translational efficiency patterns, and even helps detect horizontal gene transfer events. By examining these patterns, we can better understand everything from bacterial pathogenesis to synthetic biology applications.

## Key Features

**Comprehensive RSCU Analysis**: The tool calculates Relative Synonymous Codon Usage values, the gold standard metric for quantifying codon bias. RSCU values reveal which codons are overrepresented (>1.5), underrepresented ( 1.5), blue bars show avoided codons (RSCU < 0.5), and gray bars represent codons used at expected frequencies.

### Key Biological Discoveries

Our analysis of three distinct genomic sequences revealed compelling evolutionary signatures:

**Comprehensive Scope**: We analyzed **468 individual codons** across three different organisms, providing robust statistical power for detecting meaningful bias patterns.

**Significant Bias Detection**: **37 codons** demonstrated statistically significant bias (RSCU values either greater than 1.5 or less than 0.5), indicating strong selective pressure on codon choice.

**Extreme Preferences**: **GTG coding for Valine** emerged as the most heavily preferred codon with an RSCU value of 3.643, suggesting this codon may be associated with highly expressed genes or optimal tRNA availability.

**Evolutionary Signatures**: The pronounced **preference for G/C-ending codons** strongly suggests these sequences originated from bacterial genomes, which typically show this pattern due to their codon-anticodon optimization strategies.

## Quick Start Guide

Getting started with the DNA Codon Usage Analyzer is straightforward. The tool is designed to work out of the box with sample data, while also supporting custom genomic sequences for specialized research applications.

### Installation and Setup

```bash
# Clone the repository to your local machine
git clone https://github.com/your-username/dna-codon-analyzer.git
cd dna-codon-analyzer

# Install all required Python dependencies
pip install -r requirements.txt

# Run the complete analysis pipeline
python codon_analyzer.py
```

The analysis will automatically download sample genomic sequences, perform comprehensive codon usage calculations, generate statistical reports, and create publication-ready visualizations. Results are saved in organized directories for easy access and sharing.

### What Happens During Analysis

The tool processes DNA sequences through several sophisticated steps: first cleaning and validating input sequences, then identifying optimal reading frames to avoid stop codons, extracting individual codons, calculating usage frequencies and RSCU values, detecting statistically significant bias patterns, and finally generating comprehensive reports and visualizations.

### Expected Output

Upon completion, you'll find detailed CSV reports containing all calculated metrics, a plots directory with high-resolution visualizations suitable for presentations, statistical summaries highlighting the most significant findings, and sample genomic data for testing and validation purposes.

This comprehensive analysis framework provides everything needed to understand codon usage patterns in your genomic data, whether for academic research, biotechnology applications, or synthetic biology projects.

