# RVpt (Resistance Genes and Virulence Factors Prediction of transfer)

RVpt is a bioinformatics tool designed to predict and analyze the transfer potential of resistance genes and virulence factors in bacterial genomes. By integrating multiple analysis methods, it provides comprehensive assessment of gene mobility and associated mobile genetic elements.

## Key Features

### Gene Detection and Analysis

- Antimicrobial Resistance Genes (ARGs) detection and annotation
- Virulence Factors (VFs) detection and annotation
- Mobile Genetic Elements (MGEs) identification, including:
  - Insertion Sequences (IS)
  - Transposons (Tn)
  - Integrons (IN)

### Transfer Potential Analysis

- Plasmid association analysis
- Prophage region prediction
- Genetic context analysis
  - Adjacent mobile element detection
  - Transfer potential assessment
  - Transfer evidence compilation

## Dependencies

### Python Requirements

- Python 3.8+
- pandas
- biopython

### External Tools

- BLAST+ (for sequence alignment)
- Prokka (for genome annotation)
- Platon (for plasmid analysis)
- PhiSpy (for prophage prediction)

### Required Databases

The program requires the following database files:

- IS/Tn/IN database: `RVpt/ISTnIN/istnin.fas`
- ResFinder database: `RVpt/RESfinder/resfinder.fas`
- VFDB database: `RVpt/VFDB/vfdb.fas`
- Platon database: default location `RVpt/platon/`

## Installation

1. Clone the repository:

```bash
git clone https://github.com/leejianhai/RVpt.git
cd RVpt
```

2. Create and activate conda environment:

```bash
conda create -n rvpt python=3.8
conda activate rvpt
```

3. Install dependencies:

```bash
# Install Python packages
conda install pandas biopython

# Install analysis tools
conda install -c bioconda prokka
conda install -c bioconda blast
conda install -c bioconda platon
conda install -c bioconda phispy
```

## Usage

### Basic Usage

```bash
python RVpt.py input_genome.fasta [--outdir output_directory] [--db platon_db_path] [--threads num_threads]
```

### Parameters

- `input_genome.fasta`: Input genome sequence file (required)
- `--outdir`: Output directory path (optional, default: ../result)
- `--db`: Platon database path (optional, default: ../RVpt/platon)
- `--threads`, `-t`: Number of threads to use (optional, default: 4)

## Output Files

The program generates the following files in the specified output directory:

### Gene Annotation Results

- `prokka_results/`: Genome annotation results directory

### Sequence Analysis Results

- `ISTnIN_filtered_with_locations.tsv`: Mobile genetic elements analysis results
- `RESfinder_filtered_with_locations.tsv`: Resistance genes analysis results
- `VFDB_filtered_with_locations.tsv`: Virulence factors analysis results

### Mobile Element Analysis

- `platon_results/`: Plasmid analysis results directory
- `platon_results_summary.tsv`: Plasmid analysis summary
- `phispy_results/`: Prophage prediction results directory

### Transfer Potential Analysis

- `resistance_genes_context.tsv`: Resistance genes transfer potential analysis
  - Includes gene location, context features, and transfer evidence
- `virulence_genes_context.tsv`: Virulence factors transfer potential analysis
  - Includes gene location, context features, and transfer evidence

## Result Interpretation

### Transfer Potential Criteria

A gene is considered to have transfer potential if it meets any of the following conditions:

1. Located on a plasmid
2. Located within a prophage region
3. Located within an integron
4. Located within a transposon
5. Flanked by identical IS elements
6. Flanked by identical integrons

### Output File Format

Transfer potential analysis results include the following main fields:

- gene_id: Gene identifier
- gene_name: Gene name
- contig: Source contig
- location: Gene position
- mobile_potential: Transfer potential (Yes/No)
- mobile_evidence: Transfer evidence details

## Important Notes

1. Ensure all dependency tools are properly installed and accessible from command line
2. Verify database files are in correct locations
3. Input genome sequence must be in FASTA format
4. Absolute paths are recommended for input files and output directory
5. Analysis results are predictive and should be validated experimentally

## Support

For issues or suggestions:

1. Submit an issue on GitHub
2. Contact developers via email

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details
