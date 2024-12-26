#!/usr/bin/env python3

import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=UserWarning)

import subprocess
import sys
import os
from pathlib import Path

def check_dependencies():
    """Check required software and Python packages"""
    # Check Python packages
    python_packages = ['pandas', 'biopython']
    missing_packages = []
    
    for package in python_packages:
        try:
            __import__(package)
        except ImportError:
            # Try to check package using conda list
            try:
                result = subprocess.run(['conda', 'list', package], capture_output=True, text=True)
                if package not in result.stdout:
                    missing_packages.append(package)
            except (subprocess.CalledProcessError, FileNotFoundError):
                missing_packages.append(package)
    
    if missing_packages:
        print(f"Error: Missing Python packages: {', '.join(missing_packages)}")
        print("Please install missing packages using conda:")
        print(f"conda install {' '.join(missing_packages)}")
        sys.exit(1)
    
    # Check conda tools
    conda_tools = ['blast', 'platon', 'phispy']
    missing_tools = []
    
    # Get conda environment list
    try:
        result = subprocess.run(['conda', 'env', 'list'], capture_output=True, text=True)
        env_list = result.stdout
        
        for tool in conda_tools:
            if tool not in env_list.lower():
                missing_tools.append(tool)
    except FileNotFoundError:
        print("Error: conda command not found. Please make sure conda is installed")
        sys.exit(1)
    
    if missing_tools:
        print(f"Error: Missing conda environments: {', '.join(missing_tools)}")
        print("Please create missing environments using conda:")
        for tool in missing_tools:
            print(f"conda create -n {tool} {tool}")
        sys.exit(1)
    
    print("All dependencies check passed")

# First check and install dependencies
check_dependencies()

# Then import other packages
import pandas as pd
from Bio import SeqIO
import json
import argparse

class GenomeAnalyzer:
    def __init__(self, input_genome, outdir, platon_db=None, threads=1):
        """Initialize GenomeAnalyzer
        
        Args:
            input_genome: Path to input genome file
            outdir: Path to output directory
            platon_db: Path to Platon database (optional)
            threads: Number of threads to use (default: 1)
        """
        self.input_genome = Path(input_genome)
        self.outdir = Path(outdir)
        self.threads = threads
        
        # Create output directory
        self.outdir.mkdir(parents=True, exist_ok=True)
        
        # Set up paths for output directories
        self.prokka_dir = self.outdir / 'prokka_results'
        self.platon_dir = self.outdir / 'platon_results'
        self.phispy_dir = self.outdir / 'phispy_results'
        
        # Set up paths for BLAST databases and output
        self.db_dir = Path(__file__).parent
        
        # Database paths
        self.subject_fasta_ISTnIN = self.db_dir / 'ISTnIN/istnin.fas'
        self.subject_fasta_RESfinder = self.db_dir / 'RESfinder/resfinder.fas'
        self.subject_fasta_VFDB = self.db_dir / 'VFDB/vfdb.fas'
        self.subject_fasta_ICE = self.db_dir / 'ICE/ICE.fas'
        # Output file paths
        self.output_file_ISTnIN = self.outdir / 'ISTnIN.tsv'
        self.output_file_RESfinder = self.outdir / 'RESfinder.tsv'
        self.output_file_VFDB = self.outdir / 'VFDB.tsv'
        self.output_file_ICE = self.outdir / 'ICE.tsv'
        
        # Platon database path
        self.platon_db = Path(platon_db) if platon_db else self.db_dir / 'platon'
        
        # Check if required files exist
        required_files = [
            self.subject_fasta_ISTnIN,
            self.subject_fasta_RESfinder,
            self.subject_fasta_VFDB
        ]
        
        for file in required_files:
            if not file.exists():
                raise FileNotFoundError(f"Required database file not found: {file}")

    def check_sequence_length(self, fasta_file):
        """Check sequence length in FASTA file"""
        min_length = float('inf')
        for record in SeqIO.parse(fasta_file, "fasta"):
            min_length = min(min_length, len(record.seq))
        return min_length

    def check_output_exists(self, output_path, description):
        """Check if output file exists
        
        Args:
            output_path: Output file path
            description: Analysis step description
            
        Returns:
            bool: True if file exists, False otherwise
        """
        if output_path.exists():
            print(f"{description} results already exist: {output_path}")
            return True
        return False

    def run_prokka(self):
        """Run Prokka annotation"""
        prefix = Path(self.input_genome).stem
        self.prokka_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if GBK file already exists
        prokka_gbk = self.prokka_dir / f"{prefix}.gbk"
        if self.check_output_exists(prokka_gbk, "Prokka"):
            return
        try:
            subprocess.run([
                'prokka',
                str(self.input_genome),
                '--outdir', str(self.prokka_dir),
                '--prefix', prefix,
                '--force',
                '--cpus', str(self.threads)
            ], check=True)
            print(f"Prokka results saved to: {self.prokka_dir}")
        except subprocess.CalledProcessError as e:
            print(f"Error running Prokka: {e}")

    def filter_blast_results(self, blast_result_file, query_fasta, subject_fasta, output_file):
        """Filter BLAST results
        
        Args:
            blast_result_file: BLAST output file
            query_fasta: Query sequence file
            subject_fasta: Subject sequence file
            output_file: Filtered output file
        """
        # Define BLAST output columns
        columns = [
            'query_id', 'subject_id', 'identity', 'alignment_length',
            'mismatches', 'gap_opens', 'query_start', 'query_end',
            'subject_start', 'subject_end', 'evalue', 'bit_score'
        ]
        
        # Read BLAST results
        blast_df = pd.read_csv(blast_result_file, sep='\t', names=columns)
        
        # Get query sequence lengths
        query_lengths = {}
        for record in SeqIO.parse(query_fasta, "fasta"):
            query_lengths[record.id] = len(record.seq)
            
        # Get subject sequence lengths
        subject_lengths = {}
        for record in SeqIO.parse(subject_fasta, "fasta"):
            subject_lengths[record.id] = len(record.seq)
            
        # Add sequence lengths to DataFrame
        blast_df['query_length'] = blast_df['query_id'].map(query_lengths)
        blast_df['subject_length'] = blast_df['subject_id'].map(subject_lengths)
        
        # Calculate coverage
        blast_df['query_coverage'] = (blast_df['alignment_length'] / blast_df['query_length']) * 100
        blast_df['subject_coverage'] = (blast_df['alignment_length'] / blast_df['subject_length']) * 100
        
        # Filter results - coverage and identity
        filtered_df = blast_df[
            (blast_df['subject_coverage'] >= 80) &
            (blast_df['identity'] >= 80)
        ]
        
        # Keep only highest identity result for each query_id
        filtered_df = filtered_df.sort_values('identity', ascending=False)
        filtered_df = filtered_df.loc[filtered_df.groupby('query_id')['identity'].idxmax()]
        
        # Select output columns
        output_columns = [
            'query_id', 'subject_id', 'identity', 'alignment_length',
            'query_start', 'query_end', 'query_coverage', 
            'subject_coverage',
            'evalue', 'bit_score'
        ]
        filtered_df = filtered_df[output_columns]
        
        # Sort by identity descending
        filtered_df = filtered_df.sort_values('identity', ascending=False)
        
        # Save filtered results
        filtered_df.to_csv(output_file, sep='\t', index=False)
        print(f"Filtered BLAST results saved to: {output_file}")
        
        return filtered_df

    def run_blast(self, query_fasta, data_base, blast_result_file):
        """Run BLAST analysis"""
        # Check if filtered BLAST results already exist
        filtered_output = Path(str(blast_result_file).replace('.tsv', '_filtered.tsv'))
        if self.check_output_exists(filtered_output, "BLAST"):
            return
            
        try:
            # Run BLAST
            subprocess.run([
                'blastn',
                '-query', query_fasta,
                '-db', data_base,
                '-evalue', '1e-5',
                '-out', blast_result_file,
                '-outfmt', '6',
                '-num_threads', str(self.threads)
            ], check=True)
            print(f"BLAST results saved to: {blast_result_file}")
            
            # Filter results
            self.filter_blast_results(blast_result_file, query_fasta, data_base, filtered_output)
            
        except subprocess.CalledProcessError as e:
            print(f"Error running BLAST: {e}")

    def run_platon(self):
        """Run Platon analysis"""
        if not self.platon_db.exists():
            print(f"Warning: Platon database not found at {self.platon_db}")
            print("Skipping Platon analysis")
            return
            
        prefix = Path(self.input_genome).stem
        self.platon_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            subprocess.run([
                'platon',
                str(self.input_genome),
                '--db', str(self.platon_db),
                '--output', str(self.platon_dir),
                '--prefix', prefix,
                '--threads', str(self.threads)
            ], check=True)
            print(f"Platon results saved to: {self.platon_dir}")
            
            # Check and parse JSON file
            json_file = self.platon_dir / f"{prefix}.json"
            if json_file.exists():
                self.parse_platon_json(
                    json_file,
                    self.outdir / 'platon_results_summary.tsv'
                )
            else:
                print(f"Warning: Platon JSON output file not found: {json_file}")
            
        except subprocess.CalledProcessError as e:
            print(f"Error running Platon: {e}")

    def run_phispy(self):
        """Run PhiSpy analysis"""
        prefix = Path(self.input_genome).stem
        self.phispy_dir.mkdir(parents=True, exist_ok=True)
        
        # Use Prokka's GBK file as input
        prokka_gbk = self.prokka_dir / f"{prefix}.gbk"
        if not prokka_gbk.exists():
            print(f"Error: Prokka GBK file not found: {prokka_gbk}")
            return
        
        try:
            subprocess.run([
                'phispy',
                str(prokka_gbk),
                '-o', str(self.phispy_dir),
                '-p', prefix,
                '--threads', str(self.threads)
            ], check=True)
            print(f"PhiSpy results saved to: {self.phispy_dir}")
            
            # Check and process output files
            prophage_files = list(self.phispy_dir.glob('*prophage_coordinates.tsv'))
            if prophage_files:
                prophage_file = prophage_files[0]  # Use first file found
                print(f"Found PhiSpy output file: {prophage_file}")
                
                # Read and process prophage predictions
                try:
                    # Read the original PhiSpy output
                    phage_df = pd.read_csv(prophage_file, sep='\t', header=None,
                        names=['Prophage number', 'contig', 'start', 'end', 'attL_start', 'attL_end',
                              'attR_start', 'attR_end', 'attL_seq', 'attR_seq', 'att_explanation'])
                    
                    # Select and reorder columns for summary
                    summary_df = phage_df[['contig', 'start', 'end', 'Prophage number']]
                    
                    # Save formatted version
                    output_file = self.outdir / 'phispy_results_summary.tsv'
                    summary_df.to_csv(output_file, sep='\t', index=False)
                    print(f"Prophage predictions saved to: {output_file}")
                except Exception as e:
                    print(f"Error processing PhiSpy output file: {e}")
            else:
                print("No PhiSpy output files found")
            
        except subprocess.CalledProcessError as e:
            print(f"Error running PhiSpy: {e}")

    def parse_gff_file(self, gff_file):
        """Parse GFF file to get gene location information
        
        Returns:
            dict: {
                gene_id: {
                    'contig': contig_id,
                    'start': start_pos,
                    'end': end_pos,
                    'strand': strand
                }
            }
        """
        gene_locations = {}
        try:
            with open(gff_file) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue
                        
                    if parts[2] == 'CDS':  # Only process CDS entries
                        contig = parts[0]
                        start = int(parts[3])
                        end = int(parts[4])
                        strand = parts[6]
                        
                        # Extract gene ID from attributes field
                        attributes = dict(item.split('=') for item in parts[8].split(';') if '=' in item)
                        if 'ID' in attributes:
                            gene_id = attributes['ID']
                            gene_locations[gene_id] = {
                                'contig': contig,
                                'start': start,
                                'end': end,
                                'strand': strand
                            }
            
            print(f"Extracted location information for {len(gene_locations)} genes from GFF file")
        except Exception as e:
            print(f"Error parsing GFF file: {e}")
        
        return gene_locations

    def add_genome_locations(self, blast_file, gff_file, output_file):
        """Add genome location information to BLAST results"""
        # Read BLAST results
        blast_df = pd.read_csv(blast_file, sep='\t')
        
        # Extract gene location information from GFF file
        gene_locations = self.parse_gff_file(gff_file)
        
        # Add location information to DataFrame
        def get_location(query_id):
            # Handle Prokka gene ID format (if needed)
            gene_id = query_id
            if not gene_id in gene_locations:
                # Try adding 'CDS_' prefix (may be used in Prokka GFF files)
                gene_id = f"CDS_{query_id}"
            
            loc = gene_locations.get(gene_id, {})
            if not loc:
                print(f"Warning: No location found for {query_id}")
            return loc
        
        # Add location information
        blast_df['contig'] = blast_df['query_id'].map(lambda x: get_location(x).get('contig'))
        blast_df['genome_start'] = blast_df['query_id'].map(lambda x: get_location(x).get('start'))
        blast_df['genome_end'] = blast_df['query_id'].map(lambda x: get_location(x).get('end'))
        blast_df['strand'] = blast_df['query_id'].map(lambda x: get_location(x).get('strand'))
        

            


        columns = [
                'query_id', 'subject_id', 'identity', 'alignment_length',
                'query_start', 'query_end', 'query_coverage', 'subject_coverage',
                'contig', 'genome_start', 'genome_end', 'strand',
                'evalue', 'bit_score'
            ]
        
        result_df = blast_df[columns]
        
        # Check if location information was successfully added
        location_added = result_df['genome_start'].notna().any()
        if not location_added:
            print("Warning: Failed to add location information!")
            print("Sample records:")
            print(result_df.head())
            print("\nQuery IDs from original BLAST file:")
            print(blast_df['query_id'].head())
            print("\nExample gene locations extracted:")
            print(list(gene_locations.keys())[:5])
        
        # Save results
        result_df.to_csv(output_file, sep='\t', index=False)
        print(f"Added genome location information and saved to: {output_file}")
        return result_df
    
    def analyze_resistance_genes_context(self):
        """Analyze genomic context of resistance genes"""
        try:
            # Read all necessary files
            res_file = self.outdir / 'RESfinder_filtered_with_locations.tsv'
            istnin_file = self.outdir / 'ISTnIN_filtered_with_locations.tsv'
            phage_file = self.outdir / 'phispy_results_summary.tsv'
            plasmid_file = self.outdir / 'platon_results_summary.tsv'
            ice_file = self.outdir / 'ICE_filtered_with_locations.tsv'
            
            if not all(f.exists() for f in [res_file, istnin_file]):
                print("Missing required input files")
                return
            
            # Read files
            res_df = pd.read_csv(res_file, sep='\t') if res_file.exists() else pd.DataFrame()
            istnin_df = pd.read_csv(istnin_file, sep='\t') if istnin_file.exists() else pd.DataFrame()
            ice_df = pd.read_csv(ice_file, sep='\t') if ice_file.exists() else pd.DataFrame()
            plasmid_df = pd.read_csv(plasmid_file, sep='\t') if plasmid_file.exists() else pd.DataFrame()
            # Read phage file if exists and standardize column names
            if (phage_file.exists()):
                phage_df = pd.read_csv(phage_file, sep='\t')
                # Rename columns to match expected names
                column_mapping = {
                    'contig': 'Contig',  # if column is 'contig'
                    'start': 'Start',    # if column is 'start'
                    'end': 'End'         # if column is 'end'
                }
                phage_df = phage_df.rename(columns={v: k for k, v in column_mapping.items()})
            else:
                phage_df = pd.DataFrame()
               
            
            # Create results list
            results = []
            
            # Analyze each resistance gene
            for _, res_gene in res_df.iterrows():
                result = {
                    'resistance_gene': res_gene['query_id'],
                    'resistance_gene_name': res_gene['subject_id'],
                    'contig': res_gene['contig'],
                    'start': res_gene['genome_start'],
                    'end': res_gene['genome_end'],
                    'identity': res_gene['identity'],
                    'on_plasmid': 'No',
                    'on_phage': 'No',
                    'in_integron': 'No',
                    'in_transposon': 'No',
                    'flanking_IS': 'No',
                    'flanking_IN': 'No',
                    'mobile_potential': 'No',
                    'mobile_evidence': []
                }
                
                # Check if on plasmid
                if not plasmid_df.empty:
                    for _, plasmid in plasmid_df.iterrows():
                        if (res_gene['contig'] == plasmid['contig_id'] and
                            res_gene['genome_start'] >= plasmid['hit_contig_start'] and
                            res_gene['genome_end'] <= plasmid['hit_contig_end']):
                            result['on_plasmid'] = 'Yes'
                            result['mobile_evidence'].append(f"Located on plasmid {plasmid['hit_plasmid_id']}")
                            break
                
                # Check if in prophage region
                if not phage_df.empty:
                    for _, phage in phage_df.iterrows():
                        if (res_gene['contig'] == phage['contig'] and
                            res_gene['genome_start'] >= phage['start'] and
                            res_gene['genome_end'] <= phage['end']):
                            result['on_phage'] = 'Yes'
                            result['mobile_evidence'].append(f"Located in prophage region {phage.get('Prophage number', 'unknown')}")
                            break

                # Check if within transposon range
                for _, tn in istnin_df[istnin_df['element_type'] == 'Tn'].iterrows():
                    if (res_gene['contig'] == tn['contig'] and
                        res_gene['genome_start'] >= tn['genome_start'] and
                        res_gene['genome_end'] <= tn['genome_end']):
                        result['in_transposon'] = 'Yes'
                        result['mobile_evidence'].append(f"Located within transposon {tn['subject_id']}")
                        break

                # Check if within ICE range
                for _, ice in ice_df.iterrows():
                    if (res_gene['contig'] == ice['contig'] and
                        res_gene['genome_start'] >= ice['genome_start'] and
                        res_gene['genome_end'] <= ice['genome_end']):
                        result['in_transposon'] = 'Yes'
                        result['mobile_evidence'].append(f"Located within ICE {ice['subject_id']}")
                        break
                
                # Check if within integron range
                for _, integron in istnin_df[istnin_df['element_type'] == 'IN'].iterrows():
                    if (res_gene['contig'] == integron['contig'] and
                        res_gene['genome_start'] >= integron['genome_start'] and
                        res_gene['genome_end'] <= integron['genome_end']):
                        result['in_integron'] = 'Yes'
                        result['mobile_evidence'].append(f"Located within integron {integron['subject_id']}")
                        break
                
                # Check mobile elements within 10kb range
                window = 10000
                
                # Check for same IN on both sides
                # Check upstream IN within 10kb
                left_in = istnin_df[
                    (istnin_df['element_type'] == 'IN') &
                    (istnin_df['contig'] == res_gene['contig']) &
                    (istnin_df['genome_end'] <= res_gene['genome_start']) &
                    (istnin_df['genome_end'] >= res_gene['genome_start'] - window)
                ]
                
                # Check downstream IN within 10kb
                right_in = istnin_df[
                    (istnin_df['element_type'] == 'IN') &
                    (istnin_df['contig'] == res_gene['contig']) &
                    (istnin_df['genome_start'] >= res_gene['genome_end']) &
                    (istnin_df['genome_start'] <= res_gene['genome_end'] + window)
                ]
                
                # Get upstream and downstream IN types
                left_in_types = set(left_in['subject_id'])
                right_in_types = set(right_in['subject_id'])
                
                # Check for same IN
                common_in = left_in_types.intersection(right_in_types)
                if common_in:
                    result['flanking_IN'] = 'Yes'
                    result['mobile_evidence'].append(f"Same integron found on both sides: {';'.join(common_in)}")
                    result['mobile_potential'] = 'Yes'
                
                # Check for IS elements within 10kb range
                # Check upstream IS within 10kb
                left_is = istnin_df[
                    (istnin_df['element_type'] == 'IS') &
                    (istnin_df['contig'] == res_gene['contig']) &
                    (istnin_df['genome_end'] <= res_gene['genome_start']) &
                    (istnin_df['genome_end'] >= res_gene['genome_start'] - window)
                ]
                
                # Check downstream IS within 10kb
                right_is = istnin_df[
                    (istnin_df['element_type'] == 'IS') &
                    (istnin_df['contig'] == res_gene['contig']) &
                    (istnin_df['genome_start'] >= res_gene['genome_end']) &
                    (istnin_df['genome_start'] <= res_gene['genome_end'] + window)
                ]
                
                # Get upstream and downstream IS types
                left_is_types = set(left_is['subject_id'])
                right_is_types = set(right_is['subject_id'])
                
                # Check for same IS
                common_is = left_is_types.intersection(right_is_types)
                if common_is:
                    result['flanking_IS'] = 'Yes'
                    result['mobile_evidence'].append(f"Same IS sequence found on both sides: {';'.join(common_is)}")
                    result['mobile_potential'] = 'Yes'
                
                # Determine resistance gene mobility potential
                if (result['on_plasmid'] == 'Yes' or 
                    result['on_phage'] == 'Yes' or 
                    result['in_integron'] == 'Yes' or 
                    result['in_transposon'] == 'Yes' or 
                    result['flanking_IS'] == 'Yes' or 
                    result['flanking_IN'] == 'Yes'):
                    result['mobile_potential'] = 'Yes'
                
                results.append(result)
            
            # Convert to DataFrame
            results_df = pd.DataFrame(results)
            
            # Convert list to string
            results_df['mobile_evidence'] = results_df['mobile_evidence'].apply(lambda x: ' | '.join(x) if x else '')
            
            # Save results
            output_file = self.outdir / 'resistance_genes_context.tsv'
            results_df.to_csv(output_file, sep='\t', index=False)
            print(f"Resistance gene context analysis results saved to: {output_file}")
            
            return results_df
            
        except Exception as e:
            print(f"Error analyzing resistance gene context: {e}")
            return None

    def analyze_virulence_genes_context(self):
        """Analyze mobility potential of virulence genes"""
        try:
            # Read required files
            vfdb_file = self.outdir / 'VFDB_filtered_with_locations.tsv'
            istnin_file = self.outdir / 'ISTnIN_filtered_with_locations.tsv'
            platon_file = self.outdir / 'platon_results_summary.tsv'
            phispy_file = self.outdir / 'prophage_coordinates.tsv'
            vfdb_parsed_file = self.db_dir / 'VFDB/vfdb_parsed.tsv'  
            ice_file = self.outdir / 'ICE_filtered_with_locations.tsv'

            if not vfdb_file.exists():
                print("VFDB analysis result file not found")
                return None

            # Read data
            vfdb_df = pd.read_csv(vfdb_file, sep='\t')
            istnin_df = pd.read_csv(istnin_file, sep='\t') if istnin_file.exists() else pd.DataFrame()
            plasmid_df = pd.read_csv(platon_file, sep='\t') if platon_file.exists() else pd.DataFrame()
            ice_df = pd.read_csv(ice_file, sep='\t') if ice_file.exists() else pd.DataFrame()
            
            # Read vfdb_parsed.tsv and create ID to Gene_Name mapping
            if vfdb_parsed_file.exists():
                vfdb_parsed_df = pd.read_csv(vfdb_parsed_file, sep='\t')
                id_to_gene = dict(zip(vfdb_parsed_df['VFG_ID'], vfdb_parsed_df['Gene_Name']))
            else:
                print("Warning: VFDB parsed file not found, will use original subject_id")
                id_to_gene = {}
            
            # Results list
            results = []
            
            # Analyze each virulence gene
            for _, vf_gene in vfdb_df.iterrows():
                # Get Gene_Name, use original subject_id if not found
                gene_name = id_to_gene.get(vf_gene['subject_id'], vf_gene['subject_id'])
                
                result = {
                    'gene_id': vf_gene['query_id'],
                    'gene_name': gene_name,   # Use mapped gene name
                    'contig': vf_gene['contig'],
                    'start': vf_gene['genome_start'],
                    'end': vf_gene['genome_end'],
                    'on_plasmid': 'No',
                    'on_phage': 'No',
                    'in_integron': 'No', 
                    'in_transposon': 'No',
                    'flanking_IS': 'No',
                    'flanking_IN': 'No',
                    'mobile_potential': 'No',
                    'mobile_evidence': []
                }

                # Check if on plasmid
                if not plasmid_df.empty:
                    for _, plasmid in plasmid_df.iterrows():
                        if (vf_gene['contig'] == plasmid['contig_id'] and
                            vf_gene['genome_start'] >= plasmid['hit_contig_start'] and
                            vf_gene['genome_end'] <= plasmid['hit_contig_end']):
                            result['on_plasmid'] = 'Yes'
                            result['mobile_evidence'].append(f"Located on plasmid {plasmid['hit_plasmid_id']}")
                            break

                # Check if on prophage
                if phispy_file.exists():
                    phage_regions = pd.read_csv(phispy_file, sep='\t')
                    in_phage = False
                    for _, phage in phage_regions.iterrows():
                        if (phage['contig'] == vf_gene['contig'] and
                            phage['start'] <= vf_gene['genome_start'] <= phage['end']):
                            in_phage = True
                            break
                    if in_phage:
                        result['on_phage'] = 'Yes'
                        result['mobile_evidence'].append(f"Located in prophage region {phage.get('Prophage number', 'unknown')}")
              
                
                # Check if on Tn
                for _, tn in istnin_df[istnin_df['element_type'] == 'Tn'].iterrows():
                    if (vf_gene['contig'] == tn['contig'] and
                        vf_gene['genome_start'] >= tn['genome_start'] and
                        vf_gene['genome_end'] <= tn['genome_end']):
                        result['in_transposon'] = 'Yes'
                        result['mobile_evidence'].append(f"Located within transposon {tn['subject_id']}")
                        break
                
                # Check if on ice
                for _, ice in ice_df.iterrows():
                    if (vf_gene['contig'] == ice['contig'] and
                        vf_gene['genome_start'] >= ice['genome_start'] and
                        vf_gene['genome_end'] <= ice['genome_end']):
                        result['in_transposon'] = 'Yes'
                        result['mobile_evidence'].append(f"Located within ICE {ice['subject_id']}")
                        break
                
                # Check if on IN
                for _, integron in istnin_df[istnin_df['element_type'] == 'IN'].iterrows():
                    if (vf_gene['contig'] == integron['contig'] and
                        vf_gene['genome_start'] >= integron['genome_start'] and
                        vf_gene['genome_end'] <= integron['genome_end']):
                        result['in_integron'] = 'Yes'
                        result['mobile_evidence'].append(f"Located within integron {integron['subject_id']}")
                        break

                window = 10000
                # Check for IS elements within 10kb upstream
                left_is = istnin_df[
                    (istnin_df['element_type'] == 'IS') &
                    (istnin_df['contig'] == vf_gene['contig']) &
                    (istnin_df['genome_end'] <= vf_gene['genome_start']) &
                    (istnin_df['genome_end'] >= vf_gene['genome_start'] - window)
                ]
                
                # Check for IS elements within 10kb downstream
                right_is = istnin_df[
                    (istnin_df['element_type'] == 'IS') &
                    (istnin_df['contig'] == vf_gene['contig']) &
                    (istnin_df['genome_start'] >= vf_gene['genome_end']) &
                    (istnin_df['genome_start'] <= vf_gene['genome_end'] + window)
                ]
                
                # Get IS types from both sides
                left_is_types = set(left_is['subject_id'])
                right_is_types = set(right_is['subject_id'])
                
                # Check for common IS elements
                common_is = left_is_types.intersection(right_is_types)
                if common_is:
                    result['flanking_IS'] = 'Yes'
                    result['mobile_evidence'].append(f"Same IS sequence found on both sides: {';'.join(common_is)}")
                    result['mobile_potential'] = 'Yes'

                # Check for integrons within 10kb upstream
                left_in = istnin_df[
                    (istnin_df['element_type'] == 'IN') &
                    (istnin_df['contig'] == vf_gene['contig']) &
                    (istnin_df['genome_end'] <= vf_gene['genome_start']) &
                    (istnin_df['genome_end'] >= vf_gene['genome_start'] - window)
                ]

                # Check for integrons within 10kb downstream
                right_in = istnin_df[
                    (istnin_df['element_type'] == 'IN') &
                    (istnin_df['contig'] == vf_gene['contig']) &
                    (istnin_df['genome_start'] >= vf_gene['genome_end']) &
                    (istnin_df['genome_start'] <= vf_gene['genome_end'] + window)
                ]

                # Get integron types from both sides
                left_in_types = set(left_in['subject_id'])
                right_in_types = set(right_in['subject_id'])

                # Check for common integrons
                common_in = left_in_types.intersection(right_in_types)
                if common_in:
                    result['flanking_IN'] = 'Yes'
                    result['mobile_evidence'].append(f"Same integron found on both sides: {';'.join(common_in)}")
                    result['mobile_potential'] = 'Yes'
                
                # Determine mobility potential of virulence gene
                if (result['on_plasmid'] == 'Yes' or 
                    result['on_phage'] == 'Yes' or 
                    result['in_integron'] == 'Yes' or 
                    result['in_transposon'] == 'Yes' or 
                    result['flanking_IS'] == 'Yes' or 
                    result['flanking_IN'] == 'Yes'):
                    result['mobile_potential'] = 'Yes'
                
                results.append(result)
            
            # Convert to DataFrame
            results_df = pd.DataFrame(results)
            
            # Convert lists to strings
            results_df['mobile_evidence'] = results_df['mobile_evidence'].apply(lambda x: ' | '.join(x) if x else '')
            
            # Set column order
            columns = ['gene_id', 'gene_name', 'contig', 'start', 'end', 
                      'on_plasmid', 'on_phage', 'in_integron', 'in_transposon',
                      'flanking_IS', 'flanking_IN', 'mobile_potential', 'mobile_evidence']
            results_df = results_df[columns]
            
            # Save results
            output_file = self.outdir / 'virulence_genes_context.tsv'
            results_df.to_csv(output_file, sep='\t', index=False)
            print(f"Virulence gene context analysis results saved to: {output_file}")
            
            return results_df
            
        except Exception as e:
            print(f"Error analyzing virulence gene context: {e}")
            return None

    def run_analysis(self):
        """Run complete analysis pipeline"""
        print(f"Starting genome analysis: {self.input_genome}")
        print(f"Output directory: {self.outdir}")
        
        # First run Prokka, as PhiSpy needs its output
        self.run_prokka()
        
        # Run BLAST analysis
        prokka_ffn = self.prokka_dir / f"{Path(self.input_genome).stem}.ffn"
        prokka_gff = self.prokka_dir / f"{Path(self.input_genome).stem}.gff"
        prokka_fna = self.prokka_dir / f"{Path(self.input_genome).stem}.fna"
        
        if prokka_ffn.exists() and prokka_gff.exists() and prokka_fna.exists():
            # Run BLAST and filter
            self.run_blast(str(prokka_fna), str(self.subject_fasta_ISTnIN), str(self.output_file_ISTnIN))
            self.run_blast(str(prokka_ffn), str(self.subject_fasta_RESfinder), str(self.output_file_RESfinder))
            self.run_blast(str(prokka_ffn), str(self.subject_fasta_VFDB), str(self.output_file_VFDB))
            self.run_blast(str(prokka_fna), str(self.subject_fasta_ICE), str(self.output_file_ICE))

            # Add genome location information
            istnin_filtered = Path(str(self.output_file_ISTnIN).replace('.tsv', '_filtered.tsv'))
            resfinder_filtered = Path(str(self.output_file_RESfinder).replace('.tsv', '_filtered.tsv'))
            vfdb_filtered = Path(str(self.output_file_VFDB).replace('.tsv', '_filtered.tsv'))
            ICE_filtered = Path(str(self.output_file_ICE).replace('.tsv', '_filtered.tsv'))

            if istnin_filtered.exists():
                self.add_location(
                    istnin_filtered,
                    self.outdir / 'ISTnIN_filtered_with_locations.tsv'
                )
            
            if resfinder_filtered.exists():
                self.add_genome_locations(
                    resfinder_filtered,
                    prokka_gff,
                    self.outdir / 'RESfinder_filtered_with_locations.tsv'
                )
            
            if vfdb_filtered.exists():
                self.add_genome_locations(
                    vfdb_filtered,
                    prokka_gff,
                    self.outdir / 'VFDB_filtered_with_locations.tsv'
                )
            
            if ICE_filtered.exists():
                self.add_location(
                    ICE_filtered,
                    self.outdir / 'ICE_filtered_with_locations.tsv'
                )
        else:
            print("Error: Prokka output files not found")
        
        # Run other analyses
        print("\nRunning Platon analysis...")
        self.run_platon()
        
        print("\nRunning PhiSpy analysis...")
        self.run_phispy()
        
        # After all analyses, analyze resistance gene context
        print("\nAnalyzing resistance gene context...")
        self.analyze_resistance_genes_context()
        print("\nAnalyzing virulence gene context...")
        self.analyze_virulence_genes_context()

    def parse_platon_json(self, json_file, output_file):
        """Parse Platon JSON file and extract plasmid information
        
        Args:
            json_file: Path to Platon generated JSON file
            output_file: Path to output TSV file
        """
        try:
            # Read JSON file
            with open(json_file) as f:
                data = json.load(f)
            
            # Store all extracted information
            plasmid_info = []
            
            # Process each contig
            for contig_id, contig_data in data.items():
                # Basic information
                basic_info = {
                    'contig_id': contig_id,
                    'length': contig_data.get('length', ''),
                    'gc': contig_data.get('gc', ''),
                    'coverage': contig_data.get('coverage', ''),
                    'circular': contig_data.get('circular', False),
                    'plasmid_probability': contig_data.get('plasmid_probability', '')
                }
                
                # Extract plasmid hits information
                plasmid_hits = contig_data.get('plasmid_hits', [])
                for hit in plasmid_hits:
                    hit_info = basic_info.copy()
                    
                    # Add hit specific information
                    if 'plasmid' in hit:
                        hit_info.update({
                            'hit_plasmid_id': hit['plasmid'].get('id', ''),
                            'hit_plasmid_description': hit['plasmid'].get('description', ''),
                            'hit_plasmid_length': hit['plasmid'].get('length', ''),
                            'hit_identity': hit.get('identity', ''),
                            'hit_coverage': hit.get('coverage', ''),
                            'hit_contig_start': hit.get('contig_start', ''),
                            'hit_contig_end': hit.get('contig_end', ''),
                            'hit_plasmid_start': hit.get('plasmid_start', ''),
                            'hit_plasmid_end': hit.get('plasmid_end', '')
                        })
                    
                    plasmid_info.append(hit_info)
            
            # If no plasmid information found
            if not plasmid_info:
                print("No plasmid information found in JSON file")
                return None
            
            # Convert to DataFrame
            df = pd.DataFrame(plasmid_info)
            
            # Set column order
            columns = [
                'contig_id', 'length', 'gc', 'coverage', 'circular', 
                'plasmid_probability', 'hit_plasmid_id', 'hit_plasmid_description',
                'hit_plasmid_length', 'hit_identity', 'hit_coverage',
                'hit_contig_start', 'hit_contig_end',
                'hit_plasmid_start', 'hit_plasmid_end'
            ]
            df = df[columns]
            
            # Save as TSV file
            df.to_csv(output_file, sep='\t', index=False)
            print(f"Plasmid information saved to: {output_file}")
            
            return df
            
        except Exception as e:
            print(f"Error parsing Platon JSON file: {e}")
            return None

    def add_location(self, input_file, output_file, is_istnin=False):
        """Add columns to ICE_filtered.tsv and ISTnIN_filtered.tsv"""
        # Read the input file
        df = pd.read_csv(input_file, sep='\t')
        
        # Add new columns
        df['contig'] = df['query_id']
        df['genome_start'] = df['query_start']
        df['genome_end'] = df['query_end']
        df['strand'] = ''
        
        # Reorder columns
        columns = list(df.columns)
        insert_index = columns.index('subject_coverage') + 1
        new_columns = columns[:insert_index] + ['contig', 'genome_start', 'genome_end', 'strand'] + columns[insert_index:]
        df = df[new_columns]

        # Add type marker if ISTnIN result
        if is_istnin:
            def get_element_type(subject_id):
                if 'ISFinder' in subject_id:
                    return 'IS'
                elif 'INTEGRALL' in subject_id:
                    return 'IN'
                else:
                    return 'Tn'
            
            df['element_type'] = df['subject_id'].apply(get_element_type)
            
            # Update column order to include element_type
            columns = [
                'query_id', 'subject_id', 'identity', 'alignment_length',
                'query_start', 'query_end', 'query_coverage', 'subject_coverage',
                'contig', 'genome_start', 'genome_end', 'strand',
                'evalue', 'bit_score', 'element_type'
            ]
        else:
            columns = [
                'query_id', 'subject_id', 'identity', 'alignment_length',
                'query_start', 'query_end', 'query_coverage', 'subject_coverage',
                'contig', 'genome_start', 'genome_end', 'strand',
                'evalue', 'bit_score'
            ]
        
        # Save the updated DataFrame to the output file
        df.to_csv(output_file, sep='\t', index=False)
        print(f"Updated tsv saved to: {output_file}")

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Genome analysis pipeline')
    parser.add_argument('input_genome', help='Input genome file path')
    parser.add_argument('--outdir', default='../result', help='Output directory path (default: ../result)')
    parser.add_argument('--db', help='Platon database path (default: ../RESpt/platon)')
    parser.add_argument('--threads', '-t', type=int, default=4, help='Number of threads (default: 4)')
    return parser.parse_args()

def main():
    # Check dependencies first
    check_dependencies()
    
    args = parse_args()
    analyzer = GenomeAnalyzer(args.input_genome, args.outdir, args.db, args.threads)
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
