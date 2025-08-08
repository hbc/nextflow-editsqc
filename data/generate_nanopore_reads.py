#!/usr/bin/env python
import random
import os
import sys
import gzip

def parse_gff(gff_file):
    """Parse GFF file and extract gene features"""
    genes = []
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 9 and fields[2] == 'gene':
                seqname = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8]
                
                # Extract gene ID
                gene_id = None
                for attr in attributes.split(';'):
                    if attr.startswith('ID='):
                        gene_id = attr.split('=')[1]
                        break
                
                genes.append({
                    'seqname': seqname,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'id': gene_id or f"gene_{start}_{end}"
                })
    return genes

def parse_fasta(fasta_file):
    """Parse FASTA file and return sequences"""
    sequences = {}
    current_seq = None
    current_name = None
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name and current_seq:
                    sequences[current_name] = current_seq
                current_name = line[1:]
                current_seq = ""
            else:
                current_seq += line
    
    if current_name and current_seq:
        sequences[current_name] = current_seq
    
    return sequences

def reverse_complement(seq):
    """Return reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def introduce_nanopore_errors(sequence, error_rate=0.08):
    """Introduce nanopore-like errors (insertions, deletions, substitutions)"""
    bases = ['A', 'T', 'G', 'C']
    result = []
    i = 0
    
    while i < len(sequence):
        if random.random() < error_rate:
            error_type = random.choice(['insertion', 'deletion', 'substitution'])
            
            if error_type == 'insertion':
                # Insert a random base
                result.append(random.choice(bases))
                result.append(sequence[i])
            elif error_type == 'deletion':
                # Skip this base (deletion)
                pass
            else:  # substitution
                result.append(random.choice(bases))
        else:
            result.append(sequence[i])
        i += 1
    
    return ''.join(result)

def generate_read_from_gene(gene_seq, gene_info, read_id, coverage_bias=1.0):
    """Generate a nanopore read from a gene sequence"""
    # Add some random flanking regions to simulate partial gene reads
    flank_left = random.randint(0, min(200, len(gene_seq) // 4))
    flank_right = random.randint(0, min(200, len(gene_seq) // 4))
    
    start_pos = max(0, flank_left)
    end_pos = min(len(gene_seq), len(gene_seq) - flank_right)
    
    if start_pos >= end_pos:
        start_pos = 0
        end_pos = len(gene_seq)
    
    read_seq = gene_seq[start_pos:end_pos]
    
    # Randomly decide strand (though genes have their own strand)
    if random.random() < 0.5:
        read_seq = reverse_complement(read_seq)
        strand = '-'
    else:
        strand = '+'
    
    # Introduce nanopore errors
    read_seq = introduce_nanopore_errors(read_seq)
    
    return read_seq, strand

def generate_nanopore_reads(genome_fasta, genome_gff, plasmid_fasta, plasmid_gff, num_reads=2000, output_file="nanopore_reads.fastq.gz"):
    """Generate nanopore reads from genes in GFF files"""
    
    # Parse genome data
    genome_sequences = parse_fasta(genome_fasta)
    genome_genes = parse_gff(genome_gff)
    
    # Parse plasmid data
    plasmid_sequences = parse_fasta(plasmid_fasta)
    plasmid_genes = parse_gff(plasmid_gff)
    
    # Combine all genes
    all_genes = []
    
    # Add genome genes
    for gene in genome_genes:
        seq_name = gene['seqname']
        if seq_name in genome_sequences:
            gene_seq = genome_sequences[seq_name][gene['start']-1:gene['end']]
            all_genes.append({
                'sequence': gene_seq,
                'info': gene,
                'source': 'genome'
            })
    
    # Add plasmid genes
    for gene in plasmid_genes:
        seq_name = gene['seqname']
        if seq_name in plasmid_sequences:
            gene_seq = plasmid_sequences[seq_name][gene['start']-1:gene['end']]
            all_genes.append({
                'sequence': gene_seq,
                'info': gene,
                'source': 'plasmid'
            })
    
    if not all_genes:
        print("No genes found in the provided GFF files!")
        return
    
    print(f"Found {len(all_genes)} genes to generate reads from:")
    for gene in all_genes:
        print(f"  - {gene['info']['id']} ({gene['source']}): {len(gene['sequence'])} bp")
    
    # Generate reads
    reads = []
    reads_per_gene = num_reads // len(all_genes)
    remaining_reads = num_reads % len(all_genes)
    
    read_id = 1
    
    for i, gene_data in enumerate(all_genes):
        gene_reads = reads_per_gene
        if i < remaining_reads:
            gene_reads += 1
        
        for j in range(gene_reads):
            read_seq, strand = generate_read_from_gene(
                gene_data['sequence'], 
                gene_data['info'], 
                read_id
            )
            
            # Create FASTQ entry
            read_name = f"read_{read_id:06d}_{gene_data['info']['id']}_{strand}"
            
            # Generate quality scores (simplified - using uniform quality)
            quality = ''.join(['I'] * len(read_seq))  # Quality score ~40
            
            reads.append({
                'name': read_name,
                'sequence': read_seq,
                'quality': quality
            })
            
            read_id += 1
    
    # Shuffle reads to mix them up
    random.shuffle(reads)
    
    # Write gzipped FASTQ file
    with gzip.open(output_file, 'wt') as f:
        for read in reads:
            f.write(f"@{read['name']}\n")
            f.write(f"{read['sequence']}\n")
            f.write("+\n")
            f.write(f"{read['quality']}\n")
    
    print(f"\nGenerated {len(reads)} nanopore reads")
    print(f"Output written to: {output_file}")
    
    # Print some statistics
    read_lengths = [len(read['sequence']) for read in reads]
    print("Read length statistics:")
    print(f"  - Min: {min(read_lengths)} bp")
    print(f"  - Max: {max(read_lengths)} bp")
    print(f"  - Mean: {sum(read_lengths)/len(read_lengths):.1f} bp")

if __name__ == "__main__":
    # Set random seed for reproducibility
    random.seed(42)
    
    # File paths
    genome_fasta = "./genome.fasta"
    genome_gff = "./genome.gff"
    plasmid_fasta = "./plasmid.fasta"
    plasmid_gff = "./plasmid.gff"
    output_file = "./nanopore_reads.fastq.gz"
    
    # Check if input files exist
    for file_path in [genome_fasta, genome_gff, plasmid_fasta, plasmid_gff]:
        if not os.path.exists(file_path):
            print(f"Error: File {file_path} not found!")
            print("Please run generate_example_data.py first to create the input files.")
            sys.exit(1)
    
    # Generate reads
    generate_nanopore_reads(
        genome_fasta=genome_fasta,
        genome_gff=genome_gff,
        plasmid_fasta=plasmid_fasta,
        plasmid_gff=plasmid_gff,
        num_reads=2000,
        output_file=output_file
    )
