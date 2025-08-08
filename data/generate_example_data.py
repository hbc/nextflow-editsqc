#!/usr/bin/env python
import random
import os

def generate_random_sequence(length):
    """Generate a random DNA sequence of given length"""
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choice(bases) for _ in range(length))

def write_fasta(filename, sequences):
    """Write sequences to FASTA file"""
    with open(filename, 'w') as f:
        for seq_id, sequence in sequences.items():
            f.write(f">{seq_id}\n")
            # Write sequence in 80-character lines
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")

def write_gff(filename, features):
    """Write features to GFF file"""
    with open(filename, 'w') as f:
        f.write("##gff-version 3\n")
        for feature in features:
            f.write("\t".join(map(str, feature)) + "\n")

# Set random seed for reproducibility
random.seed(42)

# Create data directory
# os.makedirs("data", exist_ok=True)

# Generate genome sequence (20kb with 5 genes of 1kb separated by 1kb intergenic regions)
genome_sequence = ""
genome_features = []

# Pattern: 1kb intergenic + 1kb gene (repeated 5 times) + 1kb intergenic
current_pos = 1

for i in range(5):
    # Add 1kb intergenic region
    intergenic = generate_random_sequence(1000)
    genome_sequence += intergenic
    current_pos += 1000
    
    # Add 1kb gene
    gene_seq = generate_random_sequence(1000)
    genome_sequence += gene_seq
    
    # Add gene feature to GFF
    gene_start = current_pos
    gene_end = current_pos + 999
    gene_id = f"gene{i+1:02d}"
    
    genome_features.append([
        "chromosome", "auto", "gene", gene_start, gene_end, ".", "+", ".", 
        f"ID={gene_id};Name={gene_id}"
    ])
    
    current_pos += 1000

# Add final 1kb intergenic region to reach 20kb total
final_intergenic = generate_random_sequence(8000)  # 8kb to reach exactly 20kb
genome_sequence += final_intergenic

# Write genome files
genome_sequences = {"chromosome": genome_sequence}
write_fasta("./genome.fasta", genome_sequences)
write_gff("./genome.gff", genome_features)

# Generate plasmid sequence (5kb with promoter and gene)
plasmid_features = []
plasmid_sequence = ""

# Add 2kb intergenic region
intergenic1 = generate_random_sequence(2000)
plasmid_sequence += intergenic1

# Add 500bp promoter
promoter_seq = generate_random_sequence(500)
plasmid_sequence += promoter_seq
plasmid_features.append([
    "plasmid", "example", "promoter", 2001, 2500, ".", "+", ".", 
    "ID=promoter01;Name=promoter01"
])

# Add 1kb gene
gene_seq = generate_random_sequence(1000)
plasmid_sequence += gene_seq
plasmid_features.append([
    "plasmid", "auto", "gene", 2501, 3500, ".", "+", ".", 
    "ID=plasmid_gene01;Name=plasmid_gene01"
])

# Add remaining sequence to reach 5kb
remaining = generate_random_sequence(1500)
plasmid_sequence += remaining

# Write plasmid files
plasmid_sequences = {"plasmid": plasmid_sequence}
write_fasta("./plasmid.fasta", plasmid_sequences)
write_gff("./plasmid.gff", plasmid_features)

print("Generated example ./:")
print(f"- Genome: {len(genome_sequence)} bp with {len(genome_features)} genes")
print(f"- Plasmid: {len(plasmid_sequence)} bp with {len(plasmid_features)} features")
print("Files created in ./ directory:")
print("  - genome.fasta, genome.gff")
print("  - plasmid.fasta, plasmid.gff")
