#!/usr/bin/env python
import sys
import re
from collections import defaultdict


def parse_attributes(attr_str):
    attrs = {}
    for attr in attr_str.strip().split(';'):
        if '=' in attr:
            key, value = attr.split('=', 1)
            attrs[key.strip()] = value.strip()
    return attrs

def convert_gff_to_gtf(gff_file, output=sys.stdout, extension=50):
    id_counts = defaultdict(int)

    with open(gff_file, 'r') as infile:
        for line in infile:
            if line.startswith('#') or not line.strip():
                continue  # Skip comments or empty lines

            cols = line.strip().split('\t')
            if len(cols) != 9:
                print(f"Skipping malformed line {len(cols)}: {line.strip()}", file=sys.stderr)
                continue  # Skip malformed lines

            seqid, source, feature_type, start, end, score, strand, phase, attributes = cols

            valid_features = {'gene', 'CDS', 'transcription_unit', 'tRNA', 'rRNA'}
            if feature_type not in valid_features:
                print(f"Skipping invalid feature type: {feature_type}", file=sys.stderr)
                continue  # Skip invalid feature types
            feature_type = 'gene'
            attr_dict = parse_attributes(attributes)
            base_id = attr_dict.get('Name') or attr_dict.get('ID')
            if not base_id:
                print(f"Skipping line with no identifier: {line.strip()}", file=sys.stderr)
                continue  # Skip if no identifier found

            # Ensure unique ID
            id_counts[base_id] += 1
            sanitized_base_id = re.sub(r'[^a-zA-Z0-9]', '_', base_id)  # Replace non-alphanumeric characters with '_'
            unique_id = f"{sanitized_base_id}.{id_counts[base_id]}"

            # Extend start and end by given number of nt
            start = max(1, int(start) - extension)
            end = int(end) + extension

            gtf_attr = f'gene_id "{unique_id}"; transcript_id "TAC.{unique_id}";'

            # for feature in ['transcript', 'exon']:
            gtf_line = '\t'.join([
                seqid, source, feature_type, str(start), str(end),
                score, strand, '.', gtf_attr
            ])
            print(gtf_line, file=output)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Convert GFF to GTF with gene_id from Name or ID, and extend coordinates by 50nt.")
    parser.add_argument("gff_file", help="Input GFF file")
    parser.add_argument("-o", "--output", help="Output GTF file (default: stdout)", default=None)
    parser.add_argument("-e", "--extension", help="Number of nucleotides to extend at start/end", type=int, default=50)

    args = parser.parse_args()
    output = open(args.output, 'w') if args.output else sys.stdout

    convert_gff_to_gtf(args.gff_file, output=output, extension=args.extension)

    if output is not sys.stdout:
        output.close()
