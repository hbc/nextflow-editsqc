#!/usr/bin/env python3
"""
Replace gene_id with ref_gene_id in a GTF/GFF file when ref_gene_id is present.
Usage:
  replace_geneid_with_ref.py input.gtf > output.gtf
or
  replace_geneid_with_ref.py input.gtf output.gtf

The script preserves comments and other columns. It looks for the attribute
`ref_gene_id "..."` and, if found, replaces the value of `gene_id` with
that ref_gene_id. If no gene_id exists but ref_gene_id does, it will add
`gene_id "<ref>";` to the attributes.
"""
import sys
import re
from pathlib import Path

def process_line(line):
    if line.startswith('#'):
        return line
    cols = line.rstrip('\n').split('\t')
    if len(cols) < 9:
        return line
    attrs = cols[8]
    # find ref_gene_id
    m_ref = re.search(r'ref_gene_id\s+"([^"]+)"', attrs)
    if not m_ref:
        return line
    ref = m_ref.group(1)
    # replace existing gene_id if present
    if re.search(r'gene_id\s+"[^"]+"', attrs):
        attrs = re.sub(r'gene_id\s+"[^"]+"', f'gene_id "{ref}"', attrs)
    else:
        # append gene_id before final semicolon (or at end)
        attrs = attrs.strip()
        if not attrs.endswith(';'):
            attrs = attrs + ';'
        attrs = attrs + ' gene_id "{}";'.format(ref)
    cols[8] = attrs
    return '\t'.join(cols) + '\n'


def main(argv):
    if len(argv) < 2:
        print('Usage: replace_geneid_with_ref.py input.gtf [output.gtf]', file=sys.stderr)
        sys.exit(2)
    inp = Path(argv[1])
    out = None
    if len(argv) > 2:
        out = Path(argv[2])
    fh_in = inp.open('r')
    if out:
        fh_out = out.open('w')
    else:
        fh_out = sys.stdout
    try:
        for line in fh_in:
            fh_out.write(process_line(line))
    finally:
        fh_in.close()
        if out:
            fh_out.close()

if __name__ == '__main__':
    main(sys.argv)
