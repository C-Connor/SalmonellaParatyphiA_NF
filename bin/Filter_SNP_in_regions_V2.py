#!/usr/bin/env python3

# This script filters SNPs by genomic location (phage, repeat_region, recombinant regions)
# and excludes them from a VCF file.
# Input: File with tab-separated start and end positions (no empty lines).
# Example:
# 6083    6085
# 16120   16465

from argparse import ArgumentParser

def parse_args():
    """Parse input arguments"""
    parser = ArgumentParser(description='Filter VCF to exclude SNPs in specified regions.')
    parser.add_argument('--vcf', type=str, required=True,
                        help='VCF file. (Required)')
    parser.add_argument('--position', type=str, required=True,
                        help='File with start and end positions. (Required)')
    parser.add_argument('--output', type=str, default='Region_filtered.vcf',
                        help='Output file name. Default: Region_filtered.vcf')
    return parser.parse_args()

def load_regions(position_file):
    """Load regions from position file and return list of (start, end) tuples"""
    regions = []
    with open(position_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue # skip blank lines
            parts = line.split('\t')
            if len(parts) != 2:
                continue  # skip malformed lines
            start, end = map(int, parts)
            regions.append((start, end))
    return regions

def is_position_excluded(pos, regions):
    """Check if the position falls within any of the excluded regions"""
    for start, end in regions:
        if start <= pos <= end:
            return True
    return False

def main():
    args = parse_args()
    regions = load_regions(args.position)

    data_lines = []
    exclude_lines = set()

    with open(args.vcf, 'r') as vcf, open(args.output, 'w') as output:
        for line in vcf:
            line = line.rstrip()
            if line.startswith("#"):
                output.write(line + '\n')
                continue

            data_lines.append(line)
            try:
                position = int(line.split("\t")[1])
            except (IndexError, ValueError):
                continue  # skip malformed lines

            if is_position_excluded(position, regions):
                exclude_lines.add(line)

        include_lines = [line for line in data_lines if line not in exclude_lines]
        output.write('\n'.join(include_lines) + '\n')

if __name__ == '__main__':
    main()
