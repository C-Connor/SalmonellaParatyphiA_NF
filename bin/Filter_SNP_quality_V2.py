#!/usr/bin/env python3

# This script filters SNPs in a VCF file by Phred-scaled quality score.
# Optionally, it can filter by ratio of ALT alleles from high-quality reads (currently commented out).

from argparse import ArgumentParser

def parse_args():
    """Parse the input arguments, use '-h' for help"""
    parser = ArgumentParser(description='Filter VCF by Phred-scaled quality score (QUAL).')
    parser.add_argument('--vcf', type=str, required=True,
                        help='Mapped raw VCF file, with DP4 scores. (Required)')
    parser.add_argument('--phr', type=int, default=20,
                        help='Minimum Phred-scaled quality score. Default: 20.')
    parser.add_argument('--output', type=str, default='Phred_Filtered.vcf',
                        help='Output VCF file. Default: Phred_Filtered.vcf')
    return parser.parse_args()

def main():
    args = parse_args()
    phrdQ = float(args.phr)

    with open(args.vcf, 'r') as vcf, open(args.output, 'w') as output:
        for line in vcf:
            if line.startswith('#'):
                output.write(line)
                continue

            fields = line.strip().split('\t')

            try:
                phred_quality = float(fields[5])
                info_field = fields[7]

                # Extract MQ and DP4 values
                mq = float(info_field.split('MQ=')[1].split(';')[0])
                dp4 = info_field.split('DP4=')[1].split(';')[0].split(',')

                dp4_counts = list(map(int, dp4))
                total_reads = sum(dp4_counts)
                alt_reads = dp4_counts[2] + dp4_counts[3]

                # Ratio of ALT reads from high-quality reads
                read_ratio = float(alt_reads) / total_reads if total_reads > 0 else 0

                # Filtering logic
                if len(fields[4]) == 1 and phred_quality > phrdQ:
                    output.write(line)

            except Exception as e:
                # Skip malformed entries
                continue

if __name__ == '__main__':
    main()
