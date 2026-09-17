#!/usr/bin/env python3
"""Extract GRCh37 Ensembl-87 transcript TSSs matching a gene BED."""

import argparse
from nucleosuite.gene_sets import read_genes
from nucleosuite.transcript_tss import ENSEMBL_GTF_URL, extract_transcript_tss


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gtf', required=True, help='Ensembl GRCh37 release-87 GTF or GTF.gz')
    parser.add_argument('--genes-bed', required=True)
    parser.add_argument('--out', default='hg19_ensembl87_transcript_tss.tsv.gz')
    args = parser.parse_args()
    transcripts, genes = extract_transcript_tss(args.gtf, read_genes(args.genes_bed), args.out)
    print(f'Wrote {transcripts:,} transcript TSSs for {genes:,} genes: {args.out}')
    print(f'Source: {ENSEMBL_GTF_URL}')


if __name__ == '__main__':
    main()
