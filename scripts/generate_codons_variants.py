#!/usr/bin/env python3
import sys
import csv
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def read_codon_table(codon_table_file):
    """read codon table and return codon->aa mapping"""
    with open(codon_table_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        return {row['codon']: row['aa'] for row in reader}

def read_fasta(fasta_file, target_id):
    """read sequence from fasta file by identifier"""
    sequences = SeqIO.parse(fasta_file, "fasta")
    
    for record in sequences:
        if record.id == target_id:
            return record.id, str(record.seq)
    
    # sequence not found
    print(f"error: sequence with identifier '{target_id}' not found in fasta file", file=sys.stderr)
    sys.exit(1)

def generate_all_codons():
    """generate all 64 possible codons"""
    bases = ['A', 'T', 'G', 'C']
    codons = []
    for b1 in bases:
        for b2 in bases:
            for b3 in bases:
                codons.append(b1 + b2 + b3)
    return codons

def main():
    parser = argparse.ArgumentParser(description='Generate codon variants for a specific position')
    parser.add_argument('--fasta', '-f', required=True, help='Input FASTA file')
    parser.add_argument('--codon-table', '-c', required=True, help='Codon table file')
    parser.add_argument('--aa-coord', '-a', type=int, required=True, help='Amino acid coordinate (1-based)')
    parser.add_argument('--seq-id', '-s', required=True, help='Sequence identifier')
    parser.add_argument('--output-fasta', '-o', required=True, help='Output FASTA file')
    parser.add_argument('--output-codon-table', '-v', required=True, help='Output codon file (original codon)')
    parser.add_argument('--left-margin', '-l', type=int, default=2000, help='Left margin')
    parser.add_argument('--right-margin', '-r', type=int, default=1000, help='Right margin')
    parser.add_argument('--gene-start', '-g', type=int, default=None, help='Gene start nucleotide position(1-based) from table')
    parser.add_argument('--gene-end', '-e', type=int, default=None, help='Gene end nucleotide position(1-based) from table')
    
    args = parser.parse_args()
    
    # read codon table
    codon_to_aa = read_codon_table(args.codon_table)
    
    # read input fasta
    header, contig_seq = read_fasta(args.fasta, args.seq_id)

    print(f'contig length: {len(contig_seq)}')
    # making gene seq and checking gene nuc & aalengths
    gene_seq = contig_seq[args.gene_start-1:args.gene_end]
    print(f"sequence length: {len(gene_seq)} nt, {len(gene_seq) // 3} aa")

    # calculate nucleotide position (1-based aa coord to 0-based nucleotide)
    nt_start = (args.aa_coord - 1) * 3
    nt_end = nt_start + 3
    
    # check bounds
    if nt_end > len(gene_seq):
        print(f"error: position {args.aa_coord} is beyond sequence length", file=sys.stderr)
        sys.exit(1)

    # extract and print original codon as sanity check
    original_codon = gene_seq[(nt_start):(nt_end)]
    print(f"original codon: {original_codon}")
    original_aa = codon_to_aa.get(original_codon, 'X')

    print(f"original codon at position {args.aa_coord}: {original_codon} -> {original_aa}")

    # write original codon to output file
    print(f"writing original codon to {args.output_codon_table}")
    with open(args.output_codon_table, 'w') as out:
        out.write(f"coord\tcodon\taa\tleft_margin\tright_margin\n")
        out.write(f"{args.aa_coord}\t{original_codon}\t{original_aa}\t{args.left_margin}\t{args.right_margin}\n")
    
    # get all possible codons
    all_codons = generate_all_codons()

    # # implementing margins
    left_margin = contig_seq[max(0, args.gene_start - 1 - args.left_margin):args.gene_start - 1]
    right_margin = contig_seq[args.gene_end:min(len(contig_seq), args.gene_end + args.right_margin)]

    #### TESTING ####
    # variant_seq = gene_seq[:nt_start] + "XXX" + gene_seq[nt_end:]
    # print(f'here is the variant seq:{variant_seq} with length: {len(variant_seq)}')
    # print(f"gene seq: {gene_seq}")
    # print(f"Original seq: {contig_seq[:args.gene_start-1]}{gene_seq}{contig_seq[args.gene_end:]} and seq length {len(contig_seq[:args.gene_start-1]+gene_seq+contig_seq[args.gene_end:])}")
    # print(f"Margin seq: {left_margin}{gene_seq}{right_margin} and margin length: {len(left_margin+gene_seq+right_margin)}") 

    
    print(f"Generating file: {args.output_fasta}")
    # generate variants
    with open(args.output_fasta, 'w') as out:
        for codon in all_codons:
            aa = codon_to_aa.get(codon, 'X')  # X for unknown

            # replace stop codon symbol (*) with Z
            aa = 'Z' if aa == '*' else aa
            
            # create forward variant
            variant_seq = gene_seq[:nt_start] + codon + gene_seq[nt_end:]
            margin_var_seq = left_margin + variant_seq + right_margin
            # Header without strand suffix for simpler downstream matching
            out.write(f">{args.aa_coord}_{aa}_{codon}\n")
            out.write(f"{margin_var_seq}\n")

if __name__ == "__main__":
    if len(sys.argv) == 1:
        # Test mode (no CLI args passed)
        sys.argv = [
            "test_run",
            "--fasta", "input/human_contigs_src.fasta",
            "--codon-table", "input/codon_table",
            "--aa-coord", "152",
            "--seq-id", "k77_114799",
            "--output-fasta", "output/test_output.fasta",
            "--output-codon-table", "output/test_codons.tsv",
            "--gene-start", "41459",
            "--gene-end", "44200",
            "--left-margin", "2000",
            "--right-margin", "1000"
        ]
    main()