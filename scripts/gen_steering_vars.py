#!/usr/bin/env python3
import sys
import os
import csv
import argparse
from typing import Optional
from Bio import SeqIO

def ensure_dir(path: str) -> None:
    parent = os.path.dirname(path)
    if parent and not os.path.exists(parent):
        os.makedirs(parent, exist_ok=True)


def process_one_entry(
    contig_seq: str,
    aid: str,
    seq_id: str,
    gene_start: int,
    gene_end: int,
    aa_coord: int,
    tgt_codon: str,
    left_margin_nts: int,
    right_margin_nts: int,
    global_fasta_path: Optional[str] = None,
    query_table_path: Optional[str] = None,
) -> None:
    gene_seq = contig_seq[gene_start - 1 : gene_end]
    nt_start = (aa_coord - 1) * 3
    nt_end = nt_start + 3

    if nt_end > len(gene_seq):
        print(f"error: position {aa_coord} is beyond sequence length", file=sys.stderr)
        sys.exit(1)

    original_codon = gene_seq[nt_start:nt_end]

    left_margin = contig_seq[max(0, gene_start - 1 - left_margin_nts) : gene_start - 1]
    right_margin = contig_seq[gene_end : min(len(contig_seq), gene_end + right_margin_nts)]

    aid_lower = aid.lower()
    # Build the two variant sequences once
    codons = [original_codon, tgt_codon]

    # Optionally also append to a global FASTA for the whole batch
    if global_fasta_path:
        ensure_dir(global_fasta_path)
        with open(global_fasta_path, "a") as gfo:
            for codon in codons:
                variant_seq = gene_seq[:nt_start] + codon + gene_seq[nt_end:]
                margin_var_seq = left_margin + variant_seq + right_margin
                gfo.write(f">{aid_lower}_{codon}\n")
                gfo.write(f"{margin_var_seq}\n")

    # Optionally append to a global query table mapping seq_id to [start,end]
    # We want 201 nt starting at the mutation nucleotide position (2nd base of codon)
    if query_table_path:
        # compute 1-based position of mutation nucleotide within margin_var_seq
        # second base of the codon: nt_start (0-based codon start) + 1
        mut_nt_pos_1b_in_variant = len(left_margin) + (nt_start + 1) + 1
        # For each codon sequence, the start/end are identical by construction
        start_pos = mut_nt_pos_1b_in_variant
        # length of any margin_var_seq with either codon is the same
        total_len = len(left_margin) + len(gene_seq) + len(right_margin)
        end_pos = min(start_pos + 200, total_len)
        with open(query_table_path, "a") as qf:
            for codon in codons:
                qf.write(f"{aid_lower}_{codon}\t{start_pos}\t{end_pos}\n")


def main():
    parser = argparse.ArgumentParser(description="Generate codon variants (batch mode from TSV)")
    parser.add_argument("--fasta", "-f", required=True, help="Input FASTA file")

    # Batch mode inputs
    parser.add_argument(
        "--updated-data",
        "-u",
        required=True,
        help="TSV with columns: aid, gene, contig, start, end, strand, flipped, src_codon, tgt_codon, mut_pos",
    )
    parser.add_argument("--left-margin", "-l", type=int, default=2000, help="Left margin")
    parser.add_argument("--right-margin", "-r", type=int, default=1000, help="Right margin")

    # Output FASTA path (used in both modes). In batch mode, sequences are appended.
    parser.add_argument("--output-fasta", "-o", help="If set: write sequences to this FASTA (append in batch mode)")
    parser.add_argument("--query-table", "-q", help="If set: write seq_id/start/end here (append in batch mode)")

    args = parser.parse_args()

    # Read FASTA (can contain multiple sequences; we index by id when needed)
    fasta_records = {rec.id: str(rec.seq) for rec in SeqIO.parse(args.fasta, "fasta")}

    # Batch mode: iterate updated_data.tsv
    tsv_path = args.updated_data
    if not os.path.isfile(tsv_path):
        print(f"error: updated data TSV not found: {tsv_path}", file=sys.stderr)
        sys.exit(1)
    # If a global output FASTA/TABLE is set, truncate and init header so we start fresh
    if args.output_fasta:
        ensure_dir(args.output_fasta)
        open(args.output_fasta, "w").close()
    if args.query_table:
        ensure_dir(args.query_table)
        with open(args.query_table, "w") as qf:
            qf.write("seq_id\tstart\tend\n")
    with open(tsv_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            aid = row.get("aid")
            seq_id = row.get("contig")
            start = int(row.get("start"))
            end = int(row.get("end"))
            tgt_codon = row.get("tgt_codon")
            aa_coord = int(row.get("mut_pos"))

            if seq_id not in fasta_records:
                print(f"warning: seq id {seq_id} not found in FASTA; skipping", file=sys.stderr)
                continue
            contig_seq = fasta_records[seq_id]

            process_one_entry(
                contig_seq=contig_seq,
                aid=aid,
                seq_id=seq_id,
                gene_start=start,
                gene_end=end,
                aa_coord=aa_coord,
                tgt_codon=tgt_codon,
                left_margin_nts=args.left_margin,
                right_margin_nts=args.right_margin,
                global_fasta_path=(os.path.abspath(args.output_fasta) if args.output_fasta else None),
                query_table_path=(os.path.abspath(args.query_table) if args.query_table else None),
            )
    return

if __name__ == "__main__":
    main()