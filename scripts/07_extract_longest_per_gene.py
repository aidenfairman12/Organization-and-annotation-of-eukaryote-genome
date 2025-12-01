#!/usr/bin/env python3
"""
Extract longest sequence per gene from a FASTA file.
By default splits sequence header on the first occurrence of the string '-R' and
uses the left-hand side as the gene id. That handles IDs like:
  GeneX-RA, GeneX-RB, GeneX-RA1, etc.

Usage:
  ./07_extract_longest_per_gene.py input.fasta output.longest.fasta
Options:
  --delimiter STR   delimiter to split header (default: -R)
  --strip-version   strip a trailing ".N" (like ".1") from IDs before matching

This script preserves full headers for the selected (longest) isoform per gene.
"""
import sys
import argparse

def parse_args():
    p = argparse.ArgumentParser(description="Extract longest sequence per gene from FASTA")
    p.add_argument('input', nargs='?', help='input fasta (if omitted, runs defaults for protein and transcript)')
    p.add_argument('output', nargs='?', help='output fasta (one sequence per gene); if omitted when input given, will be <input>.longest.fasta')
    p.add_argument('--delimiter', default='-R', help="delimiter that separates gene id and isoform (default='-R')")
    p.add_argument('--strip-version', action='store_true', help='strip trailing ".N" from id before matching')
    return p.parse_args()


def fasta_iter(fp):
    header = None
    seq_lines = []
    for line in fp:
        line = line.rstrip('\n')
        if not line:
            continue
        if line.startswith('>'):
            if header is not None:
                yield header, ''.join(seq_lines)
            header = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
    if header is not None:
        yield header, ''.join(seq_lines)


def key_from_header(header, delim, strip_version=False):
    # gene key is everything before the first occurrence of delim
    if delim in header:
        key = header.split(delim, 1)[0]
    else:
        # fallback: if delim not present, use whole header up to first whitespace
        key = header.split(None, 1)[0]
    if strip_version:
        # strip trailing .NUMBER (e.g. .1)
        if key.rfind('.') != -1:
            base, dot, tail = key.rpartition('.')
            if tail.isdigit():
                key = base
    return key


def main():
    args = parse_args()
    def process_file(infile, outfile):
        best = {}
        total = 0
        with open(infile) as fh:
            for header, seq in fasta_iter(fh):
                total += 1
                key = key_from_header(header, args.delimiter, args.strip_version)
                cur = best.get(key)
                if cur is None or len(seq) > len(cur[1]):
                    best[key] = (header, seq)
        with open(outfile, 'w') as out:
            for header, seq in best.values():
                out.write('>' + header + '\n')
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + '\n')
        print(f"Read {total} sequences; wrote {len(best)} longest-per-gene sequences to {outfile}")

    # If no input provided, run defaults for proteins and transcripts
    if not args.input:
        prot_in = '/data/users/afairman/euk_org_ann/results/Maker/final/assembly.all.maker.proteins.fasta.renamed.filtered.fasta'
        tran_in = '/data/users/afairman/euk_org_ann/results/Maker/final/assembly.all.maker.transcripts.fasta.renamed.filtered.fasta'
        prot_out = prot_in.replace('.fasta', '.longest.fasta')
        tran_out = tran_in.replace('.fasta', '.longest.fasta')
        process_file(prot_in, prot_out)
        process_file(tran_in, tran_out)
        return

    # If input provided but no output, auto-generate output name
    infile = args.input
    outfile = args.output if args.output else (infile + '.longest.fasta')
    process_file(infile, outfile)

if __name__ == '__main__':
    main()
