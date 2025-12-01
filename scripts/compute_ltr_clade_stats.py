#!/usr/bin/env python3
import csv
from collections import defaultdict
import pandas as pd

# Paths (adjust if needed)
gff_path = '/data/users/afairman/euk_org_ann/results/edta_output/assembly.fasta.mod.EDTA.raw/assembly.fasta.mod.LTR.intact.raw.gff3'
cls_path = '/data/users/afairman/euk_org_ann/results/edta_output/assembly.fasta.mod.EDTA.raw/assembly.fasta.mod.LTR.intact.raw.fa.rexdb-plant.cls.tsv'
out_csv = '/data/users/afairman/euk_org_ann/results/edta_output/edta_ltr_clade_identity_summary.csv'
out_bins_csv = '/data/users/afairman/euk_org_ann/results/edta_output/edta_ltr_clade_identity_bins.csv'

# Parse GFF: capture Name -> ltr_identity (as percent)
gff_map = {}
with open(gff_path) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        cols = line.split('\t')
        if len(cols) < 9:
            continue
        feature = cols[2]
        if feature != 'repeat_region':
            continue
        attrs = cols[8]
        ad = {}
        for part in attrs.split(';'):
            if '=' in part:
                k,v = part.split('=',1)
                ad[k]=v
        name = ad.get('Name')
        li = ad.get('ltr_identity')
        if name and li:
            try:
                # Some values are '1' or '1.0000' etc
                pct = float(li) * 100.0
            except Exception:
                try:
                    pct = float(li)
                except Exception:
                    continue
            gff_map[name] = pct

print(f'Parsed {len(gff_map)} repeat_region entries from GFF')

# Parse classifier TSV and join
rows = []
missing = 0
with open(cls_path) as fh:
    for line in fh:
        line=line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) < 4:
            continue
        te = parts[0]
        clade = parts[3]
        complete = parts[4] if len(parts) > 4 else ''
        # key before any '#'
        key = te.split('#')[0]
        if key in gff_map:
            pct = gff_map[key]
            rows.append({'element':key,'clade':clade,'percent_identity':pct,'complete':complete})
        else:
            missing += 1

print(f'Found {len(rows)} classified elements with ltr_identity; {missing} classified elements missing in GFF map')

if not rows:
    raise SystemExit('No joined rows; check paths and identifiers')

# DataFrame and aggregations
df = pd.DataFrame(rows)
# Total per clade
summary = df.groupby('clade').agg(n_elements=('element','count'),
                                   mean_identity=('percent_identity','mean'),
                                   median_identity=('percent_identity','median')).reset_index()

# Bins: 80-90,90-95,95-99,99-100
bins = [(80,90),(90,95),(95,99),(99,100)]
for lo,hi in bins:
    col = f'{lo}-{hi}'
    summary[col] = 0

# Count per bin per clade
bin_counts = defaultdict(lambda: defaultdict(int))
for _,r in df.iterrows():
    cl = r['clade']
    pct = r['percent_identity']
    placed = False
    for lo,hi in bins:
        if pct >= lo and (pct < hi or (hi==100 and pct<=hi)):
            bin_counts[cl][f'{lo}-{hi}'] += 1
            placed = True
            break
    if not placed:
        # outside bins, ignore or could record
        pass

# Fill summary with bin counts
for idx,row in summary.iterrows():
    cl = row['clade']
    for lo,hi in bins:
        summary.at[idx,f'{lo}-{hi}'] = bin_counts[cl].get(f'{lo}-{hi}',0)

# Save outputs
summary.to_csv(out_csv,index=False)
# Also save full joined table and per-clade binned table
bin_df = summary.sort_values('n_elements',ascending=False)
bin_df.to_csv(out_bins_csv,index=False)

print('\nPer-clade summary written to:')
print(out_csv)
print(out_bins_csv)

# Print a compact table for quick review
print('\nTop clades by element count:')
print(summary.sort_values('n_elements',ascending=False).head(20).to_string(index=False))

# Also write full join table
full_out = '/data/users/afairman/euk_org_ann/results/edta_output/edta_ltr_clade_identity_joined.tsv'
df.to_csv(full_out,sep='\t',index=False)
print('\nFull joined table written to:',full_out)
