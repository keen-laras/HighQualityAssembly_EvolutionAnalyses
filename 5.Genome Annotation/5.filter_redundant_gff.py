#!/usr/bin/env python3

import sys
from collections import defaultdict

if len(sys.argv) != 4:
    sys.exit(
        "Usage: python filter_redundant_gff.py keep_ids.txt input.gff output.gff"
    )

keep_ids_file, gff_file, out_file = sys.argv[1:]

# --------------------------------------------------
# Load non-redundant transcript / mRNA IDs
# --------------------------------------------------
keep_tx = set()
with open(keep_ids_file) as f:
    for line in f:
        tid = line.strip()
        if tid:
            keep_tx.add(tid)

# --------------------------------------------------
# Containers
# --------------------------------------------------
kept_genes = set()
kept_tx_parents = {}   # transcript/mRNA -> gene
gene_lines = {}
tx_lines = {}
child_lines = defaultdict(list)

# --------------------------------------------------
# First pass: parse relationships
# --------------------------------------------------
with open(gff_file) as f:
    for line in f:
        if line.startswith("#"):
            continue

        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue

        feature = cols[2]
        attrs = cols[8]

        def get_attr(key):
            for x in attrs.split(";"):
                if x.startswith(key + "="):
                    return x.replace(key + "=", "")
            return None

        # gene
        if feature == "gene":
            gid = get_attr("ID")
            if gid:
                gene_lines[gid] = line

        # transcript or mRNA
        elif feature in ("transcript", "mRNA"):
            tid = get_attr("ID")
            gid = get_attr("Parent")
            if tid:
                tx_lines[tid] = line
            if tid in keep_tx and gid:
                kept_tx_parents[tid] = gid
                kept_genes.add(gid)

        # child features
        elif "Parent=" in attrs:
            parent = get_attr("Parent")
            if parent:
                child_lines[parent].append(line)

# --------------------------------------------------
# Second pass: write clean GFF
# --------------------------------------------------
with open(out_file, "w") as out:
    # genes
    for gid in kept_genes:
        if gid in gene_lines:
            out.write(gene_lines[gid])

    # transcripts / mRNA
    for tid in keep_tx:
        if tid in tx_lines:
            out.write(tx_lines[tid])

    # child features
    for tid in keep_tx:
        for cl in child_lines.get(tid, []):
            out.write(cl)
