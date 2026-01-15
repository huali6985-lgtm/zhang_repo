#!/usr/bin/env bash
# ============================================================
# GWAS Summary Statistics QC Script (GRCh38)
# Genome build: GRCh38
# ============================================================

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input.tsv.gz> <output.QC.tsv.gz>"
    exit 1
fi

IN=$1
OUT=$2

echo "Running GWAS QC (GRCh38)..."
echo "Input : $IN"
echo "Output: $OUT"

zcat "$IN" | \
awk -F'\t' '
BEGIN { OFS = "\t" }
# ----------------------------
# Keep header
# ----------------------------
NR == 1 { print; next }

{
    # ----------------------------
    # Column mapping (1-based)
    # ----------------------------
    chr = $1                  # chromosome
    pos = $2                  # base_pair_location
    ea = $3                   # effect_allele
    oa = $4                   # other_allele
    eaf = $7                  # effect_allele_frequency
    rsid = $9                 # rsID
    vid = $13                 # variant_id
    low_conf = $10            # low_confidence

    # ----------------------------
    # (i) MAF >= 0.01
    # ----------------------------
    maf = (eaf <= 0.5 ? eaf : 1 - eaf)
    if (maf < 0.01) next

    # ----------------------------
    # (ii) Remove low confidence variants
    # ----------------------------
    if (low_conf != "False") next

    # ----------------------------
    # (iii) Biallelic SNPs only
    # ----------------------------
    if (length(ea) != 1 || length(oa) != 1) next
    if (ea !~ /^[ACGT]$/ || oa !~ /^[ACGT]$/) next

    # ----------------------------
    # (iv) Valid rsID
    # ----------------------------
    if (rsid !~ /^rs[0-9]+$/) next

    # ----------------------------
    # (v) Remove MHC region (GRCh38)
    # chr6: 28–34 Mb
    # ----------------------------
    if (chr == 6 && pos >= 28000000 && pos <= 34000000) next

    # ----------------------------
    # (vi) Remove duplicated variants
    # ----------------------------
    if (++seen[vid] > 1) next

    # Passed all QC
    print
}' > "$OUT"

echo "QC completed successfully."
