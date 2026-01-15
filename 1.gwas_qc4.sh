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
NR==1 { print; next }

{
    # ----------------------------
    # Column mapping (1-based)
    # ----------------------------
    chr = $3                   # hm_chrom
    pos = $4                   # hm_pos
    ea = $6                    # hm_effect_allele
    oa = $5                    # hm_other_allele
    eaf = $11                  # effect_allele_frequency
    rsid = $2                  # hm_rsid
    vid = $1                   # hm_variant_id 

    # ----------------------------
    # (i) MAF >= 0.01
    # ----------------------------
    maf = (eaf <= 0.5 ? eaf : 1 - eaf)
    if (maf < 0.01) next

    # ----------------------------
    # (ii) Biallelic SNPs only
    # ----------------------------
    if (length(ea) != 1 || length(oa) != 1) next
    if (ea !~ /^[ACGT]$/ || oa !~ /^[ACGT]$/) next

    # ----------------------------
    # (iii) Valid rsID
    # ----------------------------
    if (rsid !~ /^rs[0-9]+$/) next

    # ----------------------------
    # (iv) Remove MHC region (GRCh38)
    # chr6: 28–34 Mb
    # ----------------------------
    if (chr == 6 && pos >= 28000000 && pos <= 34000000) next

    # ----------------------------
    # (v) Remove duplicated variants by hm_variant_id
    # ----------------------------
    if (++seen[vid] > 1) next

    # Passed all QC
    print
}' | gzip > "$OUT"

echo "QC completed successfully."
