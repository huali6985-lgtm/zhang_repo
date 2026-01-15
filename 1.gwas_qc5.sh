#!/usr/bin/env bash
# ============================================================
# GWAS Summary Statistics QC Script (GRCh38)
# For liftOver summary stats (no INFO, no variant_id)
# ============================================================

set -euo pipefail

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input.txt> <output_qc.txt>"
    exit 1
fi

IN=$1
OUT=$2

echo "Running GWAS QC (GRCh38, liftOver data)..."
echo "Input : $IN"
echo "Output: $OUT"

awk -F'\t' '
BEGIN { OFS = "\t" }

# ----------------------------
# Keep header
# ----------------------------
NR == 1 {
    print
    next
}

{
    # ----------------------------
    # Column mapping (IgA_liftOver)
    # ----------------------------
    rsid = $1        # SNP
    chr  = $2        # CHR_hg38
    pos  = $3        # BP_hg38
    ea   = $4        # effect_allele
    oa   = $5        # other_allele
    eaf  = $6        # effect allele frequency

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
    # (v) Remove duplicated variants
    # use chr:pos:alleles as ID
    # ----------------------------
    vid = chr ":" pos ":" ea ":" oa
    if (++seen[vid] > 1) next

    # ----------------------------
    # Passed all QC
    # ----------------------------
    print
}
' "$IN" > "$OUT"

echo "QC completed successfully."
