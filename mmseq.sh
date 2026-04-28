#!/bin/sh

# Check if the correct number of arguments is provided
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <input_fasta> <output_dir> <tsv_name> <db_path>"
    exit 1
fi

INPUT=$1
OUTDIR=$2
TSV_NAME=$3
DB=$4  # Now passed as an argument instead of hardcoded

mmseqs createdb "$INPUT" "${OUTDIR}sample.contigs"

mmseqs taxonomy "${OUTDIR}sample.contigs" "$DB" "${OUTDIR}sample.assignments" "${OUTDIR}sample.tmpFolder" \
    --tax-lineage 1 --majority 0.7 --vote-mode 1 --lca-mode 3 --orf-filter 1

mmseqs createtsv "${OUTDIR}sample.contigs" "${OUTDIR}sample.assignments" "${OUTDIR}${TSV_NAME}"

# Optional: Clean up
# rm "${OUTDIR}sample.contigs"*
# rm "${OUTDIR}sample.assignments"*
# rm -r "${OUTDIR}sample.tmpFolder"

