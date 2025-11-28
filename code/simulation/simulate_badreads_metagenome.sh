#!/usr/bin/env bash
# author chatGPT and nselem84@gmail.com
# bash simulate.sh myfile.tsv <quantity in bp>
#
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <genomes.tsv> <total_bp (e.g. 1G or 500000000)> [output_prefix]"
  exit 1
fi

TSV="$1"
TOTAL_BP_RAW="$2"
base=$(basename "$1")
OUTDIR="${base%%.*}_$TOTAL_BP_RAW"
OUTPREFIX="${3:-metagenome}"
mkdir $OUTDIR

# --- convert total_bp to a numeric value (supports K/M/G suffixes) ---
if [[ "$TOTAL_BP_RAW" =~ ^[0-9]+$ ]]; then
  TOTAL_BP="$TOTAL_BP_RAW"
elif [[ "$TOTAL_BP_RAW" =~ ^([0-9]+)([KMG])$ ]]; then
  base="${BASH_REMATCH[1]}"
  suf="${BASH_REMATCH[2]}"
  case "$suf" in
    K) TOTAL_BP=$((base * 1000));;
    M) TOTAL_BP=$((base * 1000000));;
    G) TOTAL_BP=$((base * 1000000000));;
    *) echo "Unknown suffix in total_bp: $suf"; exit 1;;
  esac
else
  echo "ERROR: total_bp must be an integer (bp) or end in K/M/G, e.g. 500000000 or 1G"
  exit 1
fi

echo "# Target total bp: $TOTAL_BP"

# --- 1. Compute the sum of weights (only rows with weight > 0) ---
total_weight=0

while IFS=$'\t' read -r weight fasta desc; do
  [[ -z "${weight:-}" ]] && continue
  [[ "$weight" =~ ^# ]] && continue
  if ! [[ "$weight" =~ ^[0-9]+$ ]]; then
    echo "Warning: non-numeric weight encountered, skipping: $weight" >&2
    continue
  fi
  if (( weight > 0 )); then
    total_weight=$((total_weight + weight))
  fi
done < "$TSV"

if (( total_weight == 0 )); then
  echo "ERROR: total weight is 0 (are all proportions 0?)."
  exit 1
fi

echo "# Sum of weights (>0 only): $total_weight"

# --- 2. For each genome, compute bp and run badreads ---
tmp_list=""

while IFS=$'\t' read -r weight fasta desc; do
  [[ -z "${weight:-}" ]] && continue
  [[ "$weight" =~ ^# ]] && continue
  if ! [[ "$weight" =~ ^[0-9]+$ ]]; then
    continue
  fi

  if (( weight == 0 )); then
    echo "# Skipping $fasta (weight = 0)"
    continue
  fi

  if [[ ! -f "$fasta" ]]; then
    echo "ERROR: FASTA file not found: $fasta" >&2
    exit 1
  fi

  # bp corresponding to this genome
  quantity_bp=$(( TOTAL_BP * weight / total_weight ))

  if (( quantity_bp == 0 )); then
    echo "Warning: $fasta receives 0 bp after normalization, skipping."
    continue
  fi

  base=$(basename "$fasta")
  base_noext="${base%%.*}"
  out_fastq="${OUTPREFIX}.${base_noext}.fastq"

  echo "# Simulating for $fasta"
  echo "#   Weight: $weight  ->  ${quantity_bp} bp"
#  badread simulate  --reference /home/shared/referencedgenomes/GCF_010111755.1_ASM1011175v1_genomic.fna --quantity 25M    --length 10000,2000   --error_model nanopore2023   --qscore_model nanopore2023   > cgla_nanopore.fastq

  badread simulate \
    --reference "$fasta" \
    --quantity "${quantity_bp}" \
    --error_model nanopore2023 \
    --qscore_model nanopore2023 \
    --length 10000,2000 > "$OUTDIR/$out_fastq"

  tmp_list="$tmp_list $OUTDIR/$out_fastq"

done < "$TSV"

# --- 3. Combine FASTQ files into one ---
combined_fastq="$OUTDIR/${OUTPREFIX}.combined.fastq"

if [[ -z "$tmp_list" ]]; then
  echo "ERROR: no FASTQ files were generated (were all weights 0?)."
  exit 1
fi

echo "# Combining files into: $combined_fastq"
cat $tmp_list > "$combined_fastq"

# --- 4. If seqkit is installed, shuffle the reads ---
if command -v seqkit >/dev/null 2>&1; then
  shuffled_fastq="$OUTDIR/${OUTPREFIX}.combined.shuffled.fastq"
  echo "# seqkit found, shuffling reads into: $shuffled_fastq"
  seqkit shuffle "$combined_fastq" > "$shuffled_fastq"
  echo "# Final metagenome (shuffled): $shuffled_fastq"
else
  echo "# seqkit NOT found; leaving unshuffled metagenome: $combined_fastq"
fi

gzip $shuffled_fastq
cp $TSV $OUTDIR/.

echo "# Done."

