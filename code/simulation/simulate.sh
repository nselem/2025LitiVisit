# Author Nelly = chatGPT Nov 2025
# nselem84@gmail.com
#
#

READ_LEN=150         # read length
MODEL=HS25           # model ART: HS25=HiSeq 2500, MSv3=MiSeq v3, etc.
CPUS=8

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <genomes.tsv> <total_bp (e.g. 1G or 500000000)> [output_prefix]"
  exit 1
fi

TSV="$1"
TOTAL_BP_RAW="$2"
base=$(basename "$1")
OUTDIR="${base%%.*}_$TOTAL_BP_RAW"
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


# Total pair reads
TOTAL_PAIRS=$(( TOTAL_BP / (2 * READ_LEN) ))
echo "reads $TOTAL_PAIRS basepairs $TOTAL_BP"

mkdir -p $OUTDIR/tmp
TRUTH=$OUTDIR/truth_map.tsv
: > "$TRUTH"  # limpia/crea el archivo de verdad-terreno

while IFS=$'\t' read -r PERCENTAGE FASTA SPECIES; do
  FRAC=$(awk -v p="$PERCENTAGE" 'BEGIN{print p/100}')

  NAME=$(basename "$FASTA" | sed 's/\.fasta$//; s/[^A-Za-z0-9_.-]/_/g')
  # longitud del genome length
  GLEN=$(seqkit fx2tab -l "$FASTA" | awk '{sum+=$NF} END{print sum}')

  # cobertura para lograr FRAC de TOTAL_READS pares (2*READ_LEN bases por par)
  # coverage_i = (FRAC * TOTAL_READS * 2 * READ_LEN) / GLEN
  COV=$(awk -v f="$FRAC" -v N="$TOTAL_PAIRS" -v L="$READ_LEN" -v G="$GLEN" \
        'BEGIN{printf "%.6f", (f*N*2*L)/G}')

  echo "[*] $NAME: genome_len=$GLEN bp, frac=$FRAC, target_cov=$COV"
  art_illumina -ss "$MODEL" -i "$FASTA" -p -l "$READ_LEN" -f "$COV" -m 350 -s 50 \
    -na -qL 26 -qU 40 -o "$OUTDIR/tmp/${NAME}_"

  # Renombrar reads para incluir el nombre del genoma (ground truth)
  zcat -f $OUTDIR/tmp/${NAME}_1.fq | awk -v N="$NAME" 'NR%4==1{printf("@%s|%s\n", N, substr($0,2)); next}1' \
    | gzip -c > "$OUTDIR/tmp/${NAME}_1.fq.gz"
  zcat -f $OUTDIR/tmp/${NAME}_2.fq | awk -v N="$NAME" 'NR%4==1{printf("@%s|%s\n", N, substr($0,2)); next}1' \
    | gzip -c > "$OUTDIR/tmp/${NAME}_2.fq.gz"

  # Mapea IDs -> origen
  zcat "$OUTDIR/tmp/${NAME}_1.fq.gz" | awk 'NR%4==1{split(substr($0,2),a,"|"); print a[2]"\t"a[1]"\t1"}' >> "$TRUTH"
  zcat "$OUTDIR/tmp/${NAME}_2.fq.gz" | awk 'NR%4==1{split(substr($0,2),a,"|"); print a[2]"\t"a[1]"\t2"}' >> "$TRUTH"

done < $TSV

# Mezcla final
cat $OUTDIR/tmp/*_1.fq.gz > $OUTDIR/$OUTDIR\_short_1.fastq.gz
cat $OUTDIR/tmp/*_2.fq.gz > $OUTDIR/$OUTDIR\_short_2.fastq.gz

echo "Listo:"
echo " - Lecturas: $OUTDIR/metagenome_R1.fastq.gz, $OUTDIR/metagenome_R2.fastq.gz"
echo " - Ground truth (readID -> genoma): $TRUTH"
