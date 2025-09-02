#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# WES pipeline — thesis-style (linear, single file)
# GRCh38 + dbSNP v155 + GATK 4.6 + VEP (+ CADD v1.7, REVEL, AlphaMissense) + optional gnomAD v4 + optional VASE
# Edit the variables below, then run.
# ============================================================

# -----------------------------
# USER CONFIG (edit these)
# -----------------------------

# Mode A: COHORT — set SAMPLES_TSV to a TSV with header: sample<TAB>r1<TAB>r2
SAMPLES_TSV=""

# Mode B: SINGLE SAMPLE — set SAMPLE_ID and FASTQs; leave SAMPLES_TSV empty
SAMPLE_ID="SAMPLE01"
R1="/path/to/${SAMPLE_ID}_R1.fastq.gz"
R2="/path/to/${SAMPLE_ID}_R2.fastq.gz"

# References / known sites (GRCh38)
REF="/path/to/Homo_sapiens_assembly38.fasta"
DBSNP="/path/to/dbsnp_155.hg38.vcf.gz"
MILLS="/path/to/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_INDELS="/path/to/Homo_sapiens_assembly38.known_indels.vcf.gz"
PHASE3="/path/to/1000G_phase3_v4_20130502.sites.hg38.vcf.gz"
BED="/path/to/exome_targets.bed"   # optional; leave empty to call genome-wide

# VEP + plugins (versions used in Feb–Sep 2024)
VEP_CACHE="/path/to/.vep"   # e.g. v110/GRCh38 cache
CADD_SNV="/path/to/CADD_v1.7/whole_genome_SNVs.tsv.gz"
CADD_INDEL="/path/to/CADD_v1.7/gnomad.genomes.r4.0.indel.tsv.gz"
REVEL_TSV="/path/to/revel_grch38.tsv.gz"
ALPHAMISS_TSV="/path/to/AlphaMissense_hg38.tsv.gz"
UTR_FILE=""                 # optional (uORF_5UTR_GRCh38_PUBLIC.txt)

# Optional frequency overlay
GNOMAD_V4_SITES=""          # optional: /path/to/gnomad.v4.0.sites.vcf.bgz

# Optional VASE resources
VASE_BIN=""                 # optional: /path/to/vase
CLINVAR=""                  # optional: /path/to/clinvar_202407.vcf.gz
PED=""                      # optional: /path/to/family.ped

# Runtime and outputs
THREADS=12
JAVA_OPTS="-Xmx20g"
OUTDIR="results"
LOGDIR="logs"

# -----------------------------
# PREP
# -----------------------------
mkdir -p "$OUTDIR"/{align,gvcf,vcf,annot} "$LOGDIR"

# Quick tool presence checks (minimal)
for t in trim_galore bwa samtools picard gatk bcftools vep; do
  command -v "$t" >/dev/null || { echo "[ERR] Missing tool: $t"; exit 1; }
done

# -----------------------------
# STAGE 1 — VARIANT DISCOVERY
# -----------------------------
echo "=== Stage 1: Variant Discovery (FASTQ → VCF) ==="

if [[ -n "$SAMPLES_TSV" ]]; then
  # ---- COHORT MODE ----
  [[ -f "$SAMPLES_TSV" ]] || { echo "[ERR] samples.tsv not found"; exit 1; }

  while IFS=$'\t' read -r sid r1 r2; do
    [[ "$sid" == "sample" ]] && continue
    echo "--- Sample: $sid ---"

    # Step 1: Trim Galore
    trim_galore -q 20 --illumina --paired --gzip --length 20 \
      -o "$OUTDIR/align" "$r1" "$r2"

    # Step 2: BWA-MEM → BAM
    bwa mem -t "$THREADS" -M "$REF" \
      "$OUTDIR/align/${sid}_R1_001_val_1.fq.gz" \
      "$OUTDIR/align/${sid}_R2_001_val_2.fq.gz" \
      -R "@RG\tID:${sid}\tSM:${sid}\tPL:ILLUMINA\tLB:${sid}_exome" \
    | samtools view -Sb - > "$OUTDIR/align/${sid}.bam"

    # Step 3: Sort + MarkDuplicates
    picard SortSam -I "$OUTDIR/align/${sid}.bam" \
      -O "$OUTDIR/align/${sid}.sort.bam" -SO coordinate -CREATE_INDEX true
    picard MarkDuplicates -I "$OUTDIR/align/${sid}.sort.bam" \
      -O "$OUTDIR/align/${sid}.dedup.bam" \
      -M "$OUTDIR/align/${sid}.metrics.txt" -CREATE_INDEX true

    # Step 4: BQSR (BaseRecalibrator + ApplyBQSR)
    gatk BaseRecalibrator -R "$REF" -I "$OUTDIR/align/${sid}.dedup.bam" \
      --known-sites "$PHASE3" --known-sites "$DBSNP" \
      --known-sites "$KNOWN_INDELS" --known-sites "$MILLS" \
      -O "$OUTDIR/align/${sid}.recal.table"

    gatk ApplyBQSR -R "$REF" -I "$OUTDIR/align/${sid}.dedup.bam" \
      --bqsr-recal-file "$OUTDIR/align/${sid}.recal.table" \
      -O "$OUTDIR/align/${sid}.bqsr.bam"

    # Step 5: HaplotypeCaller (GVCF)
    EXTRA=""; [[ -n "$BED" ]] && EXTRA="-L ${BED}"
    gatk --java-options "-Xmx4g" HaplotypeCaller -R "$REF" \
      -I "$OUTDIR/align/${sid}.bqsr.bam" \
      -O "$OUTDIR/gvcf/${sid}.g.vcf.gz" -ERC GVCF \
      -D "$DBSNP" $EXTRA

    tabix -p vcf "$OUTDIR/gvcf/${sid}.g.vcf.gz" || true
  done < "$SAMPLES_TSV"

  # Step 6: CombineGVCFs → GenotypeGVCFs (joint)
  gatk CombineGVCFs -R "$REF" \
    $(for g in "$OUTDIR"/gvcf/*.g.vcf.gz; do echo "-V $g"; done) \
    -O "$OUTDIR/vcf/cohort.g.vcf.gz"

  gatk GenotypeGVCFs -R "$REF" \
    -V "$OUTDIR/vcf/cohort.g.vcf.gz" \
    -O "$OUTDIR/vcf/cohort.raw.vcf.gz"

else
  # ---- SINGLE-SAMPLE MODE ----
  sid="$SAMPLE_ID"
  echo "--- Sample: $sid ---"

  # Step 1: Trim Galore
  trim_galore -q 20 --illumina --paired --gzip --length 20 \
    -o "$OUTDIR/align" "$R1" "$R2"

  # Step 2: BWA-MEM → BAM
  bwa mem -t "$THREADS" -M "$REF" \
    "$OUTDIR/align/${sid}_R1_001_val_1.fq.gz" \
    "$OUTDIR/align/${sid}_R2_001_val_2.fq.gz" \
    -R "@RG\tID:${sid}\tSM:${sid}\tPL:ILLUMINA\tLB:${sid}_exome" \
  | samtools view -Sb - > "$OUTDIR/align/${sid}.bam"

  # Step 3: Sort + MarkDuplicates
  picard SortSam -I "$OUTDIR/align/${sid}.bam" \
    -O "$OUTDIR/align/${sid}.sort.bam" -SO coordinate -CREATE_INDEX true
  picard MarkDuplicates -I "$OUTDIR/align/${sid}.sort.bam" \
    -O "$OUTDIR/align/${sid}.dedup.bam" \
    -M "$OUTDIR/align/${sid}.metrics.txt" -CREATE_INDEX true

  # Step 4: BQSR
  gatk BaseRecalibrator -R "$REF" -I "$OUTDIR/align/${sid}.dedup.bam" \
    --known-sites "$PHASE3" --known-sites "$DBSNP" \
    --known-sites "$KNOWN_INDELS" --known-sites "$MILLS" \
    -O "$OUTDIR/align/${sid}.recal.table"

  gatk ApplyBQSR -R "$REF" -I "$OUTDIR/align/${sid}.dedup.bam" \
    --bqsr-recal-file "$OUTDIR/align/${sid}.recal.table" \
    -O "$OUTDIR/align/${sid}.bqsr.bam"

  # Step 5: HaplotypeCaller (GVCF)
  EXTRA=""; [[ -n "$BED" ]] && EXTRA="-L ${BED}"
  gatk --java-options "-Xmx4g" HaplotypeCaller -R "$REF" \
    -I "$OUTDIR/align/${sid}.bqsr.bam" \
    -O "$OUTDIR/gvcf/${sid}.g.vcf.gz" -ERC GVCF \
    -D "$DBSNP" $EXTRA

  tabix -p vcf "$OUTDIR/gvcf/${sid}.g.vcf.gz" || true

  # Step 6: GenotypeGVCFs (single sample)
  gatk GenotypeGVCFs -R "$REF" \
    -V "$OUTDIR/gvcf/${sid}.g.vcf.gz" \
    -O "$OUTDIR/vcf/cohort.raw.vcf.gz"
fi

# Step 7: Normalize multi-allelics (left-align + split) BEFORE validation
bcftools norm -f "$REF" -m -both \
  "$OUTDIR/vcf/cohort.raw.vcf.gz" \
  -Oz -o "$OUTDIR/vcf/cohort.norm.vcf.gz"
tabix -f -p vcf "$OUTDIR/vcf/cohort.norm.vcf.gz" || true

# Step 8: Validate normalized VCF (relaxed: exclude ALL checks)
gatk ValidateVariants -R "$REF" \
  -V "$OUTDIR/vcf/cohort.norm.vcf.gz" \
  --validation-type-to-exclude ALL \
  || true

# Step 9: Split → filter → merge (SNPs / INDELs) — thresholds per thesis
gatk SelectVariants -V "$OUTDIR/vcf/cohort.norm.vcf.gz" -select-type SNP   -O "$OUTDIR/vcf/snps.vcf.gz"
gatk SelectVariants -V "$OUTDIR/vcf/cohort.norm.vcf.gz" -select-type INDEL -O "$OUTDIR/vcf/indels.vcf.gz"

# SNP filters
gatk VariantFiltration -V "$OUTDIR/vcf/snps.vcf.gz" \
  --filter-name "QD2"    --filter-expression "QD < 2.0" \
  --filter-name "QUAL35" --filter-expression "QUAL < 35.0" \
  --filter-name "SOR4"   --filter-expression "SOR > 4.0" \
  --filter-name "FS55"   --filter-expression "FS > 55.0" \
  --filter-name "MQ40"   --filter-expression "MQ < 40.0" \
  --filter-name "MQRankSum-11.5"     --filter-expression "MQRankSum < -11.5" \
  --filter-name "ReadPosRankSum-7.5" --filter-expression "ReadPosRankSum < -7.5" \
  -O "$OUTDIR/vcf/snps.filtered.vcf.gz"

# INDEL filters
gatk VariantFiltration -V "$OUTDIR/vcf/indels.vcf.gz" \
  --filter-name "QD2"    --filter-expression "QD < 2.0" \
  --filter-name "QUAL30" --filter-expression "QUAL < 30.0" \
  --filter-name "FS200"  --filter-expression "FS > 200.0" \
  --filter-name "ReadPosRankSum-20" --filter-expression "ReadPosRankSum < -20.0" \
  -O "$OUTDIR/vcf/indels.filtered.vcf.gz"

gatk MergeVcfs \
  -I "$OUTDIR/vcf/snps.filtered.vcf.gz" \
  -I "$OUTDIR/vcf/indels.filtered.vcf.gz" \
  -O "$OUTDIR/vcf/cohort.filtered.vcf.gz"

echo "=== Stage 1 complete ==="
echo "   → $OUTDIR/vcf/cohort.norm.vcf.gz"
echo "   → $OUTDIR/vcf/cohort.filtered.vcf.gz"

# -----------------------------
# STAGE 2 — INTERPRETATION
# -----------------------------
echo "=== Stage 2: Interpretation (VEP → optional gnomAD v4 → optional VASE) ==="

# Step 10–12: VEP (+ plugins) → optional gnomAD v4 overlay
VEP_PLUGINS="--plugin CADD,snv=${CADD_SNV},indels=${CADD_INDEL} \
             --plugin REVEL,file=${REVEL_TSV},no_match=1 \
             --plugin AlphaMissense,file=${ALPHAMISS_TSV},cols=all"
[[ -n "$UTR_FILE" ]] && VEP_PLUGINS="${VEP_PLUGINS} --plugin UTRAnnotator,file=${UTR_FILE}"

vep --offline --vcf --assembly GRCh38 \
    --dir_cache "$VEP_CACHE" \
    -i "$OUTDIR/vcf/cohort.filtered.vcf.gz" \
    -o "$OUTDIR/annot/cohort.vep.vcf" \
    --minimal --everything --sift b --polyphen b --canonical \
    --hgvs --pubmed --fasta "$REF" \
    --af_gnomadg \
    $VEP_PLUGINS

if [[ -n "$GNOMAD_V4_SITES" ]]; then
  bcftools annotate \
    -a "$GNOMAD_V4_SITES" \
    -c INFO/AF,INFO/AN,INFO/AC \
    -Oz -o "$OUTDIR/annot/cohort.vep.gnomad4.vcf.gz" \
    "$OUTDIR/annot/cohort.vep.vcf"
  tabix -p vcf "$OUTDIR/annot/cohort.vep.gnomad4.vcf.gz" || true
else
  bgzip -c "$OUTDIR/annot/cohort.vep.vcf" > "$OUTDIR/annot/cohort.vep.gnomad4.vcf.gz"
  tabix -p vcf "$OUTDIR/annot/cohort.vep.gnomad4.vcf.gz" || true
fi

# Step 13–14: optional VASE (dominant & recessive) + reporter
if [[ -n "$VASE_BIN" ]]; then
  VCF_IN="$OUTDIR/annot/cohort.vep.gnomad4.vcf.gz"

  "$VASE_BIN" -i "$VCF_IN" -o "$OUTDIR/annot/vase_dom.vcf" \
    ${GNOMAD_V4_SITES:+--gnomad "$GNOMAD_V4_SITES"} \
    --dbsnp "$DBSNP" --build 155 \
    ${CLINVAR:+--clinvar "$CLINVAR"} \
    --csq --canonical --pass_filters --dp 5 --gq 20 \
    --vep_af gnomAD_AF --freq 0.3 \
    --cadd_directory "$(dirname "$CADD_SNV")" --cadd_phred 15.0 \
    ${PED:+--singleton_dominant vcf_sample_ids.txt -ped "$PED"} \
    --debug --clinvar_path \
    --min_families 1 \
    --csq intron_variant missense_variant splice_region_variant stop_gained frameshift_variant splice_acceptor_variant inframe_insertion inframe_deletion || true

  "$VASE_BIN" -i "$VCF_IN" -o "$OUTDIR/annot/vase_rec.vcf" \
    ${GNOMAD_V4_SITES:+--gnomad "$GNOMAD_V4_SITES"} \
    --dbsnp "$DBSNP" --build 155 \
    ${CLINVAR:+--clinvar "$CLINVAR"} \
    --csq --canonical --pass_filters --dp 5 --gq 20 \
    --vep_af gnomAD_AF --freq 0.3 \
    --cadd_directory "$(dirname "$CADD_SNV")" --cadd_phred 15.0 \
    ${PED:+--singleton_recessive vcf_sample_ids.txt -ped "$PED"} \
    --debug --clinvar_path \
    --min_families 1 \
    --csq intron_variant missense_variant splice_region_variant stop_gained frameshift_variant splice_acceptor_variant inframe_insertion inframe_deletion || true

  if command -v vase_reporter >/dev/null; then
    picard MergeVcfs I="$OUTDIR/annot/vase_dom.vcf" I="$OUTDIR/annot/vase_rec.vcf" O="$OUTDIR/annot/AI_final_variants.vcf" || true
    vase_reporter "$OUTDIR/annot/AI_final_variants.vcf" "$OUTDIR/annot/AI_final_variants.xlsx" \
      ${PED:+--ped "$PED"} --force --quiet --debug || true
  fi
fi

echo "=== Pipeline complete ==="
echo "Outputs:"
echo "  Stage 1 → $OUTDIR/vcf/cohort.norm.vcf.gz ; $OUTDIR/vcf/cohort.filtered.vcf.gz"
echo "  Stage 2 → $OUTDIR/annot/cohort.vep.gnomad4.vcf.gz (plus VASE outputs if run)"
