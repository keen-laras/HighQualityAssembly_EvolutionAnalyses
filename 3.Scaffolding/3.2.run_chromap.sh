#!/usr/bin/env bash
set -euo pipefail

# 输入路径
CONTIGS="/path/to/{input}.fa"
HIC_DIR="/path/to/HicFiles"
THREADS=32

# 工具路径
SAMTOOLS="/samtools"

# 输出前缀
OUT_PREFIX="{output}"

echo "Indexing contig assembly"
$SAMTOOLS faidx $CONTIGS
/chromap -i -r $CONTIGS -o ${OUT_PREFIX}.index

    R1="${HIC_DIR}/{HiC}_R1.fastq"
    R2="${HIC_DIR}/{HiC}_R2.fastq"
    
    OUT="${OUT_PREFIX}}"

    echo "  -> Reading ${R1} and ${R2}"
    chromap \
      --preset hic \
      -r $CONTIGS \
      -x ${OUT_PREFIX}.index \
      --remove-pcr-duplicates \
      -1 $R1 \
      -2 $R2 \
      --SAM \
      -o ${OUT}.sam \
      -t $THREADS

    echo "  -> Converting SAM to BAM"
    $SAMTOOLS view -bh ${OUT}.sam | \
      $SAMTOOLS sort -@ $THREADS -n -o ${OUT}.bam

    rm ${OUT}.sam
    echo "=== Output written to ${OUT}.bam ==="
done

echo "Convert bed file from Bam file for downstream yahs analysis"

/samtools view -bh -F 0xF0C -q 0 ${OUT}.bam | \
/bedtools bamtobed | \
  awk -v OFS='\t' '{$4=substr($4,1,length($4)-2); print}' \
  > {output}_chromap.bed

echo "Bed file has been generated for downstream yahs analysis"

echo "=== Processing Completed ==="
