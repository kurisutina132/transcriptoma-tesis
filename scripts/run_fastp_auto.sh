
#!/bin/bash

cd /mnt/c/Users/kuris/OneDrive/Escritorio/tesis/trancriptoma

mkdir -p trimmed QC_post

for r1 in Reads/Reads/*_R1_001.fastq.gz
do
    base=$(basename ${r1} _R1_001.fastq.gz)
    r2=Reads/Reads/${base}_R2_001.fastq.gz

    echo "Procesando ${base}"

    fastp \
      -i ${r1} \
      -I ${r2} \
      -o trimmed/${base}_R1.trim.fastq.gz \
      -O trimmed/${base}_R2.trim.fastq.gz \
      --detect_adapter_for_pe \
      --cut_mean_quality 20 \
      --length_required 50 \
      --thread 8 \
      --html QC_post/${base}_fastp.html \
      --json QC_post/${base}_fastp.json

done

