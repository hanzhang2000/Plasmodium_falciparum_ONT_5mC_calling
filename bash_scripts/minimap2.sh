#!/bin/bash

module load ceuadmin/samtools/1.2

GENOME="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/Pfalciparum_ONT_data/Plasmodium_falciparum.ASM276v2.fasta"
TARGETDIR="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_11082022" 
QUERY="${TARGETDIR}/reads_barcode01.fastq"
MINIMAP2="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/minimap2-2.17_x64-linux/minimap2"
OUTSAM="${TARGETDIR}/alignments.sam"
OUTBAM="${TARGETDIR}/alignments.bam"
OUTPREF="${TARGETDIR}/alignments.sorted"

$MINIMAP2 -L -ax map-ont -t 152 -a -o $OUTSAM $GENOME $QUERY
samtools view -Sb -o $OUTBAM $OUTSAM
samtools sort -@ 152 $OUTBAM $OUTPREF
samtools index ${OUTPREF}.bam