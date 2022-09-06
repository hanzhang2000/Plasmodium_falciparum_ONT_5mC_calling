#!/bin/bash
# Nanopolish call methylation

GENOME="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/Pfalciparum_ONT_data/Plasmodium_falciparum.ASM276v2.fasta"
TARGETDIR="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_11082022" 
QUERY="${TARGETDIR}/reads_barcode01.fastq"
Nanopolish="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/nanopolish/nanopolish"
OUTPUT="${TARGETDIR}/methylation_calls.tsv"

$Nanopolish call-methylation --threads 76 --reads $QUERY --bam ${TARGETDIR}/alignments.sorted.bam --genome $GENOME > $OUTPUT