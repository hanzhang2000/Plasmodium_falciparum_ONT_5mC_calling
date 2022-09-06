#!/bin/bash
# make QC plots (read lengths, alignment quality, read quality etc.) of ONT reads using Nanoplot

source /home/hz395/rds/rds-partition_4-KWID6wBsuFo/Pfalciparum_23062022/python3/bin/activate

TARGETDIR="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/Pfalciparum_23062022/barcode0${SLURM_ARRAY_TASK_ID}" #! ${SLURM_ARRAY_TASK_ID}"
QUERY="${TARGETDIR}/alignments.sorted.bam" 
NanoPlot="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/Pfalciparum_23062022/python3/bin/NanoPlot"

python3 $NanoPlot --bam $QUERY --N50 --dpi 300 --outdir $TARGETDIR -t 76