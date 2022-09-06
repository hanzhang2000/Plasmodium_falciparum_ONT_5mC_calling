#!/bin/bash
# index fast5 files

TARGETDIR="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/Pfalciparum_23062022/barcode${SLURM_ARRAY_TASK_ID}" 
QUERY="${TARGETDIR}/reads_barcode${SLURM_ARRAY_TASK_ID}.fastq" 
Nanopolish="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/nanopolish/nanopolish"

$Nanopolish index --verbose -d /home/hz395/rds/rds-partition_4-KWID6wBsuFo/Pfalciparum_ONT_data/2021_11_01_FT_ONT_Plasmodium_BrdU_EdU_CellCycle/fast5/fast5/ -s /home/hz395/rds/rds-partition_4-KWID6wBsuFo/Pfalciparum_ONT_data/2021_11_01_FT_ONT_Plasmodium_BrdU_EdU_CellCycle/fast5/sequencing_summary_FAQ47712_762588b8.txt $QUERY
