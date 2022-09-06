#!/bin/bash
# summarise methylation frequency using helper script

source /home/hz395/rds/rds-partition_4-KWID6wBsuFo/Pfalciparum_23062022/python3/bin/activate

TARGETDIR="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_11082022" 
Nanopolish="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/nanopolish/scripts/calculate_methylation_frequency.py"

$Nanopolish  ${TARGETDIR}/methylation_calls.tsv > ${TARGETDIR}/methylation_frequency.tsv