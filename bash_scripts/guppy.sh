#!/bin/bash

GUPPY="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/guppy_v6.1.7/ont-guppy/bin/guppy_basecaller"
INPUT="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_11082022/fast5_partition/"
OUTPUT="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_11082022/fast5_basecalled"
CONFIG="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/guppy_v5.0.16/ont-guppy/data/dna_r9.4.1_450bps_fast.cfg"

$GUPPY -i $INPUT -s $OUTPUT -x 'auto' -r -c $CONFIG --cpu_threads_per_caller 32 --barcode_kits "EXP-NBD104"
# --barcode_kits "EXP-NBD104"

