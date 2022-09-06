#!/bin/bash
# Call 5mC in CpG sites 

GUPPY="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/guppy_v6.1.7/ont-guppy/bin/guppy_basecaller"
INPUT="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_11082022/fast5_partition/"
OUTPUT="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_11082022_Guppy/fast5_basecalled"
CONFIG="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/guppy_v6.1.7/ont-guppy/data/dna_r9.4.1_450bps_modbases_5hmc_5mc_cg_hac.cfg"
GENOME='/home/hz395/rds/rds-partition_4-KWID6wBsuFo/Pfalciparum_ONT_data/Plasmodium_falciparum.ASM276v2.fasta'

$GUPPY -i $INPUT -s $OUTPUT -x 'auto' -c $CONFIG --bam_out --recursive --align_ref $GENOME --cpu_threads_per_caller 32 --barcode_kits "EXP-NBD104"
# --barcode_kits "EXP-NBD104"