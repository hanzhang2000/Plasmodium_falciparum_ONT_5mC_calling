#!/bin/bash
# summarise methylation information from modified BAM files (Guppy)

module load ceuadmin/samtools/1.2

GENOME="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/Pfalciparum_ONT_data/Plasmodium_falciparum.ASM276v2.fasta"
TARGETDIR="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/output_11082022_Guppy"
MODBAM2BED="/home/hz395/rds/rds-partition_4-KWID6wBsuFo/software/modbam2bed/modbam2bed"

for number in {}
do
samtools merge -@ 152 ${TARGETDIR}/fast5_basecalled/pass/barcode0${SLURM_ARRAY_TASK_ID}/merged.bam ${TARGETDIR}/fast5_basecalled/pass/barcode0${SLURM_ARRAY_TASK_ID}/*.bam
samtools merge -@ 152 ${TARGETDIR}/fast5_basecalled/fail/barcode0${SLURM_ARRAY_TASK_ID}/merged.bam ${TARGETDIR}/fast5_basecalled/fail/barcode0${SLURM_ARRAY_TASK_ID}/*.bam
samtools merge -@ 152  ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/merged.bam ${TARGETDIR}/fast5_basecalled/pass/barcode0${SLURM_ARRAY_TASK_ID}/merged.bam ${TARGETDIR}/fast5_basecalled/fail/barcode0${SLURM_ARRAY_TASK_ID}/merged.bam
samtools sort -@ 152 ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/merged.bam ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/merged.sorted
samtools index ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/merged.sorted.bam
$MODBAM2BED -e -m 5mC --cpg -t 76 $GENOME ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/merged.sorted.bam > ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/guppy.5mC.cpg.bam
$MODBAM2BED -e -m 5mC --chg -t 76 $GENOME ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/merged.sorted.bam > ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/guppy.5mC.chg.bam
$MODBAM2BED -e -m 5mC --chh -t 76 $GENOME ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/merged.sorted.bam > ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/guppy.5mC.chh.bam
$MODBAM2BED -e -m 5hmC --cpg -t 76 $GENOME ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/merged.sorted.bam > ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/guppy.5hmC.cpg.bam
$MODBAM2BED -e -m 5hmC --chg -t 76 $GENOME ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/merged.sorted.bam > ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/guppy.5hmC.chg.bam
$MODBAM2BED -e -m 5hmC --chh -t 76 $GENOME ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/merged.sorted.bam > ${TARGETDIR}/barcode0${SLURM_ARRAY_TASK_ID}/guppy.5hmC.chh.bam
done

echo All done
