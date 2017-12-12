#! /usr/local/bin/bash

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>This is a pipe line from BMA to vcf using GATK unified genotyper<<<<<<<<<<<<<<<<<<<<<<<<<<<
## input files are passed in through bash

# input specification
core=1
pic_path="/home/yanjun/soft/picard.jar"
SAMTOOLS="samtools"
GATK="java -Xmx10g -jar /home/yanjun/soft/gatk/3.7/GenomeAnalysisTK.jar"
REFERENCE="/home/yanjun/projects/AIL_seq/ref.mette/doublechecked/galgal5_usethis.fa"
outdir="/tmp/nas/yanjun/AIL/1711recall.vcf"
KNOWN="/home/yanjun/projects/AIL_seq/ref.mette/doublechecked/dbSNP_from.AIL.1.vcf"
POS="/tmp/nas/yanjun/AIL/1711recall.vcf/171102.informative.siteschr1.txt"
POS2="/tmp/nas/yanjun/AIL/1711recall.vcf/1711.1keep.txt"
if [ ! -d "$outdir"/""vcffile"" ];then
mkdir -p "$outdir"/""vcffile""
fi
cd  "$outdir"/""vcffile""

prefix=`basename $1 `     
 $GATK       -T UnifiedGenotyper \
                -R $REFERENCE \
                -nt $core \
                -I $1 \
                --output_mode EMIT_ALL_SITES \
                -o "${prefix%.bam}.recall.1.171101.vcf" 
vcftools --vcf "${prefix%.bam}.recall.1.171101.vcf" --positions $POS2 --recode --out "${prefix%.bam}.recall.171101_all.vcf"
vcftools --vcf "${prefix%.bam}.recall.171101_all.vcf.recode.vcf" --positions $POS --recode --out "${prefix%.bam}.recall.171101_informative.vcf"
rm "${prefix%.bam}.recall.1.171101.vcf"

    
