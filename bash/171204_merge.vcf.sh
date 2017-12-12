GATK="java -Xmx10g -jar /home/yanjun/soft/gatk/3.7/GenomeAnalysisTK.jar"
REFERENCE="/home/yanjun/projects/AIL_seq/ref.mette/doublechecked/galgal5_usethis.fa"
ls | grep "informative.vcf.recode.vcf$" > vcfs.list
$GATK  -T CombineVariants -R $REFERENCE --variant vcfs.list -genotypeMergeOptions UNIQUIFY  -o output.vcf
$GATK  -T CombineVariants -R $REFERENCE --variant /tmp/nas/yanjun/AIL/founder.hap/mette/AILcohort_0620fil_redo_selected.vcf --variant ./output.vcf -genotypeMergeOptions UNIQUIFY  -o 171115_output.vcf
# once it is doen scp back to my own computer
scp ./171115_output.vcf yanjunzan@130.238.46.184:~/Documents/impute/data/