GATK="java -Xmx20g -jar /home/yanjun/soft/gatk/3.7/GenomeAnalysisTK.jar"
REFERENCE="/home/yanjun/projects/AIL_seq/ref.mette/doublechecked/galgal5_usethis.fa"
KNOWN="/home/yanjun/projects/AIL_seq/ref.mette/Gallus_gallus_incl_cons_sorted.vcf"
$GATK  -T ValidateVariants -R $REFERENCE -V $1