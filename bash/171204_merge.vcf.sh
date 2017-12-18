
GATK="java -Xmx200g -jar /home/yanjun/soft/gatk/3.7/GenomeAnalysisTK.jar"
REFERENCE="/home/yanjun/projects/AIL_seq/ref.mette/doublechecked/galgal5_usethis.fa"


ls | grep "informative.vcf.recode.vcf$" > vcfs.list
$GATK  -T CombineVariants -R $REFERENCE --variant vcfs.list -genotypeMergeOptions UNIQUIFY  -o output.vcf
$GATK  -T CombineVariants -R $REFERENCE --variant /tmp/nas/yanjun/AIL/founder.hap/mette/AILcohort_0620fil_redo_selected.vcf --variant ./output.vcf -genotypeMergeOptions UNIQUIFY  -o 171115_output.vcf
# once it is doen scp back to my own computer
scp ./171115_output.vcf yanjunzan@130.238.46.184:~/Documents/impute/data/

# 171215
ls | grep "all.vcf.recode.vcf$" > vcfs.all.list
$GATK  -T CombineVariants -R $REFERENCE --variant vcfs.all.list -genotypeMergeOptions UNIQUIFY  -o /home/yanjun/171215.all.780.vcf
$GATK  -T CombineVariants -R $REFERENCE --variant /tmp/nas/yanjun/AIL/founder.hap/mette/AILcohort_0620fil_redo_selected.vcf --variant ./171215.all.780.vcf -genotypeMergeOptions UNIQUIFY  -o 171215_all.780.F0.output.vcf
# gz it 

gzip 171215_all.780.F0.output.vcf
#once it is doen scp back to my own computer
scp ./171115_output.vcf yanjunzan@130.238.46.184:~/Documents/impute/data/

GATK="java -Xmx200g -jar /home/yanjun/soft/gatk/3.7/GenomeAnalysisTK.jar"
REFERENCE="/home/yanjun/projects/AIL_seq/ref.mette/doublechecked/galgal5_usethis.fa"
KNOWN="/home/yanjun/projects/AIL_seq/ref.mette/Gallus_gallus_incl_cons_sorted.vcf"


for i in `cat buildindex.txt`; do
$GATK  -T ValidateVariants -R $REFERENCE -V $i 
done

GATK="java -Xmx200g -jar /home/yanjun/soft/gatk/3.7/GenomeAnalysisTK.jar"
REFERENCE="/home/yanjun/projects/AIL_seq/ref.mette/doublechecked/galgal5_usethis.fa"

$GATK  -T ValidateVariants -R $REFERENCE -V F2_1105_S41_L004.sorted_dups_marked.recall.171101_all.vcf.recode.vcf

$GATK  -T ValidateVariants -R $REFERENCE -V F2_682_S48_L002.sorted_dups_marked.recall.171101_all.vcf.recode.vcf



