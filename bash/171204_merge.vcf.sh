
GATK="java -Xmx200g -jar /home/yanjun/soft/gatk/3.7/GenomeAnalysisTK.jar"
REFERENCE="/home/yanjun/projects/AIL_seq/ref.mette/doublechecked/galgal5_usethis.fa"


ls | grep "informative.vcf.recode.vcf$" > vcfs.list
$GATK  -T CombineVariants -R $REFERENCE --variant vcfs.list -genotypeMergeOptions UNIQUIFY  -o output.vcf
$GATK  -T CombineVariants -R $REFERENCE --variant /tmp/nas/yanjun/AIL/founder.hap/mette/AILcohort_0620fil_redo_selected.vcf --variant ./output.vcf -genotypeMergeOptions UNIQUIFY  -o 171115_output.vcf
# once it is doen scp back to my own computer
scp ./171115_output.vcf yanjunzan@130.238.46.184:~/Documents/impute/data/

# 171215
ls | grep "all.vcf.recode.vcf$" > vcfs.all.list
<<<<<<< Updated upstream
$GATK  -T CombineVariants -R $REFERENCE --variant vcfs.all.list -genotypeMergeOptions UNIQUIFY  -o /home/yanjun/projects/F2_seq/data/171215.all.780.vcf
$GATK  -T CombineVariants -R $REFERENCE --variant /tmp/nas/yanjun/AIL/founder.hap/mette/AILcohort_0620fil_redo_selected.vcf --variant /home/yanjun/projects/F2_seq/data/171215.all.780.vcf \ 
-genotypeMergeOptions UNIQUIFY  -o /home/yanjun/projects/F2_seq/data/171215_all.780.F0.output.vcf
=======
$GATK  -T CombineVariants -R $REFERENCE --variant vcfs.all.list -genotypeMergeOptions UNIQUIFY  -o /home/yanjun/171215.all.780.vcf
$GATK  -T CombineVariants -R $REFERENCE --variant /tmp/nas/yanjun/AIL/founder.hap/mette/AILcohort_0620fil_redo_selected.vcf --variant /home/yanjun/171215.all.780.vcf -genotypeMergeOptions UNIQUIFY  -o /home/yanjun/171215_all.780.F0.output.vcf
>>>>>>> Stashed changes
# gz it 
# extract the fixed sites to reduece computing time.
POS="/tmp/nas/yanjun/AIL/1711recall.vcf/171102.informative.siteschr1.txt"

vcftools --vcf 171215_all.780.F0.output.vcf --positions $POS --recode --out 171215_all.780.F0.output

gzip 171215_all.780.F0.output.vcf
#once it is doen scp back to my own computer
scp ./171115_output.vcf yanjunzan@130.238.46.184:~/Documents/impute/data/
<<<<<<< Updated upstream
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> a few files has problematic index, update it. 
GATK="java -Xmx200g -jar /home/yanjun/soft/gatk/3.7/GenomeAnalysisTK.jar"
REFERENCE="/home/yanjun/projects/AIL_seq/ref.mette/doublechecked/galgal5_usethis.fa"
KNOWN="/home/yanjun/projects/AIL_seq/ref.mette/Gallus_gallus_incl_cons_sorted.vcf"

=======

### kept crash due to some individual do not have the right index
GATK="java -Xmx200g -jar /home/yanjun/soft/gatk/3.7/GenomeAnalysisTK.jar"
REFERENCE="/home/yanjun/projects/AIL_seq/ref.mette/doublechecked/galgal5_usethis.fa"
KNOWN="/home/yanjun/projects/AIL_seq/ref.mette/Gallus_gallus_incl_cons_sorted.vcf"
>>>>>>> Stashed changes
for i in `cat buildindex.txt`; do
$GATK  -T ValidateVariants -R $REFERENCE -V $i 
done
## those two still do not have right index
GATK="java -Xmx200g -jar /home/yanjun/soft/gatk/3.7/GenomeAnalysisTK.jar"
REFERENCE="/home/yanjun/projects/AIL_seq/ref.mette/doublechecked/galgal5_usethis.fa"
$GATK  -T ValidateVariants -R $REFERENCE -V F2_1105_S41_L004.sorted_dups_marked.recall.171101_all.vcf.recode.vcf
<<<<<<< Updated upstream
$GATK  -T ValidateVariants -R $REFERENCE -V F2_682_S48_L002.sorted_dups_marked.recall.171101_all.vcf.recode.vcf
=======
$GATK  -T ValidateVariants -R $REFERENCE -V F2_1105_S41_L004.sorted_dups_marked.recall.171101_all.vcf.recode.vcf


>>>>>>> Stashed changes

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> the remaining 200 F2

GATK="java -Xmx200g -jar /home/yanjun/soft/gatk/3.7/GenomeAnalysisTK.jar"
REFERENCE="/home/yanjun/projects/AIL_seq/ref.mette/doublechecked/galgal5_usethis.fa"
ls /tmp/nas/yanjun/AIL/171220.final/vcffile/ | grep "informative.vcf.recode.vcf$" > vcfs.all.list
$GATK  -T CombineVariants -R $REFERENCE --variant vcfs.all.list -genotypeMergeOptions UNIQUIFY  -o /home/yanjun/projects/F2_seq/data/180109.all.223.f2.vcf
$GATK  -T CombineVariants -R $REFERENCE --variant /tmp/nas/yanjun/AIL/founder.hap/mette/AILcohort_0620fil_redo_selected.vcf --variant /home/yanjun/projects/F2_seq/data/180109.all.223.f2.vcf -genotypeMergeOptions UNIQUIFY  -o /home/yanjun/projects/F2_seq/data/180109.all.223.f2.P60.vcf
# join the 700+P60 and 200
$GATK  -T CombineVariants -R $REFERENCE --variant /home/yanjun/projects/F2_seq/data/171215_all.780.F0.output.vcf --variant /home/yanjun/projects/F2_seq/data/180109.all.223.f2.vcf -genotypeMergeOptions UNIQUIFY  -o /home/yanjun/projects/F2_seq/data/180109.all.223+700.f2.P60.vcf
## extract a few postion from the final vcf

POS="/tmp/nas/yanjun/AIL/1711recall.vcf/171102.informative.siteschr1.txt"

vcftools --vcf /home/yanjun/projects/F2_seq/data/180109.all.223.f2.P60.vcf --positions $POS --recode --out /home/yanjun/projects/F2_seq/data/180109.all.223.f2.P60.informative
gzip 180109.all.223.f2.P60.informative.recode.vcf
