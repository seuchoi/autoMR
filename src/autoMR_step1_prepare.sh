#!/bin/bash

while getopts ":r:e:f:v:p:o:m:i:d:b:c:s:" option
do
case "${option}"
in
r) reffile=${OPTARG};;
e) expname=${OPTARG};;
f) expfile=${OPTARG};;
v) expvarid=${OPTARG};;
p) exppvalue=${OPTARG};;
o) outcome=${OPTARG};;
m) outfile=${OPTARG};;
i) outvarid=${OPTARG};;
b) kbp=${OPTARG};;
c) cpval=${OPTARG};;
s) rsq=${OPTARG};;
esac
done

echo "Reference plink file: ${reffile}"
echo "Exposure name: ${expname}"
echo "Exposure summary statistics: ${expfile}"
echo "Exposure variant ID(chr:pos:ref:alt) column number: ${expvarid}"
echo "Exposure p-value column number: ${exppvalue}"
echo "Outcome name: ${outcome}"
echo "Outcome summary staitsitcs: ${outfile}"
echo "Outcome variant ID(chr:pos:ref:alt) column number: ${outvarid}"
echo "Clumping Kbp: ${kbp}"
echo "Clumping pvalue cutoff: ${cpval}"
echo "Clumping Rsquare: ${rsq}"


####
echo "prepare for MR start"

#### download plink
wget http://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip
unzip -o plink_linux_x86_64_20201019.zip

#### 1000G PLINK file
#varid1000g=varid1000g.txt
awk '{print $2}' ${reffile}.bim > refvarid.txt

#### remove ambiguous snps
egrep -v ":A:T|:T:A|:G:C|:C:G" refvarid.txt > refvarid2.txt
mv refvarid2.txt refvarid.txt

#### LA GWAS
#expvaradd=${expname}_maf0.01_add.tsv
#expvarid=${expname}_maf0.01_varid.tsv
#awk '($7>=0.01)&&($7<=0.99)&&($8>=0.3){print $0"\t"$2":"$3":"$5":"$6}' ${expfile} > ${expname}_maf0.01_add.tsv
awk -F '\t' -v varid="$expvarid" '(NR>1){print $varid}' ${expfile} > ${expname}_varid.tsv

#### AFib GWAS
#outvarid=${outcome}_maf0.01_varid.tsv
#awk '($10>=0.01)&&($10<=0.99)&&($12>=0.3) {print $13}' ${outfile} > ${outcome}_maf0.01_varid.tsv
awk -F '\t' -v varid="$outvarid" '(NR>1){print $varid}' ${outfile} > ${outcome}_varid.tsv

#### sort files
sort -o refvarid.txt refvarid.txt
sort -o ${expname}_varid.tsv ${expname}_varid.tsv
sort -o ${outcome}_varid.tsv ${outcome}_varid.tsv

#### variants in common
comm -12 refvarid.txt ${expname}_varid.tsv > ${expname}_ref_comm1.txt
var1=$(wc -l < ${expname}_ref_comm1.txt)
echo "${var1} variants are in both Reference and ${expname}"
comm -12 ${expname}_ref_comm1.txt ${outcome}_varid.tsv > ${outcome}_${expname}_ref_comm2.txt
var2=$(wc -l < ${outcome}_${expname}_ref_comm2.txt)
echo "${var2} variants are in Reference, ${expname}, ${outcome}"

#### matching varids
#### grep -Fwf comm2.txt ${expvaradd} | head
awk -F '\t' -v varid="$expvarid" 'NR==FNR {id[$1]; next} $varid in id' ${outcome}_${expname}_ref_comm2.txt ${expfile} > ${expname}_match.tsv
awk -F '\t' -v varid2="$outvarid" 'NR==FNR {id[$1]; next} $varid2 in id' ${outcome}_${expname}_ref_comm2.txt ${outfile} > ${outcome}_match.tsv

#### clumping prepare for exposure
#echo -e "SNP A1 A2 freq beta se P" > ${expname}_clump_input.tsv
#awk '{print $17,$5,$6,$7,$11,$12,$10}' ${expname}_match.tsv >> ${expname}_clump_input.tsv

echo -e "SNP P" > ${expname}_clump_input.tsv
awk -v varid="$expvarid" -v pval="$exppvalue" '{print $varid,$pval}' ${expname}_match.tsv >> ${expname}_clump_input.tsv

#### perfrom clumping using plink
for chrom in {1..22}
do
echo "perfrom clumping on ${expname} chr=${chrom} kbp=${kbp}kb pvalue=${cpval} r2=${rsq}"
#/medpop/afib/software/plink1.9/20190304/plink
./plink --silent --bfile ${reffile} \
--extract ${outcome}_${expname}_ref_comm2.txt --chr ${chrom} --clump ${expname}_clump_input.tsv --clump-kb ${kbp} --clump-p1 ${cpval} --clump-r2 ${rsq} --out ${expname}_clump_chr${chrom}_${kbp}kb_${cpval}_r2_${rsq}
done

#### combine all clumping results
clump_files=(${expname}_clump_chr*_${kbp}kb_${cpval}_r2_${rsq}.clumped)
for ((i=0; i < ${#clump_files[@]}; i++))
do
clump_file=${clump_files[i]}
if [ "$i" -eq 0 ];then
head -n -2 ${clump_file} > ${expname}_clump_ALLchr_${kbp}kb_${cpval}_r2_${rsq}.clumped
else
tail -n+2 ${clump_file} | head -n -2 >> ${expname}_clump_ALLchr_${kbp}kb_${cpval}_r2_${rsq}.clumped
fi
echo "${clump_file} is combined"
done

#### seleted variants
clump1=$(wc -l < ${expname}_clump_ALLchr_${kbp}kb_${cpval}_r2_${rsq}.clumped)
echo "$((${clump1}-1)) variants were selected"

#### take the clumped snps
awk '(NR>1){print $3}' ${expname}_clump_ALLchr_${kbp}kb_${cpval}_r2_${rsq}.clumped > ${expname}_clump_ALLchr_${kbp}kb_${cpval}_r2_${rsq}.clumped_varid.txt

#### matching variants in exposure and outcome
head -n 1 ${expfile} > ${expname}_match_clump.tsv
awk -F '\t' -v varid="$expvarid" 'NR==FNR {id[$1]; next} $varid in id' ${expname}_clump_ALLchr_${kbp}kb_${cpval}_r2_${rsq}.clumped_varid.txt ${expname}_match.tsv >> ${expname}_match_clump.tsv

head -n 1 ${outfile} > ${outcome}_match_clump.tsv
awk -F '\t' -v varid2="$outvarid" 'NR==FNR {id[$1]; next} $varid2 in id' ${expname}_clump_ALLchr_${kbp}kb_${cpval}_r2_${rsq}.clumped_varid.txt ${outcome}_match.tsv >> ${outcome}_match_clump.tsv

####
echo "MR analysis is ready"
echo "Please use ${expname}_match_clump.tsv and ${outcome}_match_clump.tsv"

######## MR preparation is done
######## clean the files
rm refvarid.txt
rm ${expname}_varid.tsv
rm ${outcome}_varid.tsv
rm ${expname}_ref_comm1.txt
rm ${outcome}_${expname}_ref_comm2.txt
rm ${expname}_match.tsv
rm ${outcome}_match.tsv
rm ${expname}_clump_input.tsv
rm ${expname}_clump_chr*_${kbp}kb_${cpval}_r2_${rsq}.*
rm ${expname}_clump_ALLchr_${kbp}kb_${cpval}_r2_${rsq}.clumped_varid.txt
rm plink_linux_x86_64_20201019.zip
rm plink
rm toy.map
rm toy.ped
rm prettify
