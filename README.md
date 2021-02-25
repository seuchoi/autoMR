# autoMR
## Two sample Mendelian Randomization can be excecuted using two simple command lines
R version > 3.6, library(optparse), library(MendelianRandomization), library(data.table) need to be installed

### 1) MR preparation (step 1: before exceuting MR analysis)
./src/autoMR_step1_prepare.sh \
-r # reference panel\
-e # exposure name \
-f # exposure summary statistics\
-v # variant ID column in exposure summary statistics \
-p # p-value column in exposure summary statistics\
-o # outcome name\
-m # outcome summary statistics \
-i # variant ID column in outcome summary statistics\
-b # clumping distance kbp\
-c # clumping p-value  threshold \
-s # clumping rsquare value

#### -r : reference panel
Plink format and variant ID needs to be chr:pos:ref:alt
#### -e : exposure variable name
It is requiared for intermediate and final file names
#### -f : exposure summary statistics
Tab delimieted summary statistics for the exposure variable (header included). VariantID (chr:pos:ref:alt) is required.
#### -v : variantID column in exposure summary statistics
The column number of variantID in exposure summary statistics (ie, 1 or 10)
#### -p : p-value column in exposure summary statistics\
The column number of P-value in exposure summary statistics (ie, 1 or 10)
#### -o : outcome name
It is requiared for intermediate and final file names
#### -m : outcome summary statistics
Tab delimieted summary statistics for the outcome variable (header included). VariantID (chr:pos:ref:alt) is required.
#### -i : variant ID column in outcome summary statistics
The column number of variantID in outcome summary statistics (ie, 1 or 10)
#### -b : clumping distance kbp
It is required for clummping in Plink (500, 1000)
#### -c : clumping p-value  threshold
It is required for clummping in Plink (5e-8, 1e-6)
#### -s : clumping rsquare value
It is required for clummping in Plink (0.3, 0.5)

#### Excecution of this command line will generate three files
1) Clumping result output from Plink
2) Eposure summary statistics with clumped variants (filename will be exposurename_"match_clump.tsv")
3) Outcome summary statistics with clumped variants (filename will be outcomename_"match_clump.tsv")

### 2) MR analysis (step 2: exceuting MR analysis using MendelianRandomization R-package)
Rscript ./src/autoMR_step2_MRtest.R \
--expfile=/path/to/exposure/summarystatistisc/withclumpedvariants/exposurename_match_clump.tsv \
--expvariables= 5 variable names used in the MR analysis  \
--outfile=/path/to/outcome/summarystatistisc/withclumpedvariants/outcomename_match_clump.tsv \
--outvariables=5 variable names used in the MR analysis  \
--ressultfile= result filename

#### --expfile: exposure file
Input file for exposure (tab delimited). Please use the output from step1. The file name will be "yourexposurename"_match_clump.tsv
#### --expvariables: variable names used in the MR analysis in exposure file
5 colum names are reuqires 1) varid 2)beta 3)standard error 4) coded allele 5) noncoded allele.\
ie, --expvariables=varid,BETA,SE,ALLELE1,ALLELE0
#### --outfile: outcome files
Input file for outcome (tab delimited). Please use the output from step1. The file name will be "youroutcomename"_match_clump.tsv
#### --outvariables: variable names used in the MR analysis in outcome file
5 colum names are reuqires 1) varid 2)beta 3)standard error 4) coded allele 5) noncoded allele.\
ie, --outvariables=varid,Effect,StdErr,Allele1,Allele2
#### --ressultfile: result filename
R object will be your result. Please use the yourresultfilename.RData
