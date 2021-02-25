# autoMR
## Two sample Mendelian Randomization can be excecuted using two simple commandlines
### MR preparation
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
Tab delimieted summary statistics for the exposure variable. VariantID (chr:pos:ref:alt) is required.
#### -v : variantID column in exposure summary statistics
The column number of variantID in exposure summary statistics (ie, 1 or 10)
#### -p : p-value column in exposure summary statistics\
The column number of P-value in exposure summary statistics (ie, 1 or 10)
#### -o : outcome name
It is requiared for intermediate and final file names
#### -m : outcome summary statistics
Tab delimieted summary statistics for the outcome variable. VariantID (chr:pos:ref:alt) is required.
#### -i : variant ID column in outcome summary statistics
The column number of variantID in outcome summary statistics (ie, 1 or 10)
#### -b : clumping distance kbp
It is required for clummping in Plink (500, 1000)
#### -c : clumping p-value  threshold
It is required for clummping in Plink (5e-8, 1e-6)
#### -s : clumping rsquare value
It is required for clummping in Plink (0.3, 0.5)
