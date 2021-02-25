# autoMR
## two sample Mendelian Randomization can be excecuted using two simple commandlines
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
