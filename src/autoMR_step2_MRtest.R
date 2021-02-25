#### R

library(optparse)
option_list <- list(
  make_option("--expfile", type="character",help="exposure summary statistics",default=NULL),
  make_option("--expvariables", type="character",help="exposure variables: varid(chr:pos:ref:alt), beta, se, effect allele, noneffect allele", default=NULL),
  make_option("--outfile", type="character", help="outcome summary statistics", default=NULL),
  make_option("--outvariables", type="character", help="outcome variables: varid(chr:pos:ref:alt), beta, se, effect allele, noneffect allele",default=NULL),
  make_option("--resultfile", type="character", help="result file name",default="Result.RData")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
if (is.null(opt$expfile)|is.null(opt$expvariables)|is.null(opt$outfile)|is.null(opt$outvariables)){
#if (is.null(opt$expfile)){
  print_help(parser)
  stop("At least one argument must be supplied.", call.=FALSE)
}

#print(opt)


expfile <- opt$expfile
expvariables <- opt$expvariables
outfile <- opt$outfile
outvariables <- opt$outvariables
resultfile <- opt$resultfile

cat("exposure summary statistics file is",expfile,"\n")
cat("the variables from exposure summary statistics file are",expvariables,"\n")
cat("outcome summary statistics file is",outfile,"\n")
cat("the variables from outcome summary statistics file are",outvariables,"\n")

expvariables0<-unlist(strsplit(expvariables,","))
outvariables0<-unlist(strsplit(outvariables,","))

print(expvariables0)
print(outfile)
print(outvariables0)



library(MendelianRandomization)
library(data.table)

exp0<-fread(expfile,header=T,sep="\t",data.table=F)
exp1<-exp0[,c(expvariables0)]
colnames(exp1)<-c("varid","exposure.beta","exposure.se","exposure.allele1","exposure.allele2")

out0<-fread(outfile,header=T,sep="\t",data.table=F)
out1<-out0[,c(outvariables0)]
colnames(out1)<-c("varid","outcome.beta","outcome.se","outcome.allele1","outcome.allele2")

comb1<-merge(exp1,out1,by="varid")
comb1$outcome.beta<-ifelse(comb1$exposure.allele1==comb1$outcome.allele1,comb1$outcome.beta,-comb1$outcome.beta)

comb2<-comb1[,c("varid","exposure.beta","exposure.se","outcome.beta","outcome.se")]
names(comb2)[1]<-"SNP"

comb2$outcome.beta<-ifelse(comb2$exposure.beta<0,-comb2$outcome.beta,comb2$outcome.beta)
comb2$exposure.beta<-ifelse(comb2$exposure.beta<0,-comb2$exposure.beta,comb2$exposure.beta)


MRInputObject <- mr_input(bx = comb2$exposure.beta,
bxse = comb2$exposure.se,
by = comb2$outcome.beta,
byse = comb2$outcome.se)


MRAllObject_all <- mr_allmethods(MRInputObject, method = "all")

MRAllObject_all
save(MRAllObject_all,file=resultfile)

sessionInfo()
quit("no")

######
cd /medpop/afib/schoi/projects/left_atrium/script/

Rscript autoMR_step2_MRtest.R --expfile=/medpop/afib/schoi/projects/left_atrium/data/invnorm_LA_poisson_min_match_clump.tsv \
--expvariables=varid,BETA,SE,ALLELE1,ALLELE0 \
--outfile=/medpop/afib/schoi/projects/left_atrium/data/AF_match_clump.tsv \
--outvariables=varid,Effect,StdErr,Allele1,Allele2 \
--ressultfile=invnorm_LA_poisson_min_AF_MR_results.RData
expfile="/medpop/afib/schoi/projects/left_atrium/data/invnorm_LA_poisson_min_match_clump.tsv"
expvariables="varid,BETA,SE,ALLELE1,ALLELE0"
outfile="/medpop/afib/schoi/projects/left_atrium/data/AF_match_clump.tsv"
outvariables="varid,Effect,StdErr,Allele1,Allele2"
ressultfile="invnorm_LA_poisson_min_AF_MR_results.RData"
