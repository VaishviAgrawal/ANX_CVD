setwd("own_path_name")
library(GenomicSEM)

#use clean GWAS summary statistics
files <- c("aragam_2022_CAD_clean.txt",
           "mishra_2022_stroke_clean.txt", 
           "Henry_2025_HF_EUR_clean.txt",
           "Barc_2022_BrS_clean.txt",
           "Zhu_2019_HYP_clean.txt",
           "Friligkou_2024_ANX_EUR_clean.txt"
           )
trait.names <- c("CAD", "STR", "HF", "BrS", "HYP", "ANX")
n_cases <- c(120788, 73652, 153174, 2820,144793, 87517)
n_controls <- c(859531, 1234808, 1793175,7181, 313761, 1008941)
N <- 4 / (1/n_cases + 1/n_controls)
print(N)
#N <- c(980319,1308460,1946349,12821,458554,1096458)
se.logit <- c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE)
linprob <-c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE) 
#ref <- "w_hm3.snplist"
ref <- "reference.1000G.maf.0.005.txt" 

sumstats_out <- sumstats(
  files = files,
  ref = ref,
  trait.names = trait.names,
  N = N,
  se.logit = se.logit,  # Binary traits
  OLS = NULL,
  linprob=linprob,
  betas=NULL,
  info.filter = 0.6,
  maf.filter = 0.01,
  keep.indel=FALSE,
  parallel=FALSE,#chatgptrec:FALSE; but accroding to wiki its TRUE
  cores=NULL,
  )

save(sumstats_out, file = "sumstats_output.RData")
View(sumstats_out)
#write.table(sumstats_out, file = "sumstats_output.RData", quote = FALSE, row.names = FALSE, col.names = TRUE, sep ="\t")
load("sumstats_output.RData")
head(sumstats_out)
colnames(sumstats_out)
