#Here we clean and munge the hypertension summary statistics from a GWAS of 2019
setwd("write_your_path_name")
require(GenomicSEM)
require(data.table)
library(tidyr)
require(tidyverse)

HYP_2019 <- fread("zhu_HT_2019") #write the name of your own file

colnames(HYP_2019)

#check the column names before changing the names
#setnames(HYP_2019, old = 'rsID', new = 'SNP')
#setnames(HYP_2019, old = 'A1', new = 'A1')
setnames(HYP_2019, old = 'A0', new = 'A2')
#setnames(HYP_2019, old = 'Freq1', new = 'EAF')
#setnames(HYP_2019, old = 'Zscore', new = 'Z')
#setnames(HYP_2019, old = 'P-value', new = 'P')
print(paste("HYP_2019 table rename(no difference!):", nrow(HYP_2019)))

#checks what values have NA
colSums(is.na(HYP_2019))

# if there are SNPs with more than one character (ATG etc) in A1
if (any(nchar(HYP_2019$A1) > 1)) {
  # make a new column where those SNPs will be highlighted true
  HYP_2019$A1_multiple <- nchar(HYP_2019$A1) > 1
  # delete those SNPs
  HYP_2019 <- HYP_2019[HYP_2019$A1_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("HYP_2019 table delete A1_multiple:", nrow(HYP_2019)))
  HYP_2019 <- HYP_2019 %>% select(-A1_multiple)
  print(paste("HYP_2019 table delete A1_multiple(no difference!):", nrow(HYP_2019)))
}

# if there are SNPs with more than one character (ATG etc) in A2
if (any(nchar(HYP_2019$A2) > 1)) {
  # make a new column where those SNPs will be highlighted true
  HYP_2019$A2_multiple <- nchar(HYP_2019$A2) > 1
  # delete those SNPs
  HYP_2019 <- HYP_2019[HYP_2019$A2_multiple != TRUE, ]
  # delete the column, there is no more use to it. make sure no changes are made
  print(paste("HYP_2019 table delete A2_multiple:", nrow(HYP_2019)))
  HYP_2019 <- HYP_2019 %>% select(-A2_multiple)
  print(paste("HYP_2019 table delete A2_multiple(no difference!):", nrow(HYP_2019)))
}

print(paste("Number of SNPs with more than one character in A1:", sum(nchar(HYP_2019$A1) > 1)))
print(paste("Number of SNPs with more than one character in A2:", sum(nchar(HYP_2019$A2) > 1)))

# Delete SNPs that are duplicates
HYP_2019 <- HYP_2019[!duplicated(HYP_2019$SNP), ]

#removing SNP with NA
HYP_2019 <- HYP_2019[!is.na(HYP_2019$SNP), ]

sum(is.na(HYP_2019$SNP))  # should return 0

#make MAF
#HYP_2019$EAF <- ifelse(HYP_2019$EAF > 0.5, 1 - HYP_2019$EAF, HYP_2019$EAF)
#setnames(HYP_2019, old = 'EAF', new = 'MAF')

HYP_2019$Z <- HYP_2019$BETA/ HYP_2019$SE


# Estimate N per SNP (if BETA and SE are available)
HYP_2019$N_est <- (HYP_2019$Z * HYP_2019$SE / HYP_2019$BETA)^2

# Remove infinite or NA values due to divide-by-zero
HYP_2019$N_est[!is.finite(HYP_2019$N_est)] <- NA


HYP_2019 <- HYP_2019 %>%
  select(-c(CHR, BP, HWEP, SE, BETA))
print(paste("HYP_2019 table delete columns (no difference!):", nrow(HYP_2019)))

# Check column names before reordering
print(colnames(HYP_2019))

#reorder columns
print(paste("HYP_2019 table reorder:", nrow(HYP_2019)))
HYP_2019 <- HYP_2019[, c('SNP', 'A1', 'A2', 'Z', 'P','N', 'MAF')]
print(paste("HYP_2019 table reorder:(no difference!)", nrow(HYP_2019)))


#save clean file in txt
write.table(HYP_2019, "Zhu_2019_HYP_clean.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

##########################################################################
###################----------MUNGING----------############################
##########################################################################
#-------------------------------------------------------------------------
data <- data.frame(
  file = c("Zhu_2019_HYP_clean.txt"), #name of your own file
  traitnames = c("HYP_2019"),
  N = 458554,  # Use NA if you want GenomicSEM to infer it, or provide actual number
  cases = 144793
)

munge(
  file = data$file,
  hm3 = "w_hm3.snplist",
  trait.names = data$traitnames,
  N = data$N,
  info.filter = 0.9,
  maf.filter = 0.01
)
#file saved as .sumstasts in your directory and you also get a .log file
