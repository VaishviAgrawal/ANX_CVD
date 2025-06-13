library(devtools)
library(ggplot2)
library(reshape2)
library(gridExtra)
require(GenomicSEM)

setwd("write_your_own_path_name")


# --------------------------------------------
# PANEL (A) FULL MODEL
# --------------------------------------------
#write the name of the files according to your own data
traits_full <- c("ANX.sumstats"  ,    #Binary
  "CAD.sumstats",      #Binary
  "STR.sumstats",        #Binary
  "HF.sumstats",          #Binary
  "BrS.sumstats",       #Binary
  "HYP_2019.sumstats", #Binary"
  "AF.sumstats",
  "QT.sumstats",
  "JT.sumstats",
  "QRS.sumstats",
  "PR.sumstats",
  "brady_snd_inc.sumstats",
  "brady_snd_rest.sumstats",
  "brady_dist_inc.sumstats",
  "brady_dist_rest.sumstats",
  "brady_pacer.sumstats"
  
)

#Prevalence
#writ the prev accordingly
sample.prev_full <- c( 0.07981,
                      0.123, 
                      0.056, 
                      0.082, 
                      0.219,
                      0.316,
                      0.112,
                      NA,NA,NA,NA,
                      0.0079,
                      0.0041,
                      0.0349,
                      0.0070,
                      0.0235
                     )
population.prev_full <- c(0.05,0.036, 0.092, 0.017, 0.0001,0.22,0.02,NA,NA,NA,NA,0.038,0.038,0.038,0.038,0.038)
trait.names_full <- c("ANX","CAD","STR", "HF","BrS","HYP","AF","QT","JT","QRS","PR","Brdy_SND_inc","Brdy_SND_res","Brdy_DIST_inc","Brdy_DIST_res","Brdy_pacer")


LDSC_full <-ldsc(
  traits = traits_full,
  sample.prev = sample.prev_full,
  population.prev = population.prev_full,
  ld = "eur_w_ld_chr/",
  wld = "eur_w_ld_chr/",
  trait.names = trait.names_full
)

save(LDSC_full, file = "LDSCoutput_full.RData")
load("LDSCoutput_full.RData")

#OPTION-1
#THIS WILL fill only the lower triangle of the SE_full matrix.
m <- nrow(LDSC_full$S)
#SE_full <- matrix(0,m,m)
#SE_full[lower.tri(SE_full, diag=TRUE)] <- sqrt(diag(LDSC_full$V))
#Z_full <- LDSC_full$S / SE_full 
#diag(Z_full)  

#OPTION-2
#THIS WILL GIVE A SYMMETRICAL VIEW IN THE HEATMAP
SE_full <- matrix(0, m, m)
SE_full[lower.tri(SE_full, diag = TRUE)] <- sqrt(diag(LDSC_full$V))
SE_full <- SE_full + t(SE_full) - diag(diag(SE_full))  # mirror lower triangle to upper
Z_full <- LDSC_full$S / SE_full
cor_full <- cov2cor(LDSC_full$S)
p_full <- 2 * (1- pnorm(abs(Z_full)))

#to get ANX as first trait
trait.names_full <- c("ANX", setdiff(trait.names_full, "ANX"))

plot_matrix <- function(mat, signif_criteria, title, trait.names, digits = 2) {
  # Generate label text with specific number of decimals, no rounding for actual data
lab_matrix <- sprintf(paste0("%.", digits, "f"), mat)
lab_matrix[is.infinite(mat)] <- "Inf" 
  
  
highlight <- signif_criteria(mat)
  lab_matrix[highlight] <- paste0(lab_matrix[highlight], "*")
  
df <- data.frame(
    Var1 = rep(trait.names, each = length(trait.names)),
    Var2 = rep(trait.names, times = length(trait.names)),
    value = as.vector(mat),    # keep original matrix values here
    label = as.vector(lab_matrix),
    highlight = as.vector(highlight)
  )
  df$Var1 <- factor(df$Var1, levels = trait.names)
  df$Var2 <- factor(df$Var2, levels = trait.names)
  
  ggplot(df, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 3) +
    #geom_text(aes(label = label, fontface = ifelse(highlight, "bold", "plain")), size = 3) +
    scale_fill_gradient(low = "white", high = "brown", na.value = "white") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(title = title, x = "", y = "", fill = "Value")
}

# Step 3: Define significance rules
z_signif <- function(mat) abs(mat) > 4
cor_signif <- function(mat) abs(mat) > 0.10
p_signif <- function(mat) mat < 0.05
cov_signif <- function(mat) matrix(FALSE, nrow(mat), ncol(mat))

# Step 4: Plot all matrices
p1 <- plot_matrix(Z_full, z_signif, "Z-Score Heatmap", trait.names_full, digits = 2)
p2 <- plot_matrix(cor_full, cor_signif, "Correlation Heatmap", trait.names_full, digits = 2)
p3 <- plot_matrix(LDSC_full$S, cov_signif, "Covariance Heatmap", trait.names_full, digits = 2)
p4 <- plot_matrix(p_full, p_signif, "P-Value Heatmap", trait.names_full, digits = 2)
#grid.arrange(p1, p2,p3,p4 ,ncol = 2)
grid.arrange(p1,p2, ncol=2)
