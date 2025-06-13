setwd("own_path_name")

# Load required packages
library(dplyr)
library(ggplot2)

# Set your path to bivariate LAVA results
lava_folder <- "own_path"
lava_files <- list.files(lava_folder, pattern = "\\.lava$", full.names = TRUE)

# Function to read and format each .lava file
read_lava <- function(file_path) {
  df <- read.table(file_path, header = TRUE)
  colnames(df) <- c("locus", "chr", "start", "stop", "n_snps", "n_pcs",
                    "phen1", "phen2", "rho", "rho_lower", "rho_upper",
                    "r2", "r2_lower", "r2_upper", "p")
  # Create a label for facet titles
  df$trait_pair <- paste0(df$phen1[1], " vs ", df$phen2[1])
  return(df)
}

# Read and combine all .lava files
lava_all <- bind_rows(lapply(lava_files, read_lava))

View (lava_all)

# Extract just the CVD trait from phen2
lava_all$trait_label <- lava_all$phen2

lava_all$trait_group <- case_when(
  lava_all$trait_label %in% c("CAD", "HYP", "HF", "STR") ~ "Atherosclerosis",
  lava_all$trait_label %in% c("QT", "PR", "JT", "QRS", "AF", "Brdy_pacer", "Brdy_SND_inc", "Brdy_SND_res", "Brdy_DIST_inc", "Brdy_DIST_res", "BrS") ~ "Rhythm",
  TRUE ~ "Other"
)

# Apply BH FDR correction within each trait pair
lava_all <- lava_all %>%
  group_by(trait_pair) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup()

fdr_lines <- lava_all %>%
  filter(p_adj < 0.05) %>%
  group_by(trait_label) %>%
  summarise(fdr_thresh = max(p)) %>%
  mutate(y_pos = -log10(fdr_thresh))

# Merge back into main data
#lava_all <- left_join(lava_all, fdr_lines, by = "trait_pair") %>%
 # mutate(significant = p_adj < 0.05)
lava_all <- lava_all %>%
  mutate(significant = p_adj < 0.05)


ggplot(lava_all, aes(x = rho, y = -log10(p), color = trait_group)) +
  geom_point(aes(shape = significant), size = 1.2, alpha = 0.8) +
  geom_hline(data = fdr_lines, aes(yintercept = y_pos), linetype = "dashed", color = "black") +  # ADD the '+' here
  geom_vline(xintercept = 0, color = "black") +
  facet_wrap(~trait_label, scales = "free") +
  scale_color_manual(values = c("Atherosclerosis" = "orange", 
                                "Rhythm" = "dodgerblue", 
                                "Other" = "grey70")) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
  theme_bw() +
  labs(title = "Local Genetic Correlation",
       x = "Local Genetic Correlation (ρ)",
       y = expression(-log[10](italic(p))))

#------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)

# Step 1: Filter for Anxiety as one phenotype
anx_data <- lava_all %>%
  filter(phen1 == "ANX")  # or phen2 == "ANX" if that's where ANX is

# Step 2: Get loci with ≥2 significant traits
sig_loci <- anx_data %>%
  filter(p_adj < 0.05) %>%
  group_by(locus) %>%
  summarise(n_sig_traits = n_distinct(trait_label)) %>%
  filter(n_sig_traits >= 2) %>%
  pull(locus)

# Step 3: Keep all traits at those loci
anx_subset <- anx_data %>%
  filter(locus %in% sig_loci) %>%
  mutate(is_significant = p_adj < 0.05)

# Step 4: Plot with color + shape
ggplot(anx_subset, aes(x = factor(locus), y = rho, color = trait_label, shape = is_significant)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "longdash", color = "black") +
  scale_color_manual(
    name = "Cardio_trait",
    values = c(
      "AF" = "#E41A1C","JT" = "cornflowerblue", "PR" = "#FFFF33",
      "QRS" = "#A65628", "QT" = "deeppink", "BrS" = "chartreuse",
      "Brdy_pacer" = "#66C2A5", "Brdy_SND_incl" = "darkgoldenrod",
      "Brdy_SND_res" = "#8DA0CB", "Brdy_DIST_incl" = "darkkhaki",
      "Brdy_DIST_res" = "darkblue", "CAD" = "cyan", "HF" = "forestgreen",
      "HYP" = "darkorchid", "STR" = "darkorange"
    )
  ) +
  scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 1), name = "Significant") +
  labs(
    title = "Regions with ≥2 Cardio Traits Significantly Correlated with Anxiety",
    x = "Locus",
    y = "Correlation"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold")
  )

##########
#------------------------FULL SCRIPT: FDR vs Raw P Comparison------------------------
rm(list=ls()) 
setwd("D:/LAVA_results")

# Load required packages
library(dplyr)
library(ggplot2)

# Set your path to bivariate LAVA results
lava_folder <- "D:/LAVA_results"
lava_files <- list.files(lava_folder, pattern = "\\.lava$", full.names = TRUE)

# Function to read and format each .lava file
read_lava <- function(file_path) {
  df <- read.table(file_path, header = TRUE)
  colnames(df) <- c("locus", "chr", "start", "stop", "n_snps", "n_pcs",
                    "phen1", "phen2", "rho", "rho_lower", "rho_upper",
                    "r2", "r2_lower", "r2_upper", "p")
  df$trait_pair <- paste0(df$phen1[1], " vs ", df$phen2[1])
  return(df)
}

# Read and combine all .lava files
lava_all <- bind_rows(lapply(lava_files, read_lava))

# Trait group categorization
lava_all$trait_label <- lava_all$phen2
lava_all$trait_group <- case_when(
  lava_all$trait_label %in% c("CAD", "HYP", "HF", "STR") ~ "Atherosclerosis",
  lava_all$trait_label %in% c("QT", "PR", "JT", "QRS", "AF", "Brdy_pacer", "Brdy_SND_inc", "Brdy_SND_res", "Brdy_DIST_inc", "Brdy_DIST_res", "BrS") ~ "Rhythm",
  TRUE ~ "Other"
)

# Apply BH FDR correction within each trait pair
lava_all <- lava_all %>%
  group_by(trait_pair) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup()

# Determine significance by adjusted and raw p-values
lava_all <- lava_all %>%
  mutate(significant_fdr = p_adj < 0.05,
         significant_raw = p < 0.05)

# Dashed lines for FDR thresholds
fdr_lines <- lava_all %>%
  filter(p_adj < 0.05) %>%
  group_by(trait_label) %>%
  summarise(fdr_thresh = max(p)) %>%
  mutate(y_pos = -log10(fdr_thresh))

# -------- Plot A: FDR-corrected significance --------
ggplot(lava_all, aes(x = rho, y = -log10(p), color = trait_group)) +
  geom_point(aes(shape = significant_fdr), size = 1.2, alpha = 0.8) +
  geom_hline(data = fdr_lines, aes(yintercept = y_pos), linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  facet_wrap(~trait_label, scales = "free") +
  scale_color_manual(values = c("Atherosclerosis" = "darkorange", 
                                "Rhythm" = "deepskyblue2", 
                                "Other" = "grey70")) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 1), name = "FDR < 0.05") +
  theme_bw() +
  labs(title = "Local Genetic Correlation (FDR Significance)",
       x = "Local Genetic Correlation (ρ)",
       y = expression(-log[10](italic(p))))

# -------- Plot B: Raw p-value significance --------
ggplot(lava_all, aes(x = rho, y = -log10(p), color = trait_group)) +
  geom_point(aes(shape = significant_raw), size = 1.2, alpha = 0.8) +
  geom_vline(xintercept = 0, color = "black") +
  facet_wrap(~trait_label, scales = "free") +
  scale_color_manual(values = c("Atherosclerosis" = "chocolate1", 
                                "Rhythm" = "cornflowerblue", 
                                "Other" = "grey70")) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 1), name = "Raw p < 0.05") +
  theme_bw() +
  labs(title = "Local Genetic Correlation (Raw p-value Significance)",
       x = "Local Genetic Correlation (ρ)",
       y = expression(-log[10](italic(p))))
