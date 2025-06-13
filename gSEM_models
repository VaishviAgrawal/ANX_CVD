# ========================
# File: 
# ========================

library(GenomicSEM)
library(DiagrammeR)
library(igraph)

# --- Load LDSC outputs ---
load("LDSCoutput_full.RData")
# -------------------------
# MODEL 1: Panel (a) 
# -------------------------
model1 <- '
F1 =~ NA*CAD +  STR + HF +  BrS + HYP
F2 =~ A*F1 + A*ANX
F2 =~ 1*F2
'
run_model_1 <- usermodel(
  covstruc = LDSC_full,
  estimation = "DWLS",
  model = model1,
  std.lv = TRUE,
  CFIcalc = TRUE,
  imp_cov = FALSE
)
run_model_1
# -------------------------
# MODEL 2: Panel (b) 
# -------------------------
model2 <- '
F1 =~ NA* CAD + STR + HF + HYP
F2 =~ A* HF + A* BrS 
F3 =~ B* F1 + B*ANX
F1 ~~ F2
F4 =~ C*F2 + C*ANX
F3 ~~ 1*F3
F4 ~~ 1*F4
'
run_model_2 <- usermodel(
  covstruc = LDSC_full,
  estimation = "DWLS",
  model = model2,
  std.lv = TRUE,
  CFIcalc = TRUE,
  imp_cov = FALSE
)
summary(run_model_2)
run_model_2$modelfit
run_model_2$results
