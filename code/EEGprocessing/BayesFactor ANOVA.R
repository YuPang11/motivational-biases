## ============================================================
## 0) Load required R packages
##    If already installed, you can directly use library()
## ============================================================
pkgs <- c("BayesFactor", "dplyr", "tidyr", "readr", "afex")
to_install <- setdiff(pkgs, rownames(installed.packages()))
if (length(to_install) > 0) {
  install.packages(to_install, dependencies = TRUE)
}

library(BayesFactor)
library(dplyr)
library(tidyr)
library(readr)
library(afex)

set.seed(123)  # Make BayesFactor results reproducible
afex::set_sum_contrasts()  # Contrast setting (more consistent with SPSS Type III)

## ============================================================
## 1) Read the same wide-format data used in SPSS
## ============================================================
fb_path <- "H:/MT_all(afterqujizhi)/cueoriandbu/feedback/TF/OCD/fb4spss.csv"
dat_wide <- read_csv(fb_path)
print(colnames(dat_wide))  # Check column names

## ============================================================
## 2) Build a 2×2 within-subject design:
##    Outcome (Win vs Pun) × Response (Go vs NoGo)
##    Using four columns: goWin, nogoWin, goPun, nogoPun
## ============================================================
dat_long <- dat_wide %>%
  select(subject, goPun, nogoPun, goWin, nogoWin) %>%
  pivot_longer(
    cols      = -subject,
    names_to  = "cond",
    values_to = "theta"           # Midfrontal theta value
  ) %>%
  mutate(
    subject = factor(subject),
    
    ## Outcome factor: Win vs Pun
    outcome = case_when(
      cond %in% c("goWin",  "nogoWin")  ~ "Win",
      cond %in% c("goPun",  "nogoPun")  ~ "Pun",
      TRUE ~ NA_character_
    ),
    
    ## Response factor: Go vs NoGo
    response = case_when(
      cond %in% c("goWin",  "goPun")    ~ "Go",
      cond %in% c("nogoWin","nogoPun")  ~ "NoGo",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    outcome  = factor(outcome,  levels = c("Win", "Pun")),
    response = factor(response, levels = c("Go",  "NoGo"))
  )

print(head(dat_long))

## ============================================================
## 3) Classical repeated-measures ANOVA (frequentist F and p)
##    Corresponds to SPSS:
##    GLM goPun nogoPun goWin nogoWin
##      /WSFACTOR=outcome 2 response 2 ...
## ============================================================
anova_res <- aov_ez(
  id     = "subject",
  dv     = "theta",
  within = c("outcome", "response"),
  data   = dat_long,
  type   = 3              # Type III, aligned with SPSS
)

tab <- as.data.frame(anova_res$anova_table)
print(tab)

## ---- F and p for the Outcome × Response interaction
##      (biased vs unbiased learning) ----
F_int   <- as.numeric(tab["outcome:response", "F"])
p_int   <- as.numeric(tab["outcome:response", "p.value"])
df1_int <- as.numeric(tab["outcome:response", "num Df"])
df2_int <- as.numeric(tab["outcome:response", "den Df"])

## ---- F and p for the Response main effect
##      (strong vs weak learning) ----
F_resp   <- as.numeric(tab["response", "F"])
p_resp   <- as.numeric(tab["response", "p.value"])
df1_resp <- as.numeric(tab["response", "num Df"])
df2_resp <- as.numeric(tab["response", "den Df"])

## ============================================================
## 4) Bayesian repeated-measures ANOVA (Bayes factor BF01)
##    Model: theta ~ outcome * response + subject
## ============================================================
bf_full <- anovaBF(
  formula     = theta ~ outcome * response + subject,
  data        = dat_long,
  whichRandom = "subject"
)

print(bf_full)
print(names(bf_full))  # Check model names for indexing below

## ---- BF01 for the Outcome × Response interaction
##      (biased vs unbiased learning) ----
## With interaction: outcome + response + outcome:response + subject
bf_with_int <- bf_full["outcome + response + outcome:response + subject"]

## Without interaction: outcome + response + subject
bf_no_int   <- bf_full["outcome + response + subject"]

## Divide the two models to get BF10 for "with interaction vs without interaction"
bf_int <- bf_with_int / bf_no_int

## extractBF() gets the numeric BF from a BFBayesFactor object
bf10_int_num <- extractBF(bf_int, logbf = FALSE)$bf
bf01_int     <- 1 / bf10_int_num

## ---- BF01 for the Response main effect
##      (strong vs weak learning) ----
## With response main effect: outcome + response + subject
bf_with_resp <- bf_full["outcome + response + subject"]

## Without response main effect: outcome + subject
bf_no_resp   <- bf_full["outcome + subject"]

bf_resp       <- bf_with_resp / bf_no_resp
bf10_resp_num <- extractBF(bf_resp, logbf = FALSE)$bf
bf01_resp     <- 1 / bf10_resp_num

## ============================================================
## 5) Output two result lines: F, p, BF01
##    1) Biased vs unbiased learning (Outcome × Response)
##    2) Strong vs weak learning (Response main effect)
## ============================================================
cat(sprintf(
  "Biased vs unbiased learning (Outcome × Response): F(%d,%d) = %.3f, p = %.3f, BF01 = %.3f\n",
  df1_int, df2_int, F_int, p_int, bf01_int
))

cat(sprintf(
  "Strong vs weak learning (Response main effect): F(%d,%d) = %.3f, p = %.3f, BF01 = %.3f\n",
  df1_resp, df2_resp, F_resp, p_resp, bf01_resp
))