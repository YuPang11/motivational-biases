## =====================================================================
## SINGLE-TRIAL MIDFRONTAL THETA / ICPC · OCD group
## - Based on the Swart (2017) framework + simple effects + Excel export
## - Exclude subjects in the OCD group: sID = 201, 218, 232, 238, 243, 256
## =====================================================================

rm(list = ls())

## ------------------ 1. Load required packages ------------------
## Install first if missing: install.packages("xxxx")

library(emmeans)
library(lmerTest)   # Extends lmer (Satterthwaite df, etc.)
library(lme4)       # (g)lmer
library(Matrix)

library(car)        # Anova, qqPlot, scatterplot
library(pbkrtest)
library(lattice)    # densityplot
library(doBy)       # esticon
library(ggplot2)
library(psych)
library(R.matlab)
library(readxl)
if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}
library(writexl)

## simpef.R (kept for compatibility, although esticon is used below)
source("H:/mt_behavalluse/EEGmixmodel/simpef.R")

## ------------------ 2. Read data, keep OCD only, and exclude specified subjects ------------------

data <- as.data.frame(
  read_excel("H:/mt_behavalluse/EEGmixmodel/EEGdatamixmodel_withgroup5.xlsx",
             sheet = "EEGdatamixmodel_withgroup")
)

## Keep OCD group only: group == 1
data <- subset(data, group == 1)

## Exclude subjects with sID = 201, 218, 232, 238, 243, 256
exclude_subs <- c(201, 218, 232, 238, 243, 256)
data <- subset(data, !(sID %in% exclude_subs))

## Convert EEG-related columns to numeric
data$MFpower     <- as.numeric(data$MFpower)
data$ICPCpfc     <- as.numeric(data$ICPCpfc)
data$l_ICPCmotor <- as.numeric(data$l_ICPCmotor)
data$r_ICPCmotor <- as.numeric(data$r_ICPCmotor)

## Build MFpower data frame for compatibility with the original script
## (data already contains only the filtered subjects)
MFpower <- data.frame(
  trial2use    = data$trial2use,
  MFpower      = data$MFpower,
  ICPCpfc      = data$ICPCpfc,
  l_ICPCmotor  = data$l_ICPCmotor,
  r_ICPCmotor  = data$r_ICPCmotor,
  IVleft       = data$IVleft,
  IVright      = data$IVright
)

## Recode predictors to [-1, +1]
data$valence[data$valence == 0] <- -1
data$action[data$action == 0]   <- -1
data$accuracy <- data$DVacc
data$accuracy[data$accuracy == 0] <- -1
data$conflict <- data$valence * data$action

## Make sure trial2use comes from MFpower
data$trial2use <- MFpower$trial2use

## Data structure check (optional)
str(data)
head(data)

## =====================================================================
## 3. SINGLE-TRIAL MIDFRONTAL THETA POWER (MFpower)
## =====================================================================

## ------------------ 3.1 Correct trials: Valence × Action ------------------

## Reassign MFpower in case it was overwritten
data$MFpower <- MFpower$MFpower

## Mean-correct within each subject for correct trials with trial2use == 1
for (iSub in sort(unique(data$sID))) {
  tmp <- data$MFpower[data$sID == iSub & data$DVacc == 1 & data$trial2use == 1]
  tmp <- tmp - mean(tmp, na.rm = TRUE)
  data$MFpower[data$sID == iSub & data$DVacc == 1 & data$trial2use == 1] <- tmp
}

## Plot distribution (optional)
dev.new(); xlab <- "MFpower"
densityplot(~MFpower | as.factor(sID),
            data[data$DVacc == 1 & data$trial2use == 1, ], xlab = xlab)

## Inverse transform to improve normality
data$invPow <- data$MFpower
data$invPow[!is.na(data$MFpower)] <-
  1 - 1 / (data$MFpower[!is.na(data$MFpower)] +
             (1 - min(data$MFpower[!is.na(data$MFpower)])))

dev.new(); xlab <- "invPow"
densityplot(~invPow | as.factor(sID),
            data[data$DVacc == 1 & data$trial2use == 1, ], xlab = xlab)

## Mixed model: invPow ~ action * valence (correct trials)
Mpower_correct <- lmer(
  invPow ~ 1 + action * valence + (1 + action * valence | sID),
  data    = data[data$DVacc == 1 & data$trial2use == 1, ],
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl   = list(maxfun = 196000))
)
summary(Mpower_correct)

## Diagnostics (optional)
densityplot(resid(Mpower_correct))
qqPlot(resid(Mpower_correct))
plot(Mpower_correct, type = c("p", "smooth"))
car::scatterplot(
  fitted(Mpower_correct),
  data$invPow[data$DVacc == 1 & data$trial2use == 1],
  boxplots = FALSE, smooth = FALSE
)

## Wald tests (fixed effects)
esticon(Mpower_correct, diag(length(fixef(Mpower_correct))))

## ------------------ 3.2 Accuracy main effect (correct + error) ------------------

## Mean-correct within each subject across all trials with trial2use == 1
data$MFpower <- MFpower$MFpower

for (iSub in sort(unique(data$sID))) {
  tmp <- data$MFpower[data$sID == iSub & data$trial2use == 1]
  tmp <- tmp - mean(tmp, na.rm = TRUE)
  data$MFpower[data$sID == iSub & data$trial2use == 1] <- tmp
}

dev.new(); xlab <- "MFpower"
densityplot(~MFpower | as.factor(sID),
            data[data$trial2use == 1, ], xlab = xlab)

## Inverse transform
data$invPow <- data$MFpower
data$invPow[!is.na(data$MFpower)] <-
  1 - 1 / (data$MFpower[!is.na(data$MFpower)] +
             (1 - min(data$MFpower[!is.na(data$MFpower)])))

dev.new(); xlab <- "invPow"
densityplot(~invPow | as.factor(sID),
            data[data$trial2use == 1, ], xlab = xlab)

## Accuracy main-effect model
Mpower_accuracy <- lmer(
  invPow ~ 1 + accuracy + (1 + accuracy | sID),
  data    = data[data$trial2use == 1, ],
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl   = list(maxfun = 196000))
)
summary(Mpower_accuracy)
esticon(Mpower_accuracy, diag(length(fixef(Mpower_accuracy))))

## Optional: full model with accuracy × valence × action
## (non-convergence is acceptable, consistent with the original logic)
Mpower_fullAcc <- lmer(
  invPow ~ 1 + valence * action * accuracy +
    (1 + valence * action * accuracy | sID),
  data    = data[data$trial2use == 1, ],
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl   = list(maxfun = 196000))
)
## Boundary / convergence warnings are acceptable (same as original)

## ------------------ 3.3 Error trials: MFpower ------------------

## Reset to raw MFpower, then mean-correct within each subject
## for error trials with trial2use == 1
data$MFpower <- MFpower$MFpower

for (iSub in sort(unique(data$sID))) {
  tmp <- data$MFpower[data$sID == iSub & data$DVacc == 0 & data$trial2use == 1]
  tmp <- tmp - mean(tmp, na.rm = TRUE)
  data$MFpower[data$sID == iSub & data$DVacc == 0 & data$trial2use == 1] <- tmp
}

dev.new(); xlab <- "MFpower"
densityplot(~MFpower | as.factor(sID),
            data[data$DVacc == 0 & data$trial2use == 1, ], xlab = xlab)

## Inverse transform
data$invPow <- data$MFpower
data$invPow[!is.na(data$MFpower)] <-
  1 - 1 / (data$MFpower[!is.na(data$MFpower)] +
             (1 - min(data$MFpower[!is.na(data$MFpower)])))

dev.new(); xlab <- "invPow"
densityplot(~invPow | as.factor(sID),
            data[data$DVacc == 0 & data$trial2use == 1, ], xlab = xlab)

## Error-trial model: invPow ~ action * valence
Mpower_error <- lmer(
  invPow ~ 1 + action * valence + (1 + action * valence | sID),
  data    = data[data$DVacc == 0 & data$trial2use == 1, ],
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl   = list(maxfun = 196000))
)
summary(Mpower_error)
esticon(Mpower_error, diag(length(fixef(Mpower_error))))

## =====================================================================
## 4. SINGLE-TRIAL MIDFRONTAL-PREFRONTAL ICPC
## =====================================================================

## ------------------ 4.1 PFC ICPC: Valence × Action ------------------

data$ICPCpfc <- MFpower$ICPCpfc

for (iSub in sort(unique(data$sID))) {
  tmp <- data$ICPCpfc[data$sID == iSub & data$DVacc == 1 & data$trial2use == 1]
  tmp <- tmp - mean(tmp, na.rm = TRUE)
  data$ICPCpfc[data$sID == iSub & data$DVacc == 1 & data$trial2use == 1] <- tmp
}

dev.new(); xlab <- "ICPCpfc"
densityplot(~ICPCpfc | as.factor(sID),
            data[data$DVacc == 1 & data$trial2use == 1, ], xlab = xlab)

## Inverse transform
data$invICPCpfc <- data$ICPCpfc
data$invICPCpfc[!is.na(data$ICPCpfc)] <-
  1 / (-data$ICPCpfc[!is.na(data$ICPCpfc)] +
         max(data$ICPCpfc[!is.na(data$ICPCpfc)]) + 1)

for (iSub in sort(unique(data$sID))) {
  tmp <- data$invICPCpfc[data$sID == iSub & data$DVacc == 1 & data$trial2use == 1]
  tmp <- tmp - mean(tmp, na.rm = TRUE)
  data$invICPCpfc[data$sID == iSub & data$DVacc == 1 & data$trial2use == 1] <- tmp
}

dev.new(); xlab <- "invICPCpfc"
densityplot(~invICPCpfc | as.factor(sID),
            data[data$DVacc == 1 & data$trial2use == 1, ], xlab = xlab)

## Mixed model: invICPCpfc ~ action * valence (no intercept)
Micpc_pfc_fit <- lmer(
  invICPCpfc ~ 0 + action * valence + (0 + action * valence | sID),
  data    = data[data$DVacc == 1 & data$trial2use == 1, ],
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl   = list(maxfun = 196000))
)
summary(Micpc_pfc_fit)
esticon(Micpc_pfc_fit, diag(length(fixef(Micpc_pfc_fit))))

## Diagnostics (optional)
densityplot(resid(Micpc_pfc_fit))
qqPlot(resid(Micpc_pfc_fit))
plot(Micpc_pfc_fit, type = c("p", "smooth"))
car::scatterplot(
  fitted(Micpc_pfc_fit),
  data$invICPCpfc[data$DVacc == 1 & data$trial2use == 1],
  boxplots = FALSE, smooth = FALSE
)

## =====================================================================
## 5. SINGLE-TRIAL MIDFRONTAL-MOTOR ICPC
## =====================================================================

## ------------------ 5.1 Build contra / ipsi / nogo ICPC ------------------

data$ICPCnogo <- NaN
data$ICPCnogo[data$action == -1] <-
  (MFpower$l_ICPCmotor[data$action == -1] +
     MFpower$r_ICPCmotor[data$action == -1]) / 2

data$IVleft  <- MFpower$IVleft
data$IVright <- MFpower$IVright

data$ICPCcontra <- NaN
data$ICPCcontra[data$IVleft == 1]  <- MFpower$r_ICPCmotor[MFpower$IVleft == 1]
data$ICPCcontra[data$IVright == 1] <- MFpower$l_ICPCmotor[MFpower$IVright == 1]

data$ICPCipsi <- NaN
data$ICPCipsi[data$IVleft == 1]  <- MFpower$l_ICPCmotor[MFpower$IVleft == 1]
data$ICPCipsi[data$IVright == 1] <- MFpower$r_ICPCmotor[MFpower$IVright == 1]

## ------------------ 5.2 Reshape into long-format icpc ------------------

icpc <- data.frame(
  sID       = rep(data$sID, 3),
  valence   = rep(data$valence, 3),
  trial2use = rep(data$trial2use, 3),
  DVacc     = rep(data$DVacc, 3)
)
icpc$DVicpc <- c(data$ICPCcontra, data$ICPCipsi, data$ICPCnogo)

## Mean-correct within each subject for correct trials with trial2use == 1
for (iSub in sort(unique(icpc$sID))) {
  tmp <- icpc$DVicpc[icpc$sID == iSub & icpc$DVacc == 1 & icpc$trial2use == 1]
  tmp <- tmp - mean(tmp, na.rm = TRUE)
  icpc$DVicpc[icpc$sID == iSub & icpc$DVacc == 1 & icpc$trial2use == 1] <- tmp
}

icpc$execMotor <- c(rep(1, length(data$sID)),
                    rep(2, length(data$sID)),
                    rep(3, length(data$sID)))
icpc$execMotor <- factor(icpc$execMotor,
                         levels = c(1, 2, 3),
                         labels = c("contra", "ipsi", "nogo"))
summary(icpc)

dev.new(); xlab <- "DVicpc"
densityplot(~DVicpc | as.factor(sID),
            icpc[icpc$DVacc == 1 & icpc$trial2use == 1, ], xlab = xlab)
densityplot(~DVicpc | as.factor(execMotor),
            icpc[icpc$DVacc == 1 & icpc$trial2use == 1, ], xlab = xlab)

## Inverse transform
icpc$DVinvicpc <- icpc$DVicpc
icpc$DVinvicpc[!is.na(icpc$DVinvicpc)] <-
  1 / (-icpc$DVinvicpc[!is.na(icpc$DVinvicpc)] +
         max(icpc$DVinvicpc[!is.na(icpc$DVinvicpc)]) + 1)

for (iSub in sort(unique(icpc$sID))) {
  tmp <- icpc$DVinvicpc[icpc$sID == iSub & icpc$DVacc == 1 & icpc$trial2use == 1]
  tmp <- tmp - mean(tmp, na.rm = TRUE)
  icpc$DVinvicpc[icpc$sID == iSub & icpc$DVacc == 1 & icpc$trial2use == 1] <- tmp
}

dev.new(); xlab <- "DVinvicpc"
densityplot(~DVinvicpc | as.factor(sID),
            icpc[icpc$DVacc == 1 & icpc$trial2use == 1, ], xlab = xlab)
densityplot(~DVinvicpc | as.factor(execMotor),
            icpc[icpc$DVacc == 1 & icpc$trial2use == 1, ], xlab = xlab)

## ------------------ 5.3 Motor ICPC model (dummy coding) ------------------

Micpc_motor0 <- lmer(
  DVinvicpc ~ 0 + valence * execMotor +
    (0 + valence * execMotor | sID),
  data    = icpc[icpc$DVacc == 1 & icpc$trial2use == 1, ],
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl   = list(maxfun = 196000))
)

## Switch to dummy coding (consistent with the original Swart script)
icpc$win_contra   <- icpc$win_ipsi   <- icpc$win_nogo   <- rep(0, length(icpc$sID))
icpc$avoid_contra <- icpc$avoid_ipsi <- icpc$avoid_nogo <- rep(0, length(icpc$sID))

icpc$win_contra[(icpc$valence == 1)  & (icpc$execMotor == "contra")] <- 1
icpc$win_ipsi[(icpc$valence == 1)    & (icpc$execMotor == "ipsi")]   <- 1
icpc$win_nogo[(icpc$valence == 1)    & (icpc$execMotor == "nogo")]   <- 1

icpc$avoid_contra[(icpc$valence == -1) & (icpc$execMotor == "contra")] <- 1
icpc$avoid_ipsi[(icpc$valence == -1)   & (icpc$execMotor == "ipsi")]   <- 1
icpc$avoid_nogo[(icpc$valence == -1)   & (icpc$execMotor == "nogo")]   <- 1

Micpc_motor_fit <- lmer(
  DVicpc ~ 0 + win_contra + win_ipsi + win_nogo +
    avoid_contra + avoid_ipsi + avoid_nogo +
    (0 + win_contra + win_ipsi + win_nogo +
       avoid_contra + avoid_ipsi + avoid_nogo | sID),
  data    = icpc[icpc$DVacc == 1 & icpc$trial2use == 1, ],
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl   = list(maxfun = 196000))
)
summary(Micpc_motor_fit)
esticon(Micpc_motor_fit, diag(length(fixef(Micpc_motor_fit))))

## =====================================================================
## 6. Correlations between MFpower and ICPC (compact version)
## =====================================================================

data$MFpower  <- MFpower$MFpower
data$ICPCpfc  <- MFpower$ICPCpfc
data$ICPCnogo <- NaN
data$ICPCnogo[data$action == -1] <-
  (MFpower$l_ICPCmotor[data$action == -1] +
     MFpower$r_ICPCmotor[data$action == -1]) / 2

data$ICPCcontra <- NaN
data$ICPCcontra[data$IVleft == 1]  <- MFpower$r_ICPCmotor[MFpower$IVleft == 1]
data$ICPCcontra[data$IVright == 1] <- MFpower$l_ICPCmotor[MFpower$IVright == 1]

data$ICPCipsi <- NaN
data$ICPCipsi[data$IVleft == 1]  <- MFpower$l_ICPCmotor[MFpower$IVleft == 1]
data$ICPCipsi[data$IVright == 1] <- MFpower$r_ICPCmotor[MFpower$IVright == 1]

subj_ids <- sort(unique(data$sID))
nSub     <- length(subj_ids)

corrPFC   <- matrix(NA, nSub, 4)
corrMotor <- matrix(NA, nSub, 6)

for (ii in seq_along(subj_ids)) {
  iSub <- subj_ids[ii]
  
  ## PFC
  idx <- data$sID == iSub & data$trial2use == 1 & data$valence == 1 & data$action == 1
  corrPFC[ii, 1] <- suppressWarnings(
    cor.test(data$MFpower[idx], data$ICPCpfc[idx], method = "spearman")$estimate
  )
  idx <- data$sID == iSub & data$trial2use == 1 & data$valence == -1 & data$action == 1
  corrPFC[ii, 2] <- suppressWarnings(
    cor.test(data$MFpower[idx], data$ICPCpfc[idx], method = "spearman")$estimate
  )
  idx <- data$sID == iSub & data$trial2use == 1 & data$valence == 1 & data$action == -1
  corrPFC[ii, 3] <- suppressWarnings(
    cor.test(data$MFpower[idx], data$ICPCpfc[idx], method = "spearman")$estimate
  )
  idx <- data$sID == iSub & data$trial2use == 1 & data$valence == -1 & data$action == -1
  corrPFC[ii, 4] <- suppressWarnings(
    cor.test(data$MFpower[idx], data$ICPCpfc[idx], method = "spearman")$estimate
  )
  
  ## Motor
  idx <- data$sID == iSub & data$trial2use == 1 & data$valence == 1 & data$action == 1
  corrMotor[ii, 1] <- suppressWarnings(
    cor.test(data$MFpower[idx], data$ICPCcontra[idx], method = "spearman")$estimate
  )
  idx <- data$sID == iSub & data$trial2use == 1 & data$valence == -1 & data$action == 1
  corrMotor[ii, 2] <- suppressWarnings(
    cor.test(data$MFpower[idx], data$ICPCcontra[idx], method = "spearman")$estimate
  )
  idx <- data$sID == iSub & data$trial2use == 1 & data$valence == 1 & data$action == 1
  corrMotor[ii, 3] <- suppressWarnings(
    cor.test(data$MFpower[idx], data$ICPCipsi[idx], method = "spearman")$estimate
  )
  idx <- data$sID == iSub & data$trial2use == 1 & data$valence == -1 & data$action == 1
  corrMotor[ii, 4] <- suppressWarnings(
    cor.test(data$MFpower[idx], data$ICPCipsi[idx], method = "spearman")$estimate
  )
  idx <- data$sID == iSub & data$trial2use == 1 & data$valence == 1 & data$action == -1
  corrMotor[ii, 5] <- suppressWarnings(
    cor.test(data$MFpower[idx], data$ICPCnogo[idx], method = "spearman")$estimate
  )
  idx <- data$sID == iSub & data$trial2use == 1 & data$valence == -1 & data$action == -1
  corrMotor[ii, 6] <- suppressWarnings(
    cor.test(data$MFpower[idx], data$ICPCnogo[idx], method = "spearman")$estimate
  )
}

## Summary of correlation means / ranges (optional)
mean(apply(corrPFC,   2, function(x) mean(x, na.rm = TRUE)))
min(apply(corrPFC,    2, function(x) min(x,  na.rm = TRUE)))
max(apply(corrPFC,    2, function(x) max(x,  na.rm = TRUE)))

mean(apply(corrMotor, 2, function(x) mean(x, na.rm = TRUE)))
min(apply(corrMotor,  2, function(x) min(x,  na.rm = TRUE)))
max(apply(corrMotor,  2, function(x) max(x,  na.rm = TRUE)))

## Export to CSV (each row corresponds to one subject, ordered by subj_ids)
write.table(as.table(corrPFC),
            "H:/mt_behavalluse/EEGmixmodel/corrPFC_power_isps_OCD.csv",
            sep = ",", col.names = FALSE, row.names = FALSE)
write.table(as.table(corrMotor),
            "H:/mt_behavalluse/EEGmixmodel/corrMotor_power_isps_OCD.csv",
            sep = ",", col.names = FALSE, row.names = FALSE)

## =====================================================================
## 7. Extract fixed effects + simple effects and export to Excel
## =====================================================================

extract_esticon_df <- function(fit, C, effect_names,
                               model_name, effect_type) {
  tab <- esticon(fit, C)
  tab <- as.data.frame(tab)
  tab$model      <- model_name
  tab$effectType <- effect_type
  tab$effect     <- effect_names
  tab <- tab[, c("model", "effectType", "effect",
                 "estimate", "std.error", "statistic", "df", "p.value")]
  rownames(tab) <- NULL
  tab
}

## 7.1 MF power (correct trials): fixed effects
stopifnot(identical(
  names(fixef(Mpower_correct)),
  c("(Intercept)", "action", "valence", "action:valence")
))

C_power_fix <- diag(length(fixef(Mpower_correct)))
rownames(C_power_fix) <- names(fixef(Mpower_correct))

res_power_fix <- extract_esticon_df(
  fit          = Mpower_correct,
  C            = C_power_fix,
  effect_names = rownames(C_power_fix),
  model_name   = "MFpower_correct",
  effect_type  = "fixed"
)

## 7.2 MF power: simple effects
C_power_simple <- rbind(
  valence_at_Go   = c(0, 0, 1,  1),
  valence_at_NoGo = c(0, 0, 1, -1),
  action_at_Win   = c(0, 1, 0,  1),
  action_at_Avoid = c(0, 1, 0, -1)
)

res_power_simple <- extract_esticon_df(
  fit          = Mpower_correct,
  C            = C_power_simple,
  effect_names = rownames(C_power_simple),
  model_name   = "MFpower_correct",
  effect_type  = "simple"
)

## 7.3 MF power + accuracy model
C_acc_fix <- diag(length(fixef(Mpower_accuracy)))
rownames(C_acc_fix) <- names(fixef(Mpower_accuracy))

res_acc_fix <- extract_esticon_df(
  fit          = Mpower_accuracy,
  C            = C_acc_fix,
  effect_names = rownames(C_acc_fix),
  model_name   = "MFpower_accuracy",
  effect_type  = "fixed"
)

## 7.4 MF power (error trials)
C_error_fix <- diag(length(fixef(Mpower_error)))
rownames(C_error_fix) <- names(fixef(Mpower_error))

res_error_fix <- extract_esticon_df(
  fit          = Mpower_error,
  C            = C_error_fix,
  effect_names = rownames(C_error_fix),
  model_name   = "MFpower_error",
  effect_type  = "fixed"
)

## 7.5 PFC ICPC: fixed effects
stopifnot(identical(
  names(fixef(Micpc_pfc_fit)),
  c("action", "valence", "action:valence")
))

C_pfc_fix <- diag(length(fixef(Micpc_pfc_fit)))
rownames(C_pfc_fix) <- names(fixef(Micpc_pfc_fit))

res_pfc_fix <- extract_esticon_df(
  fit          = Micpc_pfc_fit,
  C            = C_pfc_fix,
  effect_names = rownames(C_pfc_fix),
  model_name   = "ICPC_pfc_correct",
  effect_type  = "fixed"
)

## 7.6 PFC ICPC: simple effects
C_pfc_simple <- rbind(
  valence_at_Go   = c(0, 1,  1),
  valence_at_NoGo = c(0, 1, -1),
  action_at_Win   = c(1, 0,  1),
  action_at_Avoid = c(1, 0, -1)
)

res_pfc_simple <- extract_esticon_df(
  fit          = Micpc_pfc_fit,
  C            = C_pfc_simple,
  effect_names = rownames(C_pfc_simple),
  model_name   = "ICPC_pfc_correct",
  effect_type  = "simple"
)

## 7.7 Motor ICPC: fixed effects
stopifnot(identical(
  names(fixef(Micpc_motor_fit)),
  c("win_contra", "win_ipsi", "win_nogo",
    "avoid_contra", "avoid_ipsi", "avoid_nogo")
))

C_motor_fix <- diag(length(fixef(Micpc_motor_fit)))
rownames(C_motor_fix) <- names(fixef(Micpc_motor_fit))

res_motor_fix <- extract_esticon_df(
  fit          = Micpc_motor_fit,
  C            = C_motor_fix,
  effect_names = rownames(C_motor_fix),
  model_name   = "ICPC_motor_correct",
  effect_type  = "fixed"
)

## 7.8 Motor ICPC: simple effects / contrasts
C_motor_simple <- rbind(
  congruency_x_exec = c(1, -1,  1, -1,  1, -1),
  congruency_exec   = c(1,  0,  0, -1,  0,  0),
  congruency_nonExec= c(0, -1,  1,  0,  1, -1),
  valence_main      = c(1,  1,  1, -1, -1, -1),
  contra_vs_nogo    = c(1,  0, -1,  1,  0, -1),
  ipsi_vs_nogo      = c(0,  1, -1,  0,  1, -1)
)

res_motor_simple <- extract_esticon_df(
  fit          = Micpc_motor_fit,
  C            = C_motor_simple,
  effect_names = rownames(C_motor_simple),
  model_name   = "ICPC_motor_correct",
  effect_type  = "simple"
)

## 7.9 Combine and write to Excel
out_list <- list(
  MFpower_fixed        = res_power_fix,
  MFpower_simple       = res_power_simple,
  MFpower_accuracy_fix = res_acc_fix,
  MFpower_error_fix    = res_error_fix,
  PFC_fixed            = res_pfc_fix,
  PFC_simple           = res_pfc_simple,
  Motor_fixed          = res_motor_fix,
  Motor_simple         = res_motor_simple
)

write_xlsx(
  out_list,
  path = "H:/mt_behavalluse/EEGmixmodel/EEG_single_trial_OCD_effects.xlsx"
)

cat("Export completed: H:/mt_behavalluse/EEGmixmodel/EEG_single_trial_OCD_effects.xlsx\n")