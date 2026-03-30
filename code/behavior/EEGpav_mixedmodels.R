# ============================================================================
# MIXED MODEL ANALYSIS OF MOTIVATIONAL GO/NOGO DATA OF THE EEG STUDY.
# 1. Basic analyses.
# By Jennifer C. Swart, 2015.
# Modified: run models separately for HC (group=0) and OCD (group=1)
# ============================================================================

# SETUP.
# ----------------------------------------------------------------------------
setwd("H:/mt_behavalluse")

# install packages if necessary.
#install.packages("lme4") 
#install.packages("car")
#install.packages("pbkrtest") 
#install.packages("lattice") 
#install.packages("doBy") 
#install.packages("ggplot2") 
#install.packages("psych") 

# load package for mixed models.
library(lme4)      # for (g)lmer.
library(car)       # for Anova.
library(pbkrtest)
library(lattice)   # for plotting.
library(doBy)      # for esticon.
library(ggplot2)   # for qplot.
library(psych)     # for describeBy.
source("simpef.R")

###############################################################################
# PART 1: HC group (group == 0)
###############################################################################
cat("============== Running models for HC group (group == 0) ==============\n")

# load the data.
data <- read.csv("EEG4mixedmodelR_withGroup.csv", header = TRUE)

# keep HC group only
data <- data[data$group == 0, ]

# recode IVs to [-1 1]:
data$valence[data$valence == 0] = -1
data$action[data$action == 0]   = -1
data$accuracy = data$DVacc
data$accuracy[data$accuracy == 0] = -1

# prepare DVs.
data$DVcorrectGo   = data$DVgo * data$DVacc
data$DVincorrectGo = data$DVgo * (1 - data$DVacc)

# inspect the smallest RTs (first 20)
sort(data$DVrt)[1:20]

# discard RTs < 100 ms.
data$DVrt[data$DVrt < .1] = NaN

# ln-transform to improve normality: make sure to add a constant such that RT > 1.
data$DVlnrt = log(data$DVrt + .9)

# perform some data checks.
summary(data)
head(data)
tail(data)

filename = "EEGpav_mm_HC.RData"

# 1. PERFORMING MIXED MODEL ANALYSES PER DATA TYPE.
# -----------------------------------------------------------------------------

# P(GO).
# ============================================================================

# MIXED-EFFECTS MODEL.
Mgo.basic = glmer(DVgo ~ valence*action + (1+valence*action|sID), family = binomial, data = data, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1960)))
summary(Mgo.basic)
save(data, Mgo.basic, file = filename)

# stats.
names(fixef(Mgo.basic))
esticon(Mgo.basic, diag(length(fixef(Mgo.basic))))
simpef(Mgo.basic, "valence:action", "valence", "action")

# Correct & Incorrect Go responses (Go cues only).
# ============================================================================

# Correct Go
Mcorrectgo = glmer(DVcorrectGo ~ valence + (1+valence|sID), family = binomial, data = data[data$DVincorrectGo != 1 & data$action == 1, ], control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 250)))
summary(Mcorrectgo)
save(data, Mgo.basic, Mcorrectgo, file = filename)
names(fixef(Mcorrectgo))
esticon(Mcorrectgo, diag(length(fixef(Mcorrectgo))))

# Incorrect Go
Mincorrectgo = glmer(DVincorrectGo ~ valence + (1+valence|sID), family = binomial, data = data[data$DVcorrectGo != 1 & data$action == 1, ], control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 250)))
summary(Mincorrectgo)
save(data, Mgo.basic, Mcorrectgo, Mincorrectgo, file = filename)
names(fixef(Mincorrectgo))
esticon(Mincorrectgo, diag(length(fixef(Mincorrectgo))))

# p(correct|Go)
Maccuracyofgo = glmer(DVacc ~ valence + (1+valence|sID), family = binomial, data = data[data$DVgo == 1 & data$action == 1, ], control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 250)))
summary(Maccuracyofgo)
save(data, Mgo.basic, Mcorrectgo, Mincorrectgo, Maccuracyofgo, file = filename)
names(fixef(Maccuracyofgo))
esticon(Maccuracyofgo, diag(length(fixef(Maccuracyofgo))))

# RESPONSE TIMES.
# ============================================================================

# DIAGNOSTICS RAW DATA.
xlab = "RT (raw)"
densityplot(~DVrt, data, xlab = xlab)
densityplot(~DVrt | as.factor(sID), data, xlab = xlab)
densityplot(~DVrt | as.factor(action), data, xlab = xlab, strip = strip.custom(factor.levels = c("NoGo", "Go")))
densityplot(~DVrt | as.factor(valence), data, xlab = xlab, strip = strip.custom(factor.levels = c("Aversive", "Appetitive")))

xlab = "RT (ln-transformed)"
densityplot(~DVlnrt, data, xlab = xlab)
densityplot(~DVlnrt | as.factor(sID), data, xlab = xlab)
densityplot(~DVlnrt | as.factor(action), data, xlab = xlab, strip = strip.custom(factor.levels = c("NoGo", "Go")))
densityplot(~DVlnrt | as.factor(valence), data, xlab = xlab, strip = strip.custom(factor.levels = c("Aversive", "Appetitive")))

# MIXED-EFFECTS MODEL.
Mrt.basic = lmer(DVlnrt ~ valence*action + (1+valence*action|sID), data = data, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1960)))
summary(Mrt.basic)
save(data, Mgo.basic, Mcorrectgo, Mincorrectgo, Maccuracyofgo, Mrt.basic, file = filename)
names(fixef(Mrt.basic))
esticon(Mrt.basic, diag(length(fixef(Mrt.basic))))

# MODEL DIAGNOSTICS: BASIC MODEL.
densityplot(resid(Mrt.basic))
qqPlot(resid(Mrt.basic))
plot(Mrt.basic, sID[!is.na(DVlnrt)] ~ resid(.))
plot(Mrt.basic, type = c("p", "smooth"))
scatterplot(fitted(Mrt.basic), data$DVlnrt[!is.na(data$DVlnrt)], boxplot = FALSE, smoother = FALSE)

# PRELIM STATISTICS.
Anova(Mrt.basic)

# RT valence effect different for correct vs incorrect RT? (Go cues only)
Mrt.accuracy = lmer(DVlnrt ~ valence*accuracy + (1+valence*accuracy|sID), data = data[data$action == 1, ], control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1960)))
summary(Mrt.accuracy)
Anova(Mrt.accuracy)
save(data, Mgo.basic, Mcorrectgo, Mincorrectgo, Maccuracyofgo, Mrt.basic, Mrt.accuracy, file = filename)
names(fixef(Mrt.accuracy))
esticon(Mrt.accuracy, diag(length(fixef(Mrt.basic))))

# RELATION VALENCE EFFECT P(GO) & RT.
ranefGo = ranef(Mgo.basic)$sID
ranefRT = ranef(Mrt.basic)$sID
cor.test(x = ranefGo$valence, y = ranefRT$valence, method = "pearson")
cor.test(x = ranefGo$valence, y = ranefRT$valence, method = "spearman")
plot(x = ranefGo$valence, y = ranefRT$valence, xlab = "p(Go)", ylab = "RT")
abline(lm(ranefRT$valence ~ ranefGo$valence))
dev.new(width = 2.3, height = 2.3)
qplot(x = ranefGo$valence, y = ranefRT$valence, xlab = "p(Go)", ylab = "RT", size = I(1.2)) + stat_smooth(method = "lm", colour = "black", alpha = 1) + geom_point(size = I(1.2)) + theme(panel.grid = element_blank()) + theme_classic()
ggsave(filename = "ranefpGo_RT_HC.eps")
dev.off()



## ==========================================================
## Extract beta / statistics / chi-square / p values
## and export them to Excel
## Save results to: H:/mt_behavalluse/behavmixmodel
## ==========================================================

library(writexl)

get_effect_table <- function(m){
  sm    <- summary(m)
  coefs <- sm$coefficients
  
  # support both glmer (z) and lmer (t)
  if ("z value" %in% colnames(coefs)) {
    stat_name <- "z value"
    p_name    <- "Pr(>|z|)"
  } else if ("t value" %in% colnames(coefs)) {
    stat_name <- "t value"
    p_name    <- "Pr(>|t|)"
  } else {
    stop("Cannot find z value or t value columns. Please check the model summary.")
  }
  
  beta <- coefs[, "Estimate"]
  stat <- coefs[, stat_name]
  p    <- coefs[, p_name]
  chi2 <- stat^2
  
  out <- data.frame(
    term = rownames(coefs),
    beta = beta,
    stat = stat,
    chi2 = chi2,
    p    = p,
    row.names = NULL
  )
  
  # intercept is usually not of interest; comment out the next line if needed
  out <- out[out$term != "(Intercept)", ]
  
  return(out)
}

## ---------- 1. Extract result tables from each model ----------

tab_Mgo_HC        <- get_effect_table(Mgo.basic)
tab_Mcorrectgo    <- get_effect_table(Mcorrectgo)
tab_Mincorrectgo  <- get_effect_table(Mincorrectgo)
tab_Maccuracyofgo <- get_effect_table(Maccuracyofgo)
tab_Mrt_basic     <- get_effect_table(Mrt.basic)
tab_Mrt_accuracy  <- get_effect_table(Mrt.accuracy)

## ---------- 2. Organize into a list, one element per sheet ----------

results_list <- list(
  HC_Mgo_basic      = tab_Mgo_HC,
  HC_McorrectGo     = tab_Mcorrectgo,
  HC_MincorrectGo   = tab_Mincorrectgo,
  HC_MaccuracyOfGo  = tab_Maccuracyofgo,
  HC_Mrt_basic      = tab_Mrt_basic,
  HC_Mrt_accuracy   = tab_Mrt_accuracy
)

## ---------- 3. Define output directory and file name ----------

out_dir  <- "H:/mt_behavalluse/behavmixmodel"

# create directory if it does not exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_file <- file.path(out_dir, "EEG_HC_mixedmodel_effects.xlsx")

## ---------- 4. Write to Excel ----------

write_xlsx(results_list, path = out_file)

cat("All model results have been written to: ", out_file, "\n")









###############################################################################
# PART 2: OCD group (group == 1)
###############################################################################
cat("============== Running models for OCD group (group == 1) =============\n")

# reload the full dataset
data <- read.csv("EEG4mixedmodelR_withGroup.csv", header = TRUE)

# keep OCD group only
data <- data[data$group == 1, ]

# recode IVs to [-1 1]:
data$valence[data$valence == 0] = -1
data$action[data$action == 0]   = -1
data$accuracy = data$DVacc
data$accuracy[data$accuracy == 0] = -1

# prepare DVs.
data$DVcorrectGo   = data$DVgo * data$DVacc
data$DVincorrectGo = data$DVgo * (1 - data$DVacc)

sort(data$DVrt)[1:20]
data$DVrt[data$DVrt < .1] = NaN
data$DVlnrt = log(data$DVrt + .9)

summary(data)
head(data)
tail(data)

filename = "EEGpav_mm_OCD.RData"

# P(GO) model
Mgo.basic = glmer(DVgo ~ valence*action + (1+valence*action|sID), family = binomial, data = data, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1960)))
summary(Mgo.basic)
save(data, Mgo.basic, file = filename)
names(fixef(Mgo.basic))
esticon(Mgo.basic, diag(length(fixef(Mgo.basic))))
simpef(Mgo.basic, "valence:action", "valence", "action")

# Correct Go
Mcorrectgo = glmer(DVcorrectGo ~ valence + (1+valence|sID), family = binomial, data = data[data$DVincorrectGo != 1 & data$action == 1, ], control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 250)))
summary(Mcorrectgo)
save(data, Mgo.basic, Mcorrectgo, file = filename)
names(fixef(Mcorrectgo))
esticon(Mcorrectgo, diag(length(fixef(Mcorrectgo))))

# Incorrect Go
Mincorrectgo = glmer(DVincorrectGo ~ valence + (1+valence|sID), family = binomial, data = data[data$DVcorrectGo != 1 & data$action == 1, ], control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 250)))
summary(Mincorrectgo)
save(data, Mgo.basic, Mcorrectgo, Mincorrectgo, file = filename)
names(fixef(Mincorrectgo))
esticon(Mincorrectgo, diag(length(fixef(Mincorrectgo))))

# p(correct|Go)
Maccuracyofgo = glmer(DVacc ~ valence + (1+valence|sID), family = binomial, data = data[data$DVgo == 1 & data$action == 1, ], control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 250)))
summary(Maccuracyofgo)
save(data, Mgo.basic, Mcorrectgo, Mincorrectgo, Maccuracyofgo, file = filename)
names(fixef(Maccuracyofgo))
esticon(Maccuracyofgo, diag(length(fixef(Maccuracyofgo))))

# RT section
xlab = "RT (raw)"
densityplot(~DVrt, data, xlab = xlab)
densityplot(~DVrt | as.factor(sID), data, xlab = xlab)
densityplot(~DVrt | as.factor(action), data, xlab = xlab, strip = strip.custom(factor.levels = c("NoGo", "Go")))
densityplot(~DVrt | as.factor(valence), data, xlab = xlab, strip = strip.custom(factor.levels = c("Aversive", "Appetitive")))

xlab = "RT (ln-transformed)"
densityplot(~DVlnrt, data, xlab = xlab)
densityplot(~DVlnrt | as.factor(sID), data, xlab = xlab)
densityplot(~DVlnrt | as.factor(action), data, xlab = xlab, strip = strip.custom(factor.levels = c("NoGo", "Go")))
densityplot(~DVlnrt | as.factor(valence), data, xlab = xlab, strip = strip.custom(factor.levels = c("Aversive", "Appetitive")))

Mrt.basic = lmer(DVlnrt ~ valence*action + (1+valence*action|sID), data = data, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1960)))
summary(Mrt.basic)
save(data, Mgo.basic, Mcorrectgo, Mincorrectgo, Maccuracyofgo, Mrt.basic, file = filename)
names(fixef(Mrt.basic))
esticon(Mrt.basic, diag(length(fixef(Mrt.basic))))

densityplot(resid(Mrt.basic))
qqPlot(resid(Mrt.basic))
plot(Mrt.basic, sID[!is.na(DVlnrt)] ~ resid(.))
plot(Mrt.basic, type = c("p", "smooth"))
scatterplot(fitted(Mrt.basic), data$DVlnrt[!is.na(data$DVlnrt)], boxplot = FALSE, smoother = FALSE)

Anova(Mrt.basic)

Mrt.accuracy = lmer(DVlnrt ~ valence*accuracy + (1+valence*accuracy|sID), data = data[data$action == 1, ], control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1960)))
summary(Mrt.accuracy)
Anova(Mrt.accuracy)
save(data, Mgo.basic, Mcorrectgo, Mincorrectgo, Maccuracyofgo, Mrt.basic, Mrt.accuracy, file = filename)
names(fixef(Mrt.accuracy))
esticon(Mrt.accuracy, diag(length(fixef(Mrt.basic))))

ranefGo = ranef(Mgo.basic)$sID
ranefRT = ranef(Mrt.basic)$sID
cor.test(x = ranefGo$valence, y = ranefRT$valence, method = "pearson")
cor.test(x = ranefGo$valence, y = ranefRT$valence, method = "spearman")
plot(x = ranefGo$valence, y = ranefRT$valence, xlab = "p(Go)", ylab = "RT")
abline(lm(ranefRT$valence ~ ranefGo$valence))
dev.new(width = 2.3, height = 2.3)
qplot(x = ranefGo$valence, y = ranefRT$valence, xlab = "p(Go)", ylab = "RT", size = I(1.2)) + stat_smooth(method = "lm", colour = "black", alpha = 1) + geom_point(size = I(1.2)) + theme(panel.grid = element_blank()) + theme_classic()
ggsave(filename = "ranefpGo_RT_OCD.eps")
dev.off()








## ==========================================================
## Extract beta / statistics / chi-square / p values
## and export them to Excel
## Save results to: H:/mt_behavalluse/behavmixmodel
## ==========================================================

library(writexl)

get_effect_table <- function(m){
  sm    <- summary(m)
  coefs <- sm$coefficients
  
  # support both glmer (z) and lmer (t)
  if ("z value" %in% colnames(coefs)) {
    stat_name <- "z value"
    p_name    <- "Pr(>|z|)"
  } else if ("t value" %in% colnames(coefs)) {
    stat_name <- "t value"
    p_name    <- "Pr(>|t|)"
  } else {
    stop("Cannot find z value or t value columns. Please check the model summary.")
  }
  
  beta <- coefs[, "Estimate"]
  stat <- coefs[, stat_name]
  p    <- coefs[, p_name]
  chi2 <- stat^2
  
  out <- data.frame(
    term = rownames(coefs),
    beta = beta,
    stat = stat,
    chi2 = chi2,
    p    = p,
    row.names = NULL
  )
  
  # intercept is usually not of interest; comment out the next line if needed
  out <- out[out$term != "(Intercept)", ]
  
  return(out)
}

## ---------- 1. Extract result tables from each model ----------

tab_Mgo_HC        <- get_effect_table(Mgo.basic)
tab_Mcorrectgo    <- get_effect_table(Mcorrectgo)
tab_Mincorrectgo  <- get_effect_table(Mincorrectgo)
tab_Maccuracyofgo <- get_effect_table(Maccuracyofgo)
tab_Mrt_basic     <- get_effect_table(Mrt.basic)
tab_Mrt_accuracy  <- get_effect_table(Mrt.accuracy)

## ---------- 2. Organize into a list, one element per sheet ----------

results_list <- list(
  HC_Mgo_basic      = tab_Mgo_HC,
  HC_McorrectGo     = tab_Mcorrectgo,
  HC_MincorrectGo   = tab_Mincorrectgo,
  HC_MaccuracyOfGo  = tab_Maccuracyofgo,
  HC_Mrt_basic      = tab_Mrt_basic,
  HC_Mrt_accuracy   = tab_Mrt_accuracy
)

## ---------- 3. Define output directory and file name ----------

out_dir  <- "H:/mt_behavalluse/behavmixmodel"

# create directory if it does not exist
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

out_file <- file.path(out_dir, "EEG_OCD_mixedmodel_effects.xlsx")

## ---------- 4. Write to Excel ----------

write_xlsx(results_list, path = out_file)

cat("All model results have been written to: ", out_file, "\n")










# ============================================================================
# MIXED MODEL ANALYSIS OF MOTIVATIONAL GO/NOGO DATA OF THE EEG STUDY.
# Joint model with group (HC vs OCD) as between-subject factor.
# By Jennifer C. Swart, 2015. (modified: joint-group version)
# ============================================================================

# WORKING DIRECTORY & DATA ----------------------------------------------------
setwd("H:/mt_behavalluse")
data <- read.csv("EEG4mixedmodelR_withGroup.csv", header = TRUE)

# SETUP -----------------------------------------------------------------------
# install.packages("lme4")
# install.packages("car")
# install.packages("pbkrtest")
# install.packages("lattice")
# install.packages("doBy")
# install.packages("ggplot2")
# install.packages("psych")

library(lme4)      # for glmer / lmer
library(car)       # for Anova
library(pbkrtest)
library(lattice)   # for plotting
library(doBy)      # for esticon
library(ggplot2)   # for qplot
library(psych)     # for describeBy
source("simpef.R")

# RECODE IVs & PREPARE DVs ----------------------------------------------------
# valence, action [-1, 1]
data$valence[data$valence == 0] <- -1
data$action[data$action == 0]   <- -1

# group: 0 = HC, 1 = OCD -> factor
data$group <- factor(data$group, levels = c(0, 1), labels = c("HC", "OCD"))

# accuracy: -1 / 1
data$accuracy <- data$DVacc
data$accuracy[data$accuracy == 0] <- -1

# DVs
data$DVcorrectGo   <- data$DVgo * data$DVacc
data$DVincorrectGo <- data$DVgo * (1 - data$DVacc)

# inspect the smallest RTs (first 20)
sort(data$DVrt)[1:20]

# discard RT < 100 ms
data$DVrt[data$DVrt < .1] <- NaN

# ln-transform (make sure RT > 1)
data$DVlnrt <- log(data$DVrt + .9)

# data checks
summary(data)
head(data)
tail(data)

# output file name (joint model)
filename <- "EEGpav_group.RData"

# 1. P(GO) model ---------------------------------------------------------------
# DV: DVgo (Go vs NoGo)
# Fixed effects: group * valence * action
# Random effects: subject intercept + random slopes for valence, action, valence:action

Mgo.joint <- glmer(DVgo ~ group*valence*action + (1+valence*action|sID), family = binomial, data = data, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1960)))
summary(Mgo.joint)
save(data, Mgo.joint, file = filename)

# stats
names(fixef(Mgo.joint))
esticon(Mgo.joint, diag(length(fixef(Mgo.joint))))
# original simple-effects call kept as valence:action only, without group
simpef(Mgo.joint, "valence:action", "valence", "action")

# 2. Correct & Incorrect Go responses (Go cues only) -------------------------
# group is included in the fixed effects

# 2.1 Correct Go: DVcorrectGo ~ group*valence
Mcorrectgo.joint <- glmer(DVcorrectGo ~ group*valence + (1+valence|sID), family = binomial, data = data[data$DVincorrectGo != 1 & data$action == 1, ], control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 250)))
summary(Mcorrectgo.joint)
save(data, Mgo.joint, Mcorrectgo.joint, file = filename)
names(fixef(Mcorrectgo.joint))
esticon(Mcorrectgo.joint, diag(length(fixef(Mcorrectgo.joint))))

# 2.2 Incorrect Go: DVincorrectGo ~ group*valence
Mincorrectgo.joint <- glmer(DVincorrectGo ~ group*valence + (1+valence|sID), family = binomial, data = data[data$DVcorrectGo != 1 & data$action == 1, ], control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 250)))
summary(Mincorrectgo.joint)
save(data, Mgo.joint, Mcorrectgo.joint, Mincorrectgo.joint, file = filename)
names(fixef(Mincorrectgo.joint))
esticon(Mincorrectgo.joint, diag(length(fixef(Mincorrectgo.joint))))

# 2.3 p(correct | Go): DVacc ~ group*valence
Maccuracyofgo.joint <- glmer(DVacc ~ group*valence + (1+valence|sID), family = binomial, data = data[data$DVgo == 1 & data$action == 1, ], control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 250)))
summary(Maccuracyofgo.joint)
save(data, Mgo.joint, Mcorrectgo.joint, Mincorrectgo.joint, Maccuracyofgo.joint, file = filename)
names(fixef(Maccuracyofgo.joint))
esticon(Maccuracyofgo.joint, diag(length(fixef(Maccuracyofgo.joint))))

# 3. RESPONSE TIMES -----------------------------------------------------------

# 3.1 Diagnostics: raw data
xlab <- "RT (raw)"
densityplot(~DVrt, data, xlab = xlab)
densityplot(~DVrt | as.factor(sID), data, xlab = xlab)
densityplot(~DVrt | as.factor(action), data, xlab = xlab, strip = strip.custom(factor.levels = c("NoGo", "Go")))
densityplot(~DVrt | as.factor(valence), data, xlab = xlab, strip = strip.custom(factor.levels = c("Aversive", "Appetitive")))

xlab <- "RT (ln-transformed)"
densityplot(~DVlnrt, data, xlab = xlab)
densityplot(~DVlnrt | as.factor(sID), data, xlab = xlab)
densityplot(~DVlnrt | as.factor(action), data, xlab = xlab, strip = strip.custom(factor.levels = c("NoGo", "Go")))
densityplot(~DVlnrt | as.factor(valence), data, xlab = xlab, strip = strip.custom(factor.levels = c("Aversive", "Appetitive")))

# 3.2 Mixed-effects model: RT basic model with group
# DV: DVlnrt
# Fixed effects: group * valence * action
# Random effects: 1 + valence*action | sID

Mrt.basic.joint <- lmer(DVlnrt ~ group*valence*action + (1+valence*action|sID), data = data, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1960)))
summary(Mrt.basic.joint)
save(data, Mgo.joint, Mcorrectgo.joint, Mincorrectgo.joint, Maccuracyofgo.joint, Mrt.basic.joint, file = filename)
names(fixef(Mrt.basic.joint))
esticon(Mrt.basic.joint, diag(length(fixef(Mrt.basic.joint))))

# 3.3 Model diagnostics: basic model
densityplot(resid(Mrt.basic.joint))
qqPlot(resid(Mrt.basic.joint))
plot(Mrt.basic.joint, sID[!is.na(DVlnrt)] ~ resid(.))
plot(Mrt.basic.joint, type = c("p", "smooth"))
scatterplot(fitted(Mrt.basic.joint), data$DVlnrt[!is.na(data$DVlnrt)], boxplot = FALSE, smoother = FALSE)

# 3.4 Preliminary statistics
Anova(Mrt.basic.joint)

# 3.5 RT valence effect different for correct vs incorrect RT? (Go cues only)
# model includes group*valence*accuracy
Mrt.accuracy.joint <- lmer(DVlnrt ~ group*valence*accuracy + (1+valence*accuracy|sID), data = data[data$action == 1, ], control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1960)))
summary(Mrt.accuracy.joint)
Anova(Mrt.accuracy.joint)
save(data, Mgo.joint, Mcorrectgo.joint, Mincorrectgo.joint, Maccuracyofgo.joint, Mrt.basic.joint, Mrt.accuracy.joint, file = filename)
names(fixef(Mrt.accuracy.joint))
esticon(Mrt.accuracy.joint, diag(length(fixef(Mrt.basic.joint))))

# 4. RELATION BETWEEN VALENCE EFFECTS IN P(GO) AND RT -------------------------
# random-effect valence slopes from the joint model (group already controlled)

ranefGo  <- ranef(Mgo.joint)$sID
ranefRT  <- ranef(Mrt.basic.joint)$sID
cor.test(x = ranefGo$valence, y = ranefRT$valence, method = "pearson")
cor.test(x = ranefGo$valence, y = ranefRT$valence, method = "spearman")

plot(x = ranefGo$valence, y = ranefRT$valence, xlab = "p(Go)", ylab = "RT")
abline(lm(ranefRT$valence ~ ranefGo$valence))

dev.new(width = 2.3, height = 2.3)
qplot(x = ranefGo$valence, y = ranefRT$valence, xlab = "p(Go)", ylab = "RT", size = I(1.2)) +
  stat_smooth(method = "lm", colour = "black", alpha = 1) +
  geom_point(size = I(1.2)) +
  theme(panel.grid = element_blank()) +
  theme_classic()
ggsave(filename = "ranefpGo_RT_group.eps")
dev.off()



















## ====================================================================
## Extract all model results from EEGpav_mm_OCD.RData
## and export them as summary tables
## Path: H:/mt_behavalluse
## Included models:
##   Mgo.basic, Mcorrectgo, Mincorrectgo, Maccuracyofgo, Mrt.basic, Mrt.accuracy
## ====================================================================

## 1) Working directory & load objects ----------------------------------
setwd("H:/mt_behavalluse")

# load RData (contains data and saved model objects)
objs <- load("EEGpav_group.RData")
print(objs)   # inspect loaded object names

## 2) Load required packages --------------------------------------------
suppressPackageStartupMessages({
  library(lme4)      # glmer / lmer
  library(lmerTest)  # t tests for lmer
  library(car)       # Anova()
})

## 3) Helper functions: summarize glmer / lmer results ------------------

# keep three decimals as character, convenient for tables/manuscripts
fmt3 <- function(x) formatC(x, digits = 3, format = "f")

# for glmer models (logistic)
summ_glmer <- function(mod, name){
  sm <- coef(summary(mod))
  df_fe <- data.frame(
    model = name,
    term  = rownames(sm),
    beta  = sm[, "Estimate"],
    se    = sm[, "Std. Error"],
    stat  = sm[, "z value"],
    p     = sm[, "Pr(>|z|)"],
    row.names = NULL
  )
  # Type-II Wald chi-square test
  aov_tab <- as.data.frame(Anova(mod))
  aov_tab$effect <- rownames(aov_tab)
  rownames(aov_tab) <- NULL
  names(aov_tab)[1:3] <- c("Chisq", "Df", "Pr_Chisq")
  list(fe = df_fe, anova = aov_tab)
}

# for lmer models (continuous DV)
summ_lmer <- function(mod, name){
  sm <- coef(summary(mod))
  df_fe <- data.frame(
    model = name,
    term  = rownames(sm),
    beta  = sm[, "Estimate"],
    se    = sm[, "Std. Error"],
    stat  = sm[, "t value"],
    p     = sm[, "Pr(>|t|)"],
    row.names = NULL
  )
  aov_tab <- as.data.frame(Anova(mod))
  aov_tab$effect <- rownames(aov_tab)
  rownames(aov_tab) <- NULL
  names(aov_tab)[1:3] <- c("Chisq", "Df", "Pr_Chisq")
  list(fe = df_fe, anova = aov_tab)
}

## 4) Summarize each model ----------------------------------------------

res_Mgo        <- summ_glmer(Mgo.basic,      "Mgo.basic_group")
res_Mcorrectgo <- summ_glmer(Mcorrectgo,     "Mcorrectgo_group")
res_Mincorrect <- summ_glmer(Mincorrectgo,   "Mincorrectgo_group")
res_Macc_go    <- summ_glmer(Maccuracyofgo,  "Maccuracyofgo_group")
res_Mrt_basic  <- summ_lmer (Mrt.basic,      "Mrt.basic_group")
res_Mrt_acc    <- summ_lmer (Mrt.accuracy,   "Mrt.accuracy_group")

## 5) Combine all fixed-effects tables ----------------------------------

fe_all <- rbind(
  res_Mgo$fe,
  res_Mcorrectgo$fe,
  res_Mincorrect$fe,
  res_Macc_go$fe,
  res_Mrt_basic$fe,
  res_Mrt_acc$fe
)

## 6) Combine all Anova (chi-square) tables -----------------------------

anova_all <- rbind(
  cbind(model = "Mgo.basic_group",      res_Mgo$anova),
  cbind(model = "Mcorrectgo_group",     res_Mcorrectgo$anova),
  cbind(model = "Mincorrectgo_group",   res_Mincorrect$anova),
  cbind(model = "Maccuracyofgo_group",  res_Macc_go$anova),
  cbind(model = "Mrt.basic_group",      res_Mrt_basic$anova),
  cbind(model = "Mrt.accuracy_group",   res_Mrt_acc$anova)
)

## 7) Keep three decimals -----------------------------------------------

fe_all[, c("beta", "se", "stat", "p")]      <- lapply(fe_all[, c("beta", "se", "stat", "p")], fmt3)
anova_all[, c("Chisq", "Pr_Chisq")]         <- lapply(anova_all[, c("Chisq", "Pr_Chisq")],     fmt3)

## 8) Export to CSV ------------------------------------------------------

write.csv(fe_all,    "EEGpav_mm_group_fixed_effects.csv", row.names = FALSE)
write.csv(anova_all, "EEGpav_mm_group_Anova.csv",         row.names = FALSE)

cat("Export completed. Output files:\n",
    " - EEGpav_mm_group_fixed_effects.csv\n",
    " - EEGpav_mm_group_Anova.csv\n")






## ====================================================================
## Extract results from the joint models with group effects
## saved in EEGpav_group.RData
## and export to two CSV files:
##   - EEGpav_mm_group_fixed_effects.csv
##   - EEGpav_mm_group_Anova.csv
## Path: H:/mt_behavalluse
## Model objects in EEGpav_group.RData:
##   Mgo.joint
##   Mcorrectgo.joint
##   Mincorrectgo.joint
##   Maccuracyofgo.joint
##   Mrt.basic.joint
##   Mrt.accuracy.joint
## ====================================================================

## 1) Clear environment & load objects ----------------------------------
rm(list = ls())

setwd("H:/mt_behavalluse")

# load RData (contains data and saved joint model objects)
objs <- load("EEGpav_group.RData")
print(objs)   # should include *.joint objects

## 2) Load required packages --------------------------------------------
suppressPackageStartupMessages({
  library(lme4)      # glmer / lmer
  library(lmerTest)  # t tests, df for lmer
  library(car)       # Anova()
})

## 3) Helper functions: summarize glmer / lmer results ------------------

# keep three decimals as character, convenient for tables/manuscripts
fmt3 <- function(x) {
  # preserve NA
  out <- ifelse(is.na(x), NA, formatC(x, digits = 3, format = "f"))
  return(out)
}

# for glmer models (logistic)
summ_glmer <- function(mod, name){
  sm <- coef(summary(mod))
  df_fe <- data.frame(
    model = name,
    term  = rownames(sm),
    beta  = sm[, "Estimate"],
    se    = sm[, "Std. Error"],
    stat  = sm[, "z value"],
    p     = sm[, "Pr(>|z|)"],
    row.names = NULL
  )
  # Type-II Wald chi-square test (car::Anova)
  aov_tab <- as.data.frame(Anova(mod))
  aov_tab$effect <- rownames(aov_tab)
  rownames(aov_tab) <- NULL
  # unified column names: Chisq, Df, Pr_Chisq
  names(aov_tab)[1:3] <- c("Chisq", "Df", "Pr_Chisq")
  list(fe = df_fe, anova = aov_tab)
}

# for lmer models (continuous DV)
summ_lmer <- function(mod, name){
  sm <- coef(summary(mod))
  df_fe <- data.frame(
    model = name,
    term  = rownames(sm),
    beta  = sm[, "Estimate"],
    se    = sm[, "Std. Error"],
    stat  = sm[, "t value"],
    p     = sm[, "Pr(>|t|)"],
    row.names = NULL
  )
  aov_tab <- as.data.frame(Anova(mod))
  aov_tab$effect <- rownames(aov_tab)
  rownames(aov_tab) <- NULL
  names(aov_tab)[1:3] <- c("Chisq", "Df", "Pr_Chisq")
  list(fe = df_fe, anova = aov_tab)
}

## 4) Summarize each joint model ----------------------------------------

# note: all *.joint objects are used here
res_Mgo_joint        <- summ_glmer(Mgo.joint,            "Mgo.joint_group")
res_Mcorrectgo_joint <- summ_glmer(Mcorrectgo.joint,     "Mcorrectgo_group")
res_Mincorrect_joint <- summ_glmer(Mincorrectgo.joint,   "Mincorrectgo_group")
res_Macc_go_joint    <- summ_glmer(Maccuracyofgo.joint,  "Maccuracyofgo_group")
res_Mrt_basic_joint  <- summ_lmer (Mrt.basic.joint,      "Mrt.basic_group")
res_Mrt_acc_joint    <- summ_lmer (Mrt.accuracy.joint,   "Mrt.accuracy_group")

## 5) Combine all fixed-effects tables ----------------------------------

fe_all <- rbind(
  res_Mgo_joint$fe,
  res_Mcorrectgo_joint$fe,
  res_Mincorrect_joint$fe,
  res_Macc_go_joint$fe,
  res_Mrt_basic_joint$fe,
  res_Mrt_acc_joint$fe
)

## 6) Combine all Anova (chi-square) tables -----------------------------

anova_all <- rbind(
  cbind(model = "Mgo.joint_group",      res_Mgo_joint$anova),
  cbind(model = "Mcorrectgo_group",     res_Mcorrectgo_joint$anova),
  cbind(model = "Mincorrectgo_group",   res_Mincorrect_joint$anova),
  cbind(model = "Maccuracyofgo_group",  res_Macc_go_joint$anova),
  cbind(model = "Mrt.basic_group",      res_Mrt_basic_joint$anova),
  cbind(model = "Mrt.accuracy_group",   res_Mrt_acc_joint$anova)
)

## 7) Keep three decimals -----------------------------------------------

fe_all[, c("beta", "se", "stat", "p")] <- lapply(
  fe_all[, c("beta", "se", "stat", "p")],
  fmt3
)

anova_all[, c("Chisq", "Pr_Chisq")] <- lapply(
  anova_all[, c("Chisq", "Pr_Chisq")],
  fmt3
)

## 8) Export to CSV ------------------------------------------------------

write.csv(fe_all,    "EEGpav_mm_group_fixed_effects.csv", row.names = FALSE)
write.csv(anova_all, "EEGpav_mm_group_Anova.csv",         row.names = FALSE)

cat("Export completed. Output files:\n",
    " - EEGpav_mm_group_fixed_effects.csv\n",
    " - EEGpav_mm_group_Anova.csv\n")









## =====================================================================
## Group-related effects: P(Go), P(correct Go), P(correct|Go), RT
## Models in EEGpav_group.RData:
##   Mgo.joint, Mcorrectgo.joint, Maccuracyofgo.joint,
##   Mrt.basic.joint, Mrt.accuracy.joint
## Outputs:
##   - CSV tables (emmeans + contrasts with p-values)
##   - Publication-style figures (PNG, 600 dpi)
## =====================================================================

rm(list = ls())

## 0) Paths & global settings -------------------------------------------
rootdir <- "H:/mt_behavalluse"
setwd(rootdir)

outdir <- file.path(rootdir, "group主效应")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

suppressPackageStartupMessages({
  library(lme4)
  library(emmeans)
  library(dplyr)
  library(ggplot2)
})

## Small helper functions ------------------------------------------------
fmt3 <- function(x){
  ifelse(is.na(x), NA, formatC(x, digits = 3, format = "f"))
}

p_to_stars <- function(p){
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01,  "**",
                ifelse(p < 0.05,  "*",  "ns")))
}

## Load models -----------------------------------------------------------
objs <- load("EEGpav_group.RData")
print(objs)

## Shared colors (consistent across all figures) ------------------------
cols_grp   <- c("HC"  = "#E69F00",   # orange
                "OCD" = "#0072B2")   # blue
shapes_grp <- c("HC" = 16, "OCD" = 17)  # HC circle, OCD triangle
pd <- position_dodge(width = 0.25)

## =====================================================================
## (1) P(Go): Mgo.joint · Group × Action
## =====================================================================

emm_go_ga <- emmeans(Mgo.joint, ~ group * action, type = "response")

emm_go_ga_df <- as.data.frame(emm_go_ga) %>%
  mutate(
    Action = factor(ifelse(action == 1, "Go", "NoGo"),
                    levels = c("Go", "NoGo"))
  )

## --- Export emmeans (with formatted columns) ---
emm_go_ga_out <- emm_go_ga_df %>%
  transmute(
    group,
    action,
    Action,
    prob,
    SE,
    asymp.LCL,
    asymp.UCL,
    prob_fmt  = fmt3(prob),
    SE_fmt    = fmt3(SE),
    LCL_fmt   = fmt3(asymp.LCL),
    UCL_fmt   = fmt3(asymp.UCL)
  )

write.csv(
  emm_go_ga_out,
  file = file.path(outdir, "Mgo_joint_groupXaction_pGo_emmeans.csv"),
  row.names = FALSE
)

## --- Simple effect: HC vs OCD within each Action + export ---
emm_group_by_action <- emmeans(Mgo.joint, ~ group | action)
contr_ga <- contrast(emm_group_by_action, method = "pairwise")
contr_ga_df <- as.data.frame(summary(contr_ga)) %>%
  mutate(
    Action  = factor(ifelse(action == 1, "Go", "NoGo"),
                     levels = c("Go", "NoGo")),
    p_stars = p_to_stars(p.value),
    estimate_fmt = fmt3(estimate),
    SE_fmt       = fmt3(SE),
    p_fmt        = fmt3(p.value)
  )

write.csv(
  contr_ga_df,
  file = file.path(outdir, "Mgo_joint_groupXaction_pGo_contrasts_HC_vs_OCD_byAction.csv"),
  row.names = FALSE
)

## --- Compute star labels and line positions ---
max_prob <- emm_go_ga_df %>%
  group_by(Action) %>%
  summarise(ymax = max(prob), .groups = "drop")

sig_df <- contr_ga_df %>%
  select(Action, label = p_stars) %>%
  left_join(max_prob, by = "Action") %>%
  mutate(
    base_y = ymax + 0.07,
    adj_y  = ifelse(Action == "NoGo", ymax + 0.15, base_y),
    y_pos  = pmin(adj_y, 0.98),
    y_line = y_pos - 0.02,
    x_num   = as.numeric(Action),
    x_start = x_num - 0.10,
    x_end   = x_num + 0.10
  )

## --- Plot ---
p_Mgo_ga_pub <- ggplot(
  emm_go_ga_df,
  aes(x = Action,
      y = prob,
      colour = group,
      shape  = group,
      group  = group)
) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                position  = pd,
                width     = 0.08,
                linewidth = 0.6) +
  geom_line(position = pd,
            linewidth = 0.8) +
  geom_point(position = pd,
             size     = 2.6,
             stroke   = 0.8) +
  geom_segment(data = sig_df,
               aes(x = x_start, xend = x_end,
                   y = y_line,  yend = y_line),
               inherit.aes = FALSE,
               colour = "black",
               linewidth = 0.6) +
  geom_text(data = sig_df,
            aes(x = Action, y = y_pos, label = label),
            inherit.aes = FALSE,
            size = 4.2) +
  scale_y_continuous(
    name   = "P(Go)",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0.01)
  ) +
  scale_x_discrete(name = "Action") +
  scale_colour_manual(values = cols_grp, name = NULL) +
  scale_shape_manual(values = shapes_grp, name = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.line        = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.5, colour = "black"),
    axis.text        = element_text(size = 11, colour = "black"),
    axis.title.x     = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y     = element_text(size = 12, margin = margin(r = 8)),
    legend.position      = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.85),
                                     colour = NA),
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    plot.margin      = margin(t = 10, r = 8, b = 8, l = 8)
  )

print(p_Mgo_ga_pub)

ggsave(
  filename = file.path(outdir, "Mgo_joint_groupXaction_pGo_pub.png"),
  plot     = p_Mgo_ga_pub,
  width    = 5.2,
  height   = 4.0,
  dpi      = 600
)

## =====================================================================
## (2) P(correct Go): Mcorrectgo.joint · Group main effect
## =====================================================================

emm_correctGo_g <- emmeans(Mcorrectgo.joint, ~ group, type = "response")
emm_correctGo_g_df <- as.data.frame(emm_correctGo_g)
print(emm_correctGo_g_df)

## emmeans table
emm_correctGo_g_out <- emm_correctGo_g_df %>%
  rename(
    pCorrectGo      = prob,
    pCorrectGo_SE   = SE,
    pCorrectGo_low  = asymp.LCL,
    pCorrectGo_high = asymp.UCL
  ) %>%
  mutate(across(c(pCorrectGo, pCorrectGo_SE,
                  pCorrectGo_low, pCorrectGo_high), fmt3))

write.csv(
  emm_correctGo_g_out,
  file = file.path(outdir, "Mcorrectgo_joint_group_main_pCorrectGo_emmeans.csv"),
  row.names = FALSE
)

## HC vs OCD contrast + p value table
emm_group_correct <- emmeans(Mcorrectgo.joint, ~ group)
contr_group_correct <- contrast(emm_group_correct, method = "pairwise")
contr_group_correct_df <- as.data.frame(summary(contr_group_correct)) %>%
  mutate(
    p_stars = p_to_stars(p.value),
    estimate_fmt = fmt3(estimate),
    SE_fmt       = fmt3(SE),
    p_fmt        = fmt3(p.value)
  )

write.csv(
  contr_group_correct_df,
  file = file.path(outdir, "Mcorrectgo_joint_group_main_pCorrectGo_contrast_HC_minus_OCD.csv"),
  row.names = FALSE
)

star_label_correct <- contr_group_correct_df$p_stars[1]

ymax_c   <- max(emm_correctGo_g_df$prob)
y_line_c <- min(ymax_c + 0.06, 0.98)
y_star_c <- min(y_line_c + 0.03, 0.99)

annot_correct_df <- data.frame(
  x_start = 1,
  x_end   = 2,
  x_star  = 1.5,
  y_line  = y_line_c,
  y_star  = y_star_c,
  label   = star_label_correct
)

p_Mcorrect_pub <- ggplot(
  emm_correctGo_g_df,
  aes(x = group,
      y = prob,
      colour = group,
      shape  = group)
) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width     = 0.10,
                linewidth = 0.6) +
  geom_point(size = 2.8, stroke = 0.8) +
  geom_segment(data = annot_correct_df,
               aes(x = x_start, xend = x_end,
                   y = y_line,  yend = y_line),
               inherit.aes = FALSE,
               colour = "black",
               linewidth = 0.6) +
  geom_text(data = annot_correct_df,
            aes(x = x_star, y = y_star, label = label),
            inherit.aes = FALSE,
            size = 4.2) +
  scale_y_continuous(
    name   = "P(correct Go)",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0.01)
  ) +
  scale_x_discrete(name = "Group") +
  scale_colour_manual(values = cols_grp, name = NULL) +
  scale_shape_manual(values = shapes_grp, name = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.line    = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.5, colour = "black"),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, margin = margin(r = 8)),
    legend.position = "none",
    plot.margin    = margin(t = 10, r = 8, b = 8, l = 8)
  )

print(p_Mcorrect_pub)

ggsave(
  filename = file.path(outdir, "Mcorrectgo_joint_group_main_pCorrectGo_pub.png"),
  plot     = p_Mcorrect_pub,
  width    = 4.0,
  height   = 4.0,
  dpi      = 600
)

## =====================================================================
## (3) P(correct | Go): Maccuracyofgo.joint · Group main effect
## =====================================================================

emm_accGo_g <- emmeans(Maccuracyofgo.joint, ~ group, type = "response")
emm_accGo_g_df <- as.data.frame(emm_accGo_g)
print(emm_accGo_g_df)

## emmeans table
emm_accGo_g_out <- emm_accGo_g_df %>%
  rename(
    pAccuracyGo      = prob,
    pAccuracyGo_SE   = SE,
    pAccuracyGo_low  = asymp.LCL,
    pAccuracyGo_high = asymp.UCL
  ) %>%
  mutate(across(c(pAccuracyGo, pAccuracyGo_SE,
                  pAccuracyGo_low, pAccuracyGo_high), fmt3))

write.csv(
  emm_accGo_g_out,
  file = file.path(outdir, "Maccuracyofgo_joint_group_main_pAccuracyGo_emmeans.csv"),
  row.names = FALSE
)

## HC vs OCD contrast + p value
emm_group_accGo <- emmeans(Maccuracyofgo.joint, ~ group)
contr_group_accGo <- contrast(emm_group_accGo, method = "pairwise")
contr_group_accGo_df <- as.data.frame(summary(contr_group_accGo)) %>%
  mutate(
    p_stars = p_to_stars(p.value),
    estimate_fmt = fmt3(estimate),
    SE_fmt       = fmt3(SE),
    p_fmt        = fmt3(p.value)
  )

write.csv(
  contr_group_accGo_df,
  file = file.path(outdir, "Maccuracyofgo_joint_group_main_pAccuracyGo_contrast_HC_minus_OCD.csv"),
  row.names = FALSE
)

star_label_accGo <- contr_group_accGo_df$p_stars[1]

ymax_a   <- max(emm_accGo_g_df$prob)
y_line_a <- min(ymax_a + 0.06, 0.98)
y_star_a <- min(y_line_a + 0.03, 0.99)

annot_accGo_df <- data.frame(
  x_start = 1,
  x_end   = 2,
  x_star  = 1.5,
  y_line  = y_line_a,
  y_star  = y_star_a,
  label   = star_label_accGo
)

p_MaccGo_pub <- ggplot(
  emm_accGo_g_df,
  aes(x = group,
      y = prob,
      colour = group,
      shape  = group)
) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width     = 0.10,
                linewidth = 0.6) +
  geom_point(size = 2.8, stroke = 0.8) +
  geom_segment(data = annot_accGo_df,
               aes(x = x_start, xend = x_end,
                   y = y_line,  yend = y_line),
               inherit.aes = FALSE,
               colour = "black",
               linewidth = 0.6) +
  geom_text(data = annot_accGo_df,
            aes(x = x_star, y = y_star, label = label),
            inherit.aes = FALSE,
            size = 4.2) +
  scale_y_continuous(
    name   = "P(correct | Go)",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0.01)
  ) +
  scale_x_discrete(name = "Group") +
  scale_colour_manual(values = cols_grp, name = NULL) +
  scale_shape_manual(values = shapes_grp, name = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.line    = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.5, colour = "black"),
    axis.text    = element_text(size = 11, colour = "black"),
    axis.title.x = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y = element_text(size = 12, margin = margin(r = 8)),
    legend.position = "none",
    plot.margin    = margin(t = 10, r = 8, b = 8, l = 8)
  )

print(p_MaccGo_pub)

ggsave(
  filename = file.path(outdir, "Maccuracyofgo_joint_group_main_pAccuracyGo_pub.png"),
  plot     = p_MaccGo_pub,
  width    = 4.0,
  height   = 4.0,
  dpi      = 600
)

## =====================================================================
## (4) RT: Mrt.basic.joint · Group × Action
## =====================================================================

emm_rt_ga <- emmeans(Mrt.basic.joint, ~ group * action)

emm_rt_ga_df <- summary(emm_rt_ga, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(
    rt_sec      = exp(emmean)    - 0.9,
    rt_sec_low  = exp(asymp.LCL) - 0.9,
    rt_sec_high = exp(asymp.UCL) - 0.9,
    Action = factor(ifelse(action == 1, "Go", "NoGo"),
                    levels = c("Go", "NoGo"))
  )

## emmeans table
emm_rt_ga_out <- emm_rt_ga_df %>%
  transmute(
    group,
    action,
    Action,
    emmean_lnRT    = fmt3(emmean),
    asymp.LCL_lnRT = fmt3(asymp.LCL),
    asymp.UCL_lnRT = fmt3(asymp.UCL),
    rt_sec         = fmt3(rt_sec),
    rt_sec_low     = fmt3(rt_sec_low),
    rt_sec_high    = fmt3(rt_sec_high)
  )

write.csv(
  emm_rt_ga_out,
  file = file.path(outdir, "Mrt_basic_joint_groupXaction_RT_emmeans.csv"),
  row.names = FALSE
)

## HC vs OCD by Action contrast + p values
emm_rt_group_by_action <- emmeans(Mrt.basic.joint, ~ group | action)
contr_rt_ga <- contrast(emm_rt_group_by_action, method = "pairwise")
contr_rt_df <- as.data.frame(summary(contr_rt_ga)) %>%
  mutate(
    Action = factor(ifelse(action == 1, "Go", "NoGo"),
                    levels = c("Go", "NoGo")),
    p_stars      = p_to_stars(p.value),
    estimate_fmt = fmt3(estimate),
    SE_fmt       = fmt3(SE),
    p_fmt        = fmt3(p.value)
  )

write.csv(
  contr_rt_df,
  file = file.path(outdir, "Mrt_basic_joint_groupXaction_RT_contrasts_HC_vs_OCD_byAction.csv"),
  row.names = FALSE
)

## star labels and horizontal lines
max_rt <- emm_rt_ga_df %>%
  group_by(Action) %>%
  summarise(ymax = max(rt_sec), .groups = "drop")

annot_rt_df <- contr_rt_df %>%
  select(Action, label = p_stars) %>%
  left_join(max_rt, by = "Action") %>%
  mutate(
    x_num   = as.numeric(Action),
    x_start = x_num - 0.18,
    x_end   = x_num + 0.18,
    x_star  = x_num,
    y_line  = ymax + 0.05,
    y_star  = ymax + 0.09
  )

p_Mrt_ga_pub <- ggplot(
  emm_rt_ga_df,
  aes(x = Action,
      y = rt_sec,
      colour = group,
      shape  = group,
      group  = group)
) +
  geom_errorbar(aes(ymin = rt_sec_low, ymax = rt_sec_high),
                position  = pd,
                width     = 0.08,
                linewidth = 0.6) +
  geom_line(position = pd,
            linewidth = 0.8) +
  geom_point(position = pd,
             size     = 2.6,
             stroke   = 0.8) +
  geom_segment(data = annot_rt_df,
               aes(x = x_start, xend = x_end,
                   y = y_line,  yend = y_line),
               inherit.aes = FALSE,
               colour = "black",
               linewidth = 0.6) +
  geom_text(data = annot_rt_df,
            aes(x = x_star, y = y_star, label = label),
            inherit.aes = FALSE,
            size = 4.2) +
  scale_y_continuous(
    name   = "RT (s)",
    breaks = scales::pretty_breaks(5),
    expand = c(0, 0.02)
  ) +
  scale_x_discrete(name = "Action") +
  scale_colour_manual(values = cols_grp, name = NULL) +
  scale_shape_manual(values = shapes_grp, name = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.line        = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.5, colour = "black"),
    axis.text        = element_text(size = 11, colour = "black"),
    axis.title.x     = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y     = element_text(size = 12, margin = margin(r = 8)),
    legend.position      = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.85),
                                     colour = NA),
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    plot.margin      = margin(t = 10, r = 8, b = 8, l = 8)
  )

print(p_Mrt_ga_pub)

ggsave(
  filename = file.path(outdir, "Mrt_basic_joint_groupXaction_RT_pub.png"),
  plot     = p_Mrt_ga_pub,
  width    = 5.2,
  height   = 4.0,
  dpi      = 600
)

## =====================================================================
## (5) RT: Mrt.accuracy.joint · Group × Accuracy (Go trials only)
## =====================================================================

emm_rt_gaAcc <- emmeans(Mrt.accuracy.joint, ~ group * accuracy)

emm_rt_gaAcc_df <- summary(emm_rt_gaAcc, infer = TRUE) %>%
  as.data.frame() %>%
  mutate(
    rt_sec      = exp(emmean)    - 0.9,
    rt_sec_low  = exp(asymp.LCL) - 0.9,
    rt_sec_high = exp(asymp.UCL) - 0.9,
    Accuracy    = factor(ifelse(accuracy == 1, "Correct", "Error"),
                         levels = c("Correct", "Error"))
  )

emm_rt_gaAcc_out <- emm_rt_gaAcc_df %>%
  transmute(
    group,
    accuracy,
    Accuracy,
    emmean_lnRT    = fmt3(emmean),
    asymp.LCL_lnRT = fmt3(asymp.LCL),
    asymp.UCL_lnRT = fmt3(asymp.UCL),
    rt_sec         = fmt3(rt_sec),
    rt_sec_low     = fmt3(rt_sec_low),
    rt_sec_high    = fmt3(rt_sec_high)
  )

write.csv(
  emm_rt_gaAcc_out,
  file = file.path(outdir, "Mrt_accuracy_joint_groupXaccuracy_RT_emmeans.csv"),
  row.names = FALSE
)

## HC vs OCD by Accuracy + p values
emm_rt_group_by_acc <- emmeans(Mrt.accuracy.joint, ~ group | accuracy)
contr_rt_acc <- contrast(emm_rt_group_by_acc, method = "pairwise")
contr_rt_acc_df <- as.data.frame(summary(contr_rt_acc)) %>%
  mutate(
    Accuracy    = factor(ifelse(accuracy == 1, "Correct", "Error"),
                         levels = c("Correct", "Error")),
    p_stars      = p_to_stars(p.value),
    estimate_fmt = fmt3(estimate),
    SE_fmt       = fmt3(SE),
    p_fmt        = fmt3(p.value)
  )

write.csv(
  contr_rt_acc_df,
  file = file.path(outdir, "Mrt_accuracy_joint_groupXaccuracy_RT_contrasts_HC_vs_OCD_byAccuracy.csv"),
  row.names = FALSE
)

max_rt_acc <- emm_rt_gaAcc_df %>%
  group_by(Accuracy) %>%
  summarise(ymax = max(rt_sec), .groups = "drop")

annot_rt_acc_df <- contr_rt_acc_df %>%
  select(Accuracy, label = p_stars) %>%
  left_join(max_rt_acc, by = "Accuracy") %>%
  mutate(
    x_num   = as.numeric(Accuracy),
    x_start = x_num - 0.18,
    x_end   = x_num + 0.18,
    x_star  = x_num,
    y_line  = ymax + 0.05,
    y_star  = ymax + 0.09
  )

p_Mrt_gaAcc_pub <- ggplot(
  emm_rt_gaAcc_df,
  aes(x = Accuracy,
      y = rt_sec,
      colour = group,
      shape  = group,
      group  = group)
) +
  geom_errorbar(aes(ymin = rt_sec_low, ymax = rt_sec_high),
                position  = pd,
                width     = 0.08,
                linewidth = 0.6) +
  geom_line(position = pd,
            linewidth = 0.8) +
  geom_point(position = pd,
             size     = 2.6,
             stroke   = 0.8) +
  geom_segment(data = annot_rt_acc_df,
               aes(x = x_start, xend = x_end,
                   y = y_line,  yend = y_line),
               inherit.aes = FALSE,
               colour = "black",
               linewidth = 0.6) +
  geom_text(data = annot_rt_acc_df,
            aes(x = x_star, y = y_star, label = label),
            inherit.aes = FALSE,
            size = 4.2) +
  scale_y_continuous(
    name   = "RT (s)",
    breaks = scales::pretty_breaks(5),
    expand = c(0, 0.02)
  ) +
  scale_x_discrete(name = "Accuracy") +
  scale_colour_manual(values = cols_grp, name = NULL) +
  scale_shape_manual(values = shapes_grp, name = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.line        = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks       = element_line(linewidth = 0.5, colour = "black"),
    axis.text        = element_text(size = 11, colour = "black"),
    axis.title.x     = element_text(size = 12, margin = margin(t = 8)),
    axis.title.y     = element_text(size = 12, margin = margin(r = 8)),
    legend.position      = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.85),
                                     colour = NA),
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    plot.margin      = margin(t = 10, r = 8, b = 8, l = 8)
  )

print(p_Mrt_gaAcc_pub)

ggsave(
  filename = file.path(outdir, "Mrt_accuracy_joint_groupXaccuracy_RT_pub.png"),
  plot     = p_Mrt_gaAcc_pub,
  width    = 5.2,
  height   = 4.0,
  dpi      = 600
)

## =====================================================================
cat("Group-related main effects and interaction visualizations completed.\nOutput directory:\n",
    outdir, "\n")