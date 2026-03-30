* ============================================================
* 0) Read the original long-format Excel file
* ============================================================.

GET DATA
  /TYPE=XLSX
  /FILE="H:\mt_behavalluse\EEGmixmodel\EEGdatamixmodel_withgroup.xlsx"
  /SHEET=index 1
  /CELLRANGE=FULL
  /READNAMES=ON.
CACHE.
EXECUTE.

DATASET NAME rawdata WINDOW=FRONT.

* ============================================================
* 1) Check variable types
* ============================================================.
DISPLAY DICTIONARY.

* ============================================================
* 2) Force key EEG variables to numeric
*    This is also safe if they are already numeric
* ============================================================.
ALTER TYPE MFpower ICPCpfc l_ICPCmotor r_ICPCmotor (F12.6).
EXECUTE.

* ============================================================
* A) HC analysis
*    Assume group = 0 is HC
* ============================================================.

DATASET COPY hcdata.
DATASET ACTIVATE hcdata.

SELECT IF (group = 0).
SELECT IF (trial2use = 1).
EXECUTE.

* -------- Create trial-level temporary variables for 8 conditions --------.
COMPUTE GW_corr_tmp  = $SYSMIS.
COMPUTE GA_corr_tmp  = $SYSMIS.
COMPUTE NGW_corr_tmp = $SYSMIS.
COMPUTE NGA_corr_tmp = $SYSMIS.
COMPUTE GW_err_tmp   = $SYSMIS.
COMPUTE GA_err_tmp   = $SYSMIS.
COMPUTE NGW_err_tmp  = $SYSMIS.
COMPUTE NGA_err_tmp  = $SYSMIS.

* Correct trials.
IF (action = 1 AND valence = 1 AND DVacc = 1) GW_corr_tmp  = MFpower.
IF (action = 1 AND valence = 0 AND DVacc = 1) GA_corr_tmp  = MFpower.
IF (action = 0 AND valence = 1 AND DVacc = 1) NGW_corr_tmp = MFpower.
IF (action = 0 AND valence = 0 AND DVacc = 1) NGA_corr_tmp = MFpower.

* Error trials.
IF (action = 1 AND valence = 1 AND DVacc = 0) GW_err_tmp   = MFpower.
IF (action = 1 AND valence = 0 AND DVacc = 0) GA_err_tmp   = MFpower.
IF (action = 0 AND valence = 1 AND DVacc = 0) NGW_err_tmp  = MFpower.
IF (action = 0 AND valence = 0 AND DVacc = 0) NGA_err_tmp  = MFpower.
EXECUTE.

* -------- Aggregate to wide format by subject --------.
DATASET DECLARE agg_hc.

AGGREGATE
  /OUTFILE=agg_hc
  /BREAK=sID
  /GW_corr  = MEAN(GW_corr_tmp)
  /GA_corr  = MEAN(GA_corr_tmp)
  /NGW_corr = MEAN(NGW_corr_tmp)
  /NGA_corr = MEAN(NGA_corr_tmp)
  /GW_err   = MEAN(GW_err_tmp)
  /GA_err   = MEAN(GA_err_tmp)
  /NGW_err  = MEAN(NGW_err_tmp)
  /NGA_err  = MEAN(NGA_err_tmp).

DATASET ACTIVATE agg_hc.

FORMATS GW_corr GA_corr NGW_corr NGA_corr GW_err GA_err NGW_err NGA_err (F9.6).
VARIABLE LEVEL GW_corr GA_corr NGW_corr NGA_corr GW_err GA_err NGW_err NGA_err (SCALE).
EXECUTE.

* -------- Save aggregated HC wide data to Excel --------.
SAVE TRANSLATE
  /OUTFILE="H:\mt_behavalluse\EEGmixmodel\HC_MFpower_wide.xlsx"
  /TYPE=XLSX
  /VERSION=12
  /FIELDNAMES
  /REPLACE.

* -------- Create a new HC output window --------.
OUTPUT NEW NAME=HC_Output.
OUTPUT ACTIVATE HC_Output.

* -------- HC: 2×2×2 repeated-measures GLM --------.
GLM GW_corr GA_corr NGW_corr NGA_corr GW_err GA_err NGW_err NGA_err
  /WSFACTOR=accuracy 2 Polynomial action 2 Polynomial valence 2 Polynomial
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(valence*action*accuracy)
  /PRINT=DESCRIPTIVE ETASQ
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=accuracy action valence
            accuracy*action accuracy*valence action*valence
            accuracy*action*valence.

* -------- HC: simple effect - correct trials --------.
GLM GW_corr GA_corr NGW_corr NGA_corr
  /WSFACTOR=action 2 Polynomial valence 2 Polynomial
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(valence*action)
  /PRINT=DESCRIPTIVE ETASQ
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=action valence action*valence.

* -------- HC: simple effect - error trials --------.
GLM GW_err GA_err NGW_err NGA_err
  /WSFACTOR=action 2 Polynomial valence 2 Polynomial
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(valence*action)
  /PRINT=DESCRIPTIVE ETASQ
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=action valence action*valence.

* -------- Export HC results to Excel --------.
OUTPUT EXPORT
  /CONTENTS EXPORT=ALL
  /XLSX DOCUMENTFILE="H:\mt_behavalluse\EEGmixmodel\HC_MFpower_results.xlsx".

OUTPUT CLOSE HC_Output.



* ============================================================
* B) OCD analysis
*    Assume group = 1 is OCD
* ============================================================.

DATASET ACTIVATE rawdata.
DATASET COPY ocddata.
DATASET ACTIVATE ocddata.

SELECT IF (group = 1).
SELECT IF (trial2use = 1).
EXECUTE.

* -------- Create trial-level temporary variables for 8 conditions --------.
COMPUTE GW_corr_tmp  = $SYSMIS.
COMPUTE GA_corr_tmp  = $SYSMIS.
COMPUTE NGW_corr_tmp = $SYSMIS.
COMPUTE NGA_corr_tmp = $SYSMIS.
COMPUTE GW_err_tmp   = $SYSMIS.
COMPUTE GA_err_tmp   = $SYSMIS.
COMPUTE NGW_err_tmp  = $SYSMIS.
COMPUTE NGA_err_tmp  = $SYSMIS.

* Correct trials.
IF (action = 1 AND valence = 1 AND DVacc = 1) GW_corr_tmp  = MFpower.
IF (action = 1 AND valence = 0 AND DVacc = 1) GA_corr_tmp  = MFpower.
IF (action = 0 AND valence = 1 AND DVacc = 1) NGW_corr_tmp = MFpower.
IF (action = 0 AND valence = 0 AND DVacc = 1) NGA_corr_tmp = MFpower.

* Error trials.
IF (action = 1 AND valence = 1 AND DVacc = 0) GW_err_tmp   = MFpower.
IF (action = 1 AND valence = 0 AND DVacc = 0) GA_err_tmp   = MFpower.
IF (action = 0 AND valence = 1 AND DVacc = 0) NGW_err_tmp  = MFpower.
IF (action = 0 AND valence = 0 AND DVacc = 0) NGA_err_tmp  = MFpower.
EXECUTE.

* -------- Aggregate to wide format by subject --------.
DATASET DECLARE agg_ocd.

AGGREGATE
  /OUTFILE=agg_ocd
  /BREAK=sID
  /GW_corr  = MEAN(GW_corr_tmp)
  /GA_corr  = MEAN(GA_corr_tmp)
  /NGW_corr = MEAN(NGW_corr_tmp)
  /NGA_corr = MEAN(NGA_corr_tmp)
  /GW_err   = MEAN(GW_err_tmp)
  /GA_err   = MEAN(GA_err_tmp)
  /NGW_err  = MEAN(NGW_err_tmp)
  /NGA_err  = MEAN(NGA_err_tmp).

DATASET ACTIVATE agg_ocd.

FORMATS GW_corr GA_corr NGW_corr NGA_corr GW_err GA_err NGW_err NGA_err (F9.6).
VARIABLE LEVEL GW_corr GA_corr NGW_corr NGA_corr GW_err GA_err NGW_err NGA_err (SCALE).
EXECUTE.

* -------- Save aggregated OCD wide data to Excel --------.
SAVE TRANSLATE
  /OUTFILE="H:\mt_behavalluse\EEGmixmodel\OCD_MFpower_wide.xlsx"
  /TYPE=XLSX
  /VERSION=12
  /FIELDNAMES
  /REPLACE.

* -------- Create a new OCD output window --------.
OUTPUT NEW NAME=OCD_Output.
OUTPUT ACTIVATE OCD_Output.

* -------- OCD: 2×2×2 repeated-measures GLM --------.
GLM GW_corr GA_corr NGW_corr NGA_corr GW_err GA_err NGW_err NGA_err
  /WSFACTOR=accuracy 2 Polynomial action 2 Polynomial valence 2 Polynomial
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(valence*action*accuracy)
  /PRINT=DESCRIPTIVE ETASQ
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=accuracy action valence
            accuracy*action accuracy*valence action*valence
            accuracy*action*valence.

* -------- OCD: simple effect - correct trials --------.
GLM GW_corr GA_corr NGW_corr NGA_corr
  /WSFACTOR=action 2 Polynomial valence 2 Polynomial
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(valence*action)
  /PRINT=DESCRIPTIVE ETASQ
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=action valence action*valence.

* -------- OCD: simple effect - error trials --------.
GLM GW_err GA_err NGW_err NGA_err
  /WSFACTOR=action 2 Polynomial valence 2 Polynomial
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(valence*action)
  /PRINT=DESCRIPTIVE ETASQ
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=action valence action*valence.

* -------- Export OCD results to Excel --------.
OUTPUT EXPORT
  /CONTENTS EXPORT=ALL
  /XLSX DOCUMENTFILE="H:\mt_behavalluse\EEGmixmodel\OCD_MFpower_results.xlsx".

OUTPUT CLOSE OCD_Output.