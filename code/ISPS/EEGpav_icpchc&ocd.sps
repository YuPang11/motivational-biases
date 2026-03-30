* Encoding: UTF-8.
GET DATA  /TYPE=TXT
  /FILE="F:\PIT\correct_all\ISPC_hc_ocd_compare.csv"
  /ENCODING='Locale'
  /DELCASE=LINE
  /DELIMITERS=","
  /ARRANGEMENT=DELIMITED
  /FIRSTCASE=2
  /IMPORTCASE=ALL
  /VARIABLES=
    subject               F3.0
    group                 F1.0
    icpc_pfc_win_contra   F11.8
    icpc_pfc_win_ipsi     F11.8
    icpc_pfc_win_nogo     F11.8
    icpc_pfc_avoid_contra F11.8
    icpc_pfc_avoid_ipsi   F11.8
    icpc_pfc_avoid_nogo   F11.8
    icpc_par_win_contra   F11.8
    icpc_par_win_ipsi     F11.8
    icpc_par_win_nogo     F11.8
    icpc_par_avoid_contra F11.8
    icpc_par_avoid_ipsi   F11.8
    icpc_par_avoid_nogo   F11.8.
CACHE.
EXECUTE.
DATASET NAME DataSet1 WINDOW=FRONT.

VALUE LABELS group
  0 'HC'
  1 'OCD'.
VARIABLE LEVEL group (NOMINAL).
EXECUTE.

COMPUTE icpc_par_win_lateralization  = icpc_par_win_contra - icpc_par_win_ipsi.
COMPUTE icpc_par_avoid_lateralization= icpc_par_avoid_contra - icpc_par_avoid_ipsi.
COMPUTE icpc_par_win_go              = (icpc_par_win_contra + icpc_par_win_ipsi) / 2.
COMPUTE icpc_par_avoid_go            = (icpc_par_avoid_contra + icpc_par_avoid_ipsi) / 2.
COMPUTE icpc_pfc_win_go              = (icpc_pfc_win_contra + icpc_pfc_win_ipsi) / 2.
COMPUTE icpc_pfc_avoid_go            = (icpc_pfc_avoid_contra + icpc_pfc_avoid_ipsi) / 2.
COMPUTE icpc_par_congr               = (icpc_par_avoid_nogo + icpc_par_win_ipsi) / 2.
COMPUTE icpc_par_incongr             = (icpc_par_avoid_ipsi + icpc_par_win_nogo) / 2.
EXECUTE.

*-------------------------------------------------------------------*.
* GLM 1: midfrontal–parietal ICPC · valence (Win/Avoid) × response (nogo/contra/ipsi) × group.
*-------------------------------------------------------------------*.

GLM icpc_par_win_nogo icpc_par_win_contra icpc_par_win_ipsi
    icpc_par_avoid_nogo icpc_par_avoid_contra icpc_par_avoid_ipsi
  BY group
  /WSFACTOR = valence 2 Simple response 3 Simple(1)
  /METHOD   = SSTYPE(3)
  /PLOT     = PROFILE(response*valence*group)
  /CRITERIA = ALPHA(.05)
  /WSDESIGN = valence response valence*response
  /DESIGN   = group.

*-------------------------------------------------------------------*.
* GLM 2: midfrontal–parietal ICPC · congruency (congr/incongr) × response (executing vs non-executing) × group.
*   - congruency:  icpc_par_congr vs icpc_par_incongr
*   - response:    executing(=contra) vs non-executing(=ipsi&nogo)
*-------------------------------------------------------------------*.

GLM icpc_par_congr icpc_par_win_contra
    icpc_par_incongr icpc_par_avoid_contra
  BY group
  /WSFACTOR = congruency 2 Simple(1) response 2 Simple
  /METHOD   = SSTYPE(3)
  /PLOT     = PROFILE(congruency*response*group)
  /CRITERIA = ALPHA(.05)
  /WSDESIGN = congruency response congruency*response
  /DESIGN   = group.

*-------------------------------------------------------------------*.
* GLM 3: midfrontal–parietal ICPC (ipsi / nogo)，valence × response × group.
*-------------------------------------------------------------------*.

GLM icpc_par_win_ipsi icpc_par_win_nogo
    icpc_par_avoid_ipsi icpc_par_avoid_nogo
  BY group
  /WSFACTOR = valence 2 Polynomial response 2 Polynomial
  /METHOD   = SSTYPE(3)
  /PLOT     = PROFILE(response*valence*group)
  /CRITERIA = ALPHA(.05)
  /WSDESIGN = valence response valence*response
  /DESIGN   = group.

*-------------------------------------------------------------------*.
* GLM 4: midfrontal–dlPFC ICPC · valence (Win/Avoid) × action (Go/NoGo) × group.
*-------------------------------------------------------------------*.

GLM icpc_pfc_win_nogo icpc_pfc_win_go
    icpc_pfc_avoid_nogo icpc_pfc_avoid_go
  BY group
  /WSFACTOR = valence 2 Polynomial action 2 Polynomial
  /METHOD   = SSTYPE(3)
  /PLOT     = PROFILE(action*valence*group)
  /PRINT    = DESCRIPTIVE ETASQ
  /CRITERIA = ALPHA(.05)
  /WSDESIGN = valence action valence*action
  /DESIGN   = group.
