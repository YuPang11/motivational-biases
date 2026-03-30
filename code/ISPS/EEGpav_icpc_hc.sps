* Encoding: UTF-8.

GET DATA  /TYPE=TXT
  /FILE="F:\PIT\correct_all\HC\results\ICPC_contra_ipsiocd0.4.65plottest.csv"
  /ENCODING='Locale'
  /DELCASE=LINE
  /DELIMITERS=","
  /ARRANGEMENT=DELIMITED
  /FIRSTCASE=2
  /IMPORTCASE=ALL
  /VARIABLES=
  subject F2.0
  icpc_pfc_win_contra A11
  icpc_pfc_win_ipsi A11
  icpc_pfc_win_nogo A11
  icpc_pfc_avoid_contra F11.8
  icpc_pfc_avoid_ipsi F11.8
  icpc_pfc_avoid_nogo F11.8
  icpc_par_win_contra A11
  icpc_par_win_ipsi A12
  icpc_par_win_nogo A12
  icpc_par_avoid_contra F11.8
  icpc_par_avoid_ipsi F11.9
  icpc_par_avoid_nogo F10.8.
CACHE.
EXECUTE.
DATASET NAME DataSet1 WINDOW=FRONT.

ALTER TYPE  icpc_pfc_win_contra(F33.9).
ALTER TYPE  icpc_pfc_win_ipsi(F33.9).
ALTER TYPE  icpc_pfc_win_nogo(F33.9).
ALTER TYPE  icpc_par_win_contra(F33.9).
ALTER TYPE  icpc_par_win_ipsi(F33.9).
ALTER TYPE  icpc_par_win_nogo(F33.9).
FORMATS  icpc_pfc_win_contra(F33.9).
FORMATS  icpc_pfc_win_ipsi(F33.9).
FORMATS  icpc_pfc_win_nogo(F33.9).
FORMATS  icpc_par_win_contra(F33.9).
FORMATS  icpc_par_win_ipsi(F33.9).
FORMATS  icpc_par_win_nogo(F33.9).
VARIABLE LEVEL  icpc_pfc_win_contra(SCALE).
VARIABLE LEVEL  icpc_pfc_win_ipsi(SCALE).
VARIABLE LEVEL  icpc_pfc_win_nogo(SCALE).
VARIABLE LEVEL  icpc_par_win_contra(SCALE).
VARIABLE LEVEL  icpc_par_win_ipsi(SCALE).
VARIABLE LEVEL  icpc_par_win_nogo(SCALE).
EXECUTE.

COMPUTE icpc_par_win_lateralization=icpc_par_win_contra - icpc_par_win_ipsi.
COMPUTE icpc_par_avoid_lateralization=icpc_par_avoid_contra - icpc_par_avoid_ipsi.
COMPUTE icpc_par_win_go=(icpc_par_win_contra + icpc_par_win_ipsi) /2.
COMPUTE icpc_par_avoid_go=(icpc_par_avoid_contra + icpc_par_avoid_ipsi) /2.
COMPUTE icpc_pfc_win_go=(icpc_pfc_win_contra + icpc_pfc_win_ipsi) /2.
COMPUTE icpc_pfc_avoid_go=(icpc_pfc_avoid_contra + icpc_pfc_avoid_ipsi) /2.
COMPUTE icpc_par_congr=(icpc_par_avoid_nogo + icpc_par_win_ipsi) /2.
COMPUTE icpc_par_incongr=(icpc_par_avoid_ipsi + icpc_par_win_nogo) /2.
EXECUTE.

* midfrontal-parietal ICPC: valence x response(nogo/contra-/ipsilateral).
GLM icpc_par_win_nogo icpc_par_win_contra icpc_par_win_ipsi  
    icpc_par_avoid_nogo icpc_par_avoid_contra icpc_par_avoid_ipsi 
  /WSFACTOR=valence 2 Simple response 3 Simple(1) 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(response*valence)
  /PRINT=DESCRIPTIVE ETASQ
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=valence response valence*response.

* rephrase in terms of congruency x response, where response indicates contralateral vs. ipsilateral&nogo.
GLM icpc_par_congr icpc_par_win_contra icpc_par_incongr icpc_par_avoid_contra
  /WSFACTOR=congruency 2 Simple(1) response 2 Simple 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(congruency*response)
  /PRINT=DESCRIPTIVE ETASQ
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=congruency response congruency*response.

* simple effect for non-executing sites.
GLM icpc_par_win_ipsi icpc_par_win_nogo icpc_par_avoid_ipsi icpc_par_avoid_nogo
  /WSFACTOR=valence 2 Polynomial response 2 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(response*valence)
  /PRINT=DESCRIPTIVE ETASQ
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=valence response valence*response.

* simple effects of valence x response interaction for barplot.
T-TEST PAIRS=icpc_par_win_contra icpc_par_win_ipsi icpc_par_win_nogo 
      WITH icpc_par_avoid_contra icpc_par_avoid_ipsi icpc_par_avoid_nogo (PAIRED)
  /CRITERIA=CI(.9500)
  /MISSING=ANALYSIS.

* midfrontal-dlPFC ICPC: valence x response(nogo/go).
GLM icpc_pfc_win_nogo icpc_pfc_win_go icpc_pfc_avoid_nogo icpc_pfc_avoid_go
  /WSFACTOR=valence 2 Polynomial action 2 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(action*valence)
  /PRINT=DESCRIPTIVE ETASQ 
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=valence action valence*action.

* simple effects of valence x response interaction.
T-TEST PAIRS=icpc_pfc_win_nogo icpc_pfc_win_go 
      WITH icpc_pfc_avoid_nogo icpc_pfc_avoid_go (PAIRED)
  /CRITERIA=CI(.9500)
  /MISSING=ANALYSIS.


* incongr vs congr.
T-TEST PAIRS=icpc_par_incongr WITH icpc_par_congr (PAIRED)
 /CRITERIA=CI(.9500)
 /MISSING=ANALYSIS.

* contra：incongr vs congr.
T-TEST PAIRS=icpc_par_avoid_contra WITH icpc_par_win_contra (PAIRED)
 /CRITERIA=CI(.9500)
 /MISSING=ANALYSIS.


* -------------------------------------------------------------------- .
* Export Output Results to Excel (Fixed for your version) .
* -------------------------------------------------------------------- .

OUTPUT EXPORT
  /CONTENTS EXPORT=ALL
  /XLSX DOCUMENTFILE="F:\PIT\correct_all\hc0.4.65.xlsx".
