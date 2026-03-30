GET DATA  /TYPE=TXT
  /FILE="H:\MT_all(afterqujizhi)\feedback\feedback2\OCD\fb4spss.csv"
  /ENCODING='Locale'
  /DELCASE=LINE
  /DELIMITERS=","
  /ARRANGEMENT=DELIMITED
  /FIRSTCASE=2
  /IMPORTCASE=ALL
  /VARIABLES=
  subject F2.0
  goNoWin F9.7
  nogoNoWin A9
  goPun F10.7
  nogoPun F9.7
  goWin F10.7
  nogoWin F9.7
  goNoPun F10.7
  nogoNoPun F8.6.
CACHE.
EXECUTE.
DATASET NAME DataSet1 WINDOW=FRONT.
ALTER TYPE  nogoNoWin(F27.9).
FORMATS  nogoNoWin(F27.9).
VARIABLE LEVEL  nogoNoWin(SCALE).
EXECUTE.

GLM goNoWin nogoNoWin goPun nogoPun goWin nogoWin goNoPun nogoNoPun
  /WSFACTOR=outcome 2 Polynomial valence 2 Polynomial response 2 Polynomial 
  /METHOD=SSTYPE(3)
  /PLOT=PROFILE(valence*response*outcome)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=outcome valence response outcome*valence outcome*response valence*response 
    outcome*valence*response.

* simple effects - Win cues.
GLM goNoWin nogoNoWin goWin nogoWin
  /WSFACTOR=outcome 2 Polynomial response 2 Polynomial 
  /METHOD=SSTYPE(3)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=outcome response outcome*response.
* simple effects - Avoid cues.
GLM goPun nogoPun goNoPun nogoNoPun
  /WSFACTOR=outcome 2 Polynomial response 2 Polynomial 
  /METHOD=SSTYPE(3)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=outcome response outcome*response.

* planned contrasts:
i) biased vs. unbiased learning (outcomexresponse)
ii) enhanced vs. reduced learning (response).
GLM goPun nogoPun goWin nogoWin
  /WSFACTOR=outcome 2 Polynomial response 2 Polynomial 
  /METHOD=SSTYPE(3)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=outcome response outcome*response.

* planned contrasts: neutral outcomes for win vs. avoid cues.
GLM goNoWin nogoNoWin goNoPun nogoNoPun
  /WSFACTOR=outcome 2 Polynomial response 2 Polynomial 
  /METHOD=SSTYPE(3)
  /CRITERIA=ALPHA(.05)
  /WSDESIGN=outcome response outcome*response.

