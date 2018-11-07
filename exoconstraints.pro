PRO exoconstraints
;;A Geo-&Cosmo-model [CODE-1] for estimating rocky exoplanetary
;;composition and structures: mantle and core composition, core mass
;;fraction, core radius.
;;Created and maintained by Haiyang S. Wang
;;Citation: Wang, H. S., Liu, F., Ireland, T., Brasser, R., Yong, D., and Lineweaver, C. H. 2018. Enhanced constraints on the interior composition and structure of terrestrial exoplanets. MNRAS, in press. doi.org/10.1093/mnras/sty2749

TIC
cgCleanup ;;This procedure cleans-up and/or destroys any open graphics or widget windows on the display.
;;;>>>>>>>KEYS SET UP
COTcCorrectlable='YES'
solarHgCorrectlable='YES'
sigmapplylable='0' ;;;1sigma pattern of the DP: 0-using the extrem point within 68% contour; ELSE-using the 1sigma of the coefficients with the constraint of the T_D range
starkeyRlable='YES'
modelEarthkeyR='YES'
ratiolmtshowlable='NO'
plotlabel='NO'

Nelems=83
nan=!values.f_nan

;;;>>>INPUT NECESSARY SUPPORTING DATA
readcol, 'data/atomwttc_new.txt', F='I,A,F,F', atomid, elemid, atomwt, elemtc
readcol, 'data/protosunppmwhy.txt', F='x,d,d,d', solarppm, solarppmerrup, solarppmerrdn
readcol, 'data/protosunwhy.txt', F='x,f,f,f', protosundex, protosundexerrup, protosundexerrdn
readcol, 'data/PEppmwhy.txt', F='x,d,d,d', bulkearthppm, bulkearthppmerrup, bulkearthppmerrdn
readcol, 'data/Asplund09dex.txt', F='x,x,f,f', A09solardex, A09solardexerr, skipline=3
;;--Example of stars----
readcol, 'data/Kepler10_liu16.txt', F='x,f,f,f,f', k10HETdex, k10HETdexerr, k10CFHTdex, k10CFHTdexerr
readcol, 'data/Extrasolarabu.txt', F='x,f,f,f,f,f,f', K20dex, K20dexerr, K21dex, K21dexerr, K100dex, K100dexerr, skipline=2
A09solardex(where(A09solardex eq -11))=nan
A09solardexerr(where(A09solardexerr eq -11))=nan
solardex=A09solardex
solardexerr=A09solardexerr

K10HETdex(where(K10HETdex eq -11))=nan
K10HETdexerr(where(K10HETdexerr eq -11))=nan
K20dex(where(K20dex eq -11))=nan
K20dexerr(where(K20dexerr eq -11))=nan
K21dex(where(K21dex eq -11))=nan
K21dexerr(where(K21dexerr eq -11))=nan
K100dex(where(K100dex eq -11))=nan
K100dexerr(where(K100dexerr eq -11))=nan


;STOP
for i=0, Nelems-1 do begin
 if elemid[i] eq 'C' then nC=i
 if elemid[i] eq 'O' then nO=i
 if elemid[i] eq 'S' then nS=i
 if elemid[i] eq 'Hg' then nHg=i
 if elemid[i] eq 'Na' then nNa=i
 if elemid[i] eq 'Si' then nSi=i
 if elemid[i] eq 'Fe' then nFe=i
 if elemid[i] eq 'Mg' then nMg=i
 if elemid[i] eq 'Ca' then nCa=i
 if elemid[i] eq 'Ti' then nTi=i
 if elemid[i] eq 'Al' then nAl=i
 if elemid[i] eq 'Cr' then nCr=i
 if elemid[i] eq 'Mn' then nMn=i
 if elemid[i] eq 'Ni' then nNi=i
endfor
;;
;;;earth test
nref=nAl
earthabuN=(bulkearthppm/atomwt)/(bulkearthppm[nref]/atomwt[nref])*1e6
mantleOusedNarr=[earthabuN[nSi]*2., earthabuN[nMg], earthabuN[nAl]/2.*3., earthabuN[nCa]]
MgSiCaAl2O=total(mantleOusedNarr)/earthabuN[nO]
print, 'Earth: (Mg + 2Si + Ca + 3/2*Al) / O: ', MgSiCaAl2O

solarabuN=10^(protosundex - protosundex[nref])*1e6

earth2sunabuN=earthabuN/solarabuN
earth2sundex=alog10(earth2sunabuN)
earth2sundex_Mg2Si=earth2sundex[nMg] - earth2sundex[nSi]
earth2sundex_Fe2Si=earth2sundex[nFe] - earth2sundex[nSi]
earth2sundex_C2O=earth2sundex[nC] - earth2sundex[nO]
print, "earth to sun key ratios (in dex): [Mg/Si], [Fe/Si], [C/O]:"
print, earth2sundex_Mg2Si, earth2sundex_Fe2Si, earth2sundex_C2O
;stop

;;;;>>>>extract modelling results
;PRINT, 'CALCULATING AND EXTRACTING RESULTS FOR TARGETS...'
DPelemshow=['O', 'Mg', 'Si', 'Fe', 'Ca', 'Al', 'Na', 'Ni', 'S', 'C']

starlable='protosun'
solardexnorm=protosundex - protosundex
solardexnormerr=protosundexerrup
Result_proto=exomodel_revised(solardexnorm, solardexnormerr, starlable=starlable)
;stop

starlable='K10' ;; All good
Result_k10=exomodel_revised(K10HETdex, K10HETdexerr, starlable=starlable)
;stop

starlable='K20' ;; Star's O is starved
Result_k20=exomodel_revised(K20dex, K20dexerr, starlable=starlable)
;stop

starlable='K21' ;; NO S
Result_k21=exomodel_revised(K21dex, K21dexerr, starlable=starlable)
;stop

starlable='K100' ;; All good
Result_k100=exomodel_revised(K100dex, K100dexerr, starlable=starlable)

;;;;>>>>>>USE THE SAME FORMAT AS ABOVE TO ADD OTHER STARS OF YOUR INTEREST
;;;;---YOUR STARS---;;;;;









;stop;---------------------------
TOC
print, 'excecution done'
stop
END
