PRO exoconstraints
;;THE ONLY PART OF THE CODES THAT YOU NEED TO RUN OR MODIFY
;;A Geo-&Cosmo-model [CODE-1] for estimating rocky exoplanetary
;;composition and structures: mantle and core composition, core mass fraction, core radius.
;;Created and maintained by Haiyang S. Wang
;;Citation: Wang, H. S., Liu, F., Ireland, T., Brasser, R., Yong, D.,
;;and Lineweaver, C. H. 2019. Enhanced constraints on the interior
;;composition and structure of terrestrial exoplanets. MNRAS 482:2222-2233. doi.org/10.1093/mnras/sty2749

TIC
cgCleanup ;;This procedure cleans-up and/or destroys any open graphics or widget windows on the display.
outpath='output/'
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
A09solardex(where(A09solardex eq -11))=nan
A09solardexerr(where(A09solardexerr eq -11))=nan
solardex=A09solardex
solardexerr=A09solardexerr
;;;;>>>>>>>>>>>>>>>>>>>>>>>>
;;;rename indexies
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


;;;;;--Example of stars----
readcol, 'data/Kepler10.txt', F='x,f,f', K10dex, k10dexerr ;; differential abundances for Kepler 10 as the example star
K10dex(where(K10dex eq -11))=nan
K10dexerr(where(K10dexerr eq -11))=nan

;;;------------------------------------------------------------------
;;;--------PLEASE PUT YOUR SAMPLE OF STARS AS FOLLOWS;;;;
;;notes: 
;;1. the abundances shall be differential or otherwise normalised
;;to any solar abundance of your choice (the default solar abundance
;;is from Asplund et al. 2009, ARAA)
;;2. please make sure that any denotes for invalid values in your file
;;would be reset as "!values.f_nan" in the code, as the example for
;;Kepler 10

;;>>replace the following line to your real sample of star(s)
;readcol, 'data/yoursample.txt', F='x,f,f', yoursampledex, yoursampledexerr




;;;;---------END OF SAMPLE OF STARS-----------------
;;;;;---------------------------------------------------------

;;;;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
;;;;>>>>>MODELLING OF INTERIORS: MANTLE AND CORE COMPOSITION, AS WELL AS CORE MASS FRACTION

;;;>>>>EXAMPLE
starlable='K10'
Result_k10=exomodel(K10dex, K10dexerr, starlable=starlable)
;;;;;--PARAMETERS IN THE RESULTS---
;result={mantlecompsname:fullcompsname, mantlecompsmassfra:mantlecompsmassfra_mean, mantlecompsmassfraERR:mantlecompsmassfraERR, corecompsname:fullcorecompsname, corecompsmassfra:corecompsmassfra_mean, corecompsmassfraERR:corecompsmassfraERR, fcoremass:fcoremass_mean, fcoremassERR:fcoremassERR, planetN_keyR:planetN_keyR, planetN_keyRerr:planetN_keyRerr, starN_keyR:starN_keyR, starN_keyRerr:starN_keyRerr, planetdex:planetdex, planetdexerr:planetdexerr, planetkeyR:planetkeyR, planetkeyRerr:planetkeyRerr, starkeyR:starkeyR, starkeyRerr:starkeyRerr}


print, 'output results to files...'
writecol, outpath+starlable+'_mantlecomp.txt', Result_k10.mantlecompsname, Result_k10.mantlecompsmassfra, Result_k10.mantlecompsmassfraERR, fmt='(A, F,F)'
writecol, outpath+starlable+'_corecomp.txt',  Result_k10.corecompsname, Result_k10.corecompsmassfra, Result_k10.corecompsmassfraERR, fmt='(A, F,F)' 
writecol, outpath+starlable+'_fcoremass.txt',  Result_k10.fcoremass, Result_k10.fcoremassERR, fmt='(F,F)'  
writecol, outpath+starlable+'_exopldex.txt', atomid, elemid, Result_k10.planetdex, Result_k10.planetdexerr, fmt='(i2, A, F,F)' ;; EXOPLANET ABUNDANCES (DIFFERENTIAL)
writecol, outpath+starlable+'_keyRdex.txt',  Result_k10.keyRdexlable, Result_k10.planetkeyR, Result_k10.planetkeyRerr, Result_k10.starkeyR, Result_k10.starkeyRerr, fmt='(A, F,F,F,F)' ;;;KEY ELEMENTAL RATIOS - '[Fe/Si]', '[Mg/Si]', '[C/O]' - FOR THE PLANET AND HOST STAR

;;;>>>>YOUR SAMPLE OF STARS
;;;replace the following lines to your own cases

;starlable='sample' ;;replace 'sample' to the name of the sample star(s) of your choice 
;Result_sample=exomodel(yoursampledex, yoursampledexerr, starlable=starlable)

;print, 'output results to files...'
;writecol, outpath+starlable+'_mantlecomp.txt', Result_sample.mantlecompsname, Result_sample.mantlecompsmassfra, Result_sample.mantlecompsmassfraERR, fmt='(A, F,F)'
;writecol, outpath+starlable+'_corecomp.txt',  Result_sample.corecompsname, Result_sample.corecompsmassfra, Result_sample.corecompsmassfraERR, fmt='(A, F,F)' 
;writecol, outpath+starlable+'_fcoremass.txt',  Result_sample.fcoremass, Result_sample.fcoremassERR, fmt='(F,F)'  
;writecol, outpath+starlable+'_exopldex.txt', atomid, elemid, Result_sample.planetdex, Result_sample.planetdexerr, fmt='(i2, A, F,F)' ;; EXOPLANET ABUNDANCES (DIFFERENTIAL)
;writecol, outpath+starlable+'_keyRdex.txt',  Result_sample.keyRdexlable, Result_sample.planetkeyR, Result_sample.planetkeyRerr, Result_sample.starkeyR, Result_sample.starkeyRerr, fmt='(A, F,F,F,F)' ;;;KEY ELEMENTAL RATIOS - '[Fe/Si]', '[Mg/Si]', '[C/O]' - FOR THE PLANET AND HOST STAR

;;;;;>>>>>>>>>>>>>>>>>>>>>>END OF MODELLING>>>>>>>>>>>>>>>>>>>>>>>>


;stop;---------------------------
TOC
print, 'excecution done'
stop
END
