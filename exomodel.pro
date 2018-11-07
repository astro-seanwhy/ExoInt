Function exomodel, stardex, stardexerr, nodevostar=nodevostar, starlable=starlable, xplot=xplot ; sg3level=sg3level, 
;;;;DO NOT MODIFY THIS CODE, WITH THE PERMISION OF THE CREATER
;;;;A Geo-&Cosmo-model [CODE-2] for estimating rocky exoplanetary
;;composition and structures: mantle and core composition, core mass fraction, core radius.
;;Created and maintained by Haiyang S. Wang
;;Correspondance: haiyang.wang@anu.edu.au 
;;Citation: Wang, H. S., Liu, F., Ireland, T., Brasser, R., Yong, D., and Lineweaver, C. H. 2018. Enhanced constraints on the interior composition and structure of terrestrial exoplanets. MNRAS, in press. doi.org/10.1093/mnras/sty2749
 
TIC
cgCleanup ;;This procedure cleans-up and/or destroys any open graphics or widget windows on the display.
if keyword_set(xplot) then set_plot, 'x'
outpath='outputs/'
respath='results/'
if keyword_set(nodevostar) then begin
outpath='outputs_2/'
respath='results_2/'
endif
;;;>>>>>>>KEYS SET UP
COTcCorrectlable='YES'
solarHgCorrectlable='YES'

secondaryOxidizeslable='NO' ;; if do secondary oxidizes: TiO2, Cr2O3, MnO, K2O, and P2O5. 
sunearthlable='NO' ;; if calculate the key elemental ratios (by number) of the sun and the earth.
histoforkeyRdexlable='YES' ;; histogauss cal for computing error of key R in dex 
symetriclable='YES' ;;set all lower error of inputs (model and stellar abu) equal to upper error bar 
oxidesMClable='G' ;; DISTRIBUTION OF MC FOR CALING OXIDES: G-GUASSIAN ; U-UNIFORM
outputref2protosun='NO'
 
Nelems=83
nan=!values.f_nan
Nsim=2e4 ;2e4 ;; the sim should be no less than 1e3 to make reasonably results; ideally, 

if keyword_set(starlable) then print, 'Run for '+starlable

;;;>>>INPUT NECESSARY SUPPORTING DATA
readcol, 'data/atomwttc_new.txt', F='I,A,F,F', atomid, elemid, atomwt, elemtc
readcol, 'data/protosunppmwhy.txt', F='x,d,d,d', solarppm, solarppmerrup, solarppmerrdn
readcol, 'data/PEppmwhy.txt', F='x,d,d,d', bulkearthppm, bulkearthppmerrup, bulkearthppmerrdn
readcol, 'data/Asplund09dex.txt', F='x,x,f,f', A09solardex, A09solardexerr, skipline=3
readcol, 'data/ASG15dex.txt', F='x,x,f,f', ASGsolardex, ASGsolardexerr, skipline=3
readcol, 'data/protosunwhy.txt', F='x,f,f,f', protosundex, protosundexerrup, protosundexerrdn
readcol, 'data/devolmodel_lin_pattern.txt', F='x,f,f,f', ymodelN, ymodelNsdup, ymodelNsddn
readcol, 'data/devolmodel_log_pattern.txt',  F='x,f,f,f', ymodellog, ymodeldexsdup, ymodeldexsddn
;; ENDELSE
ymodeldexsdmax=max([[ymodeldexsdup], [ymodeldexsddn]], dimension=2, /nan)
;stop

;;;process the stellar abundance with -11 replaced by nan
A09solardex(where(A09solardex eq -11))=nan
A09solardexerr(where(A09solardexerr eq -11))=nan
ASGsolardex(where(ASGsolardex eq -11))=nan
ASGsolardexerr(where(ASGsolardexerr eq -11))=nan
solardex=A09solardex
solardexerr=A09solardexerr

stardex(where(stardex eq -11, /null))=nan ;;;Use /NULL, if no elements match then array will not be modified: OR Use Count to get the number of nonzero elements
stardexerr(where(stardexerr eq -11, /null))=nan ;; Only subscript the array if it is safe
;stop

DPelemshow=['O', 'Mg', 'Si', 'Fe', 'Ca', 'Al', 'Na', 'Ni', 'S', 'C']
Nelemshow=N_ELEMENTS(DPelemshow)
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
 if elemid[i] eq 'K' then nK=i
 if elemid[i] eq 'P' then nP=i
endfor
if solarHgCorrectlable eq 'YES' then begin
 solarppm[nHg]=solarppm[nHg]*13.
 solarppmerrup[nHg]=solarppmerrup[nHg]*13.
 solarppmerrdn[nHg]=solarppmerrdn[nHg]*13.
endif
if COTcCorrectlable eq 'YES' then begin
   elemtc[nC]=305
   elemtc[nO]=875
   TcCerrup=73 ;60
   TcCerrdn=135
   TcOerrup=45 ;40
   TcOerrdn=45 ;40
endif 
;;locate elements WITH available stellar abundance
avaielemid=elemid
avaielemTc=elemtc ; make_array(Nelems, value=nan)
avaielemTc(where(finite(stardex) eq 0))=nan
avaielemid(where(finite(stardex) eq 0))=' '
IF where(avaielemid eq 'O') eq -1 OR where(avaielemid eq 'Mg') eq -1 OR where(avaielemid eq 'Si') eq -1 OR where(avaielemid eq 'Fe') eq -1 THEN BEGIN
print, 'ERROR: You must have measurements for O, Mg, Si, and Fe, AT LEAST, for estimating a reasonable planetary composition !'
stop
ENDIF

IF where(avaielemid eq 'C') eq -1 THEN print, 'WARNING: No C/O can be outputed when missing stellar C abundance'
stardex_C2O= stardex[nC] -stardex[nO]
stardex_Mg2Si= stardex[nMg] - stardex[nSi]
;stardex_Mantlecity= (stardex[nMg] + 2*stardex[nSi]

 ;;;>>>>>>>>>>>calculate the atomic abundance of the available elements, and thus
;;;the mass fractions of the probably formed compounds. 
;;first, convert solar dex to atomic abundance [N(H)=10^12]
;;SOLAR
Nmc=1e5 ;;fixed
nref=nAl
solarabuN=10^(solardex - solardex[nref]) * 100. ;to Al=100 ;1e6  ;; to Al=10^6
solarabuNerr=(10^solardexerr-1)*solarabuN ;;conservative error bars
solarabuNerrdn=(1-10^(-solardexerr))*solarabuN
if symetriclable eq 'YES' then solarabuNerrdn=solarabuNerr
solarabuNarrup=fltarr(Nmc/2., Nelems)
solarabuNarrdn=fltarr(Nmc/2., Nelems)
for i=0, Nelems-1 do begin 
solarabuNarrup[*,i]=solarabuN[i] + abs(randomn(seed, Nmc/2.)*solarabuNerr[i])
solarabuNarrdn[*,i]=solarabuN[i] - abs(randomn(seed, Nmc/2.)*solarabuNerrdn[i])
endfor
solarabuNarr=[solarabuNarrup, solarabuNarrdn]
;;Earth
earthabuN=(bulkearthppm/atomwt)/(bulkearthppm[nref]/atomwt[nref]) * 100. ;1e6
;bulkearthppmerr=bulkearthppmerrup ;; for C, O, Mg, Si, Fe, Ca, Al, Na, and S, error bars of Wang+2018 are symetric
earthabuNerr=(bulkearthppmerrup/bulkearthppm)*earthabuN ;
earthabuNerrdn=(bulkearthppmerrdn/bulkearthppm)*earthabuN 
earthabuNarrup=fltarr(Nmc/2., Nelems)
earthabuNarrdn=fltarr(Nmc/2., Nelems)
for i=0, Nelems-1 do begin
earthabuNarrup[*,i]=earthabuN[i] + abs(randomn(seed, Nmc/2.)*earthabuNerr[i])
earthabuNarrdn[*,i]=earthabuN[i] - abs(randomn(seed, Nmc/2.)*earthabuNerrdn[i])
endfor
earthabuNarr=[earthabuNarrup, earthabuNarrdn]

;;STAR
refsolardex=solardex
if starlable eq 'protosun' then refsolardex=protosundex
starabudex=stardex+refsolardex  ;;;differential dex;; irrelevant to solar abundances
starabudexerr=stardexerr
;nref=nAl ;Fe ;nAl
if finite(stardex[nAl]) eq 0 AND finite(stardex[nCa]) eq 1 then nref=nCa
if finite(stardex[nAl]) eq 0 AND finite(stardex[nCa]) eq 0 AND finite(stardex[nFe]) eq 1 then nref=nFe
starabuN=10^(starabudex - starabudex[nref]) * 100 ;1e6 
starabuNerr=(10^starabudexerr - 1)*starabuN
starabuNerrdn=(1-10^(-starabudexerr))*starabuN
if starlable eq 'protosun' then starabuNerrdn=(1-10^(-protosundexerrdn))*starabuN
if symetriclable eq 'YES' then starabuNerrdn=starabuNerr

;;Matrix
;starabuNarr=randomn(seed, Nmc)#starabuNerr + replicate(1, Nmc)#starabuN
starabuNarrup=replicate(1, Nmc/2.)#starabuN + abs(randomn(seed, Nmc/2.)#starabuNerr)
starabuNarrdn=replicate(1, Nmc/2.)#starabuN - abs(randomn(seed, Nmc/2.)#starabuNerrdn)
for i=0, Nelems-1 do begin
starabuNarrup[*,i]=starabuN[i] + abs(randomn(seed, Nmc/2.)*starabuNerr[i])
starabuNarrdn[*,i]=starabuN[i] - abs(randomn(seed, Nmc/2.)*starabuNerrdn[i])
endfor
starabuNarr=[starabuNarrup, starabuNarrdn]


;;;>>>INDEPENDENT CALCULATIONS BY ERROR PROPAGATION 
;;STAR's planet
planetabuN=starabuN*ymodelN  ;;;starabuN and ymodelN are totally independent sources
;goto, fromdexcal

planetabuNerrup=sqrt(starabuNerr^2.*ymodelN^2. + ymodelNsdup^2.*starabuN^2.) ;;up-up
planetabuNerrdn=sqrt(starabuNerrdn^2.*ymodelN^2. + ymodelNsddn^2.*starabuN^2.) ;;dn-dn
planetabuNarrup=make_array(Nmc/2., Nelems, value=nan)
planetabuNarrdn=make_array(Nmc/2., Nelems, value=nan)
for i=0, Nelems-1 do begin
planetabuNarrup[*,i]=planetabuN[i] + abs(randomn(seed, Nmc/2.)*planetabuNerrup[i])
planetabuNarrdn[*,i]=planetabuN[i] - abs(randomn(seed, Nmc/2.)*planetabuNerrdn[i])
endfor
planetabuNarr=[planetabuNarrup, planetabuNarrdn]

print, 'output planet abundance (by number normalized to Al=1) to file'
writecol, outpath+starlable+'exoplabuN.txt', atomid, elemid, planetabuN/100., planetabuNerrup/100., planetabuNerrdn/100., fmt='(i2, A, F,F,F)'

;;;;>>>>>>>calculate key elemental ratios for a rocky planet around a star
;;with known stellar abundances for these elements
;;Mg/Si, Fe/Si, C/O, Na/Al
;;the K10HETdex is in [X/H] (dex), convert the normalization to Al
;;(with error for elements unchanged.)
;;Elemental abundance, [X/Al]_p (=log ((X/Al)_p/(X/Al)_sun) = log(f) +
;;[X/Al]_star ;; if upper and lower limit, I should calculate the
;;limits directly and then subtract with the best fit to get the upper
;;and lower error bars of results. The reason is that the upper and
;;lower limits are one of the possibilities, rather than the normal
;;uncertainty (it is deemed certain in one direction), thus rather than
;;using the upper error bar


if finite(starabudex[nAl]) eq 0 and finite(starabudex[nCa]) eq 1 then nref=nCa
if finite(starabudex[nAl]) eq 0 and finite(starabudex[nCa]) eq 0 AND finite(starabudex[nFe]) eq 1 then nref=nFe
starabudex2Al=starabudex - starabudex[nref]
starabudex2Alerr=starabudexerr

ymodeldex=ymodellog
planetabudex2Al=ymodeldex + starabudex2Al
planetabudex2Alerr0=sqrt(ymodeldexsdmax^2. + starabudex2Alerr^2.)

;;ARR
ymodeldexarrup=fltarr(Nmc/2., Nelems)
ymodeldexarrdn=fltarr(Nmc/2., Nelems)
for i=0, Nelems-1 do begin 
ymodeldexarrup[*,i]=ymodellog[i] + abs(randomn(seed, Nmc/2.)*ymodeldexsdup[i])
ymodeldexarrdn[*,i]=ymodellog[i] - abs(randomn(seed, Nmc/2.)*ymodeldexsddn[i]) 
endfor
ymodeldexarr=[ymodeldexarrup, ymodeldexarrdn]

stardexarr=randomn(seed, Nmc)#starabudex2Alerr + replicate(1, Nmc)#starabudex2Al
for i=0, Nelems-1 do begin
stardexarr[*,i]=randomn(seed, Nmc)*starabudex2Alerr[i] + starabudex2Al[i]
endfor
planetdexarr=stardexarr + ymodeldexarr ;ymodellogarrc

if keyword_set(xplot) then window, /free, title='Histograms of Planet Abundance in Dex'
planetabudex2Alerr=make_array(Nelems, value=nan)
for i=0, Nelems-1 do begin
if  where(DPelemshow eq elemid[i]) NE -1 AND finite(starabudex2Al[i]) eq 1 then begin
histogauss, planetdexarr[*,i], A, /noplot
planetabudex2Alerr[i]=A[2] ;; A[0]-height, A[1]-mean, A[2]-weight of histo (standard deviation)
;stop
endif
endfor


planetdex=planetabudex2Al - (solardex - solardex[nref])  ;;2Al2sun, namely [X/Ca]
if outputref2protosun eq 'YES' then planetdex=planetabudex2Al - (protosundex - protosundex[nref])
planetdexerr=planetabudex2Alerr0
print, 'output planet abundance in dex to file'
writecol, outpath+starlable+'exopldex.txt', atomid, elemid, planetdex, planetdexerr, fmt='(i2, A, F,F)'
;STOP
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print, 'Cal key elemental ratios in both linear and logarithm...'
;;;;;;
planetC2O=planetabuN[nC] / planetabuN[nO]
print, starlable+' planetary C/O: ', planetC2O  ;;According to Unterborn+2014, C would be oxidized after Fe, and thus it is reansonable to consider C as the last. ; But this must be the case when C/O < 0.8; otherwise, diamond or carbide (e.g. FeC, CaC) mantle would form.  
Mantleicity0=(planetabuN[nMg] + 2.*(planetabuN[nSi])) / (planetabuN[nO])
print, starlable+' planetary mean: (Mg + 2Si) / O: ', Mantleicity0
Mantleicity1=total([planetabuN[nMg], 2.*(planetabuN[nSi]),  planetabuN[nCa], 3/2.*planetabuN[nAl]], /nan) / (planetabuN[nO])
print, starlable+' planetary mean: (Mg + 2Si + Ca + 3/2Al) / O: ', Mantleicity1


;stop
;;;>>>>cal the abundance ratios of the bodies themselves
;;solar
;set_plot, 'x'
ratiolable='star'
bodyabuN=starabuN
bodyabuNarr=starabuNarr
calratio:
;Nmc=1e5
C2O=bodyabuN[nC]/bodyabuN[nO]
C2ONarr=bodyabuNarr[*,nC]/bodyabuNarr[*,nO]

if keyword_set(xplot) then begin
window, 1, title=starlable+ratiolable+' C/O'
histogauss, C2ONarr, A
endif else histogauss, C2ONarr, A, /noplot
C2Oerr=A[2]
;stop

Mg2Si=bodyabuN[nMg]/bodyabuN[nSi]
Mg2SiNarr=bodyabuNarr[*,nMg]/bodyabuNarr[*,nSi]
if keyword_set(xplot) then begin
window, 2, title=starlable+ratiolable+' Mg/Si'
histogauss, Mg2SiNarr, A 
endif else histogauss, Mg2SiNarr, A ,/noplot
Mg2Sierr=A[2]

maxFeO= bodyabuN[nO] - (bodyabuN[nMg] + bodyabuN[nSi]*2.)
mantlecity=maxFeO/bodyabuN[nFe]
maxFeONarr=bodyabuNarr[*,nO] - (bodyabuNarr[*,nMg] + bodyabuNarr[*,nSi]*2.)
mantlecityNarr=maxFeONarr/bodyabuNarr[*,nFe]
if keyword_set(xplot) then begin 
window, 3, title=starlable+ratiolable+' (O-Mg-2Si)/Fe'
histogauss, mantlecityNarr, A 
endif else histogauss, mantlecityNarr, A ,/noplot
mantlecityerr=A[2]
;STOP
Fe2Ni=bodyabuN[nFe]/bodyabuN[nNi]
Fe2NiNarr=bodyabuNarr[*,nFe]/bodyabuNarr[*,nNi]
Fe2Nierr=nan
if finite(bodyabuN[nNi]) eq 1 then begin
if keyword_set(xplot) then begin
window, 4, title=starlable+ratiolable+' Fe/Ni'
histogauss, Fe2NiNarr, A
endif else histogauss, Fe2NiNarr, A,/noplot
Fe2Nierr=A[2]
endif

case ratiolable of
'star': begin
starN_keyR=[C2O, Mg2Si, Mantlecity, Fe2Ni]
starN_keyRerr=[C2Oerr, Mg2Sierr, Mantlecityerr, Fe2Nierr]

ratiolable='planet'
bodyabuN=planetabuN
bodyabuNarr=planetabuNarr
goto, calratio
end
'planet': begin
planetN_keyR=[C2O, Mg2Si, Mantlecity, Fe2Ni]
planetN_keyRerr=[C2Oerr, Mg2Sierr, Mantlecityerr, Fe2Nierr]

ratiolable='solar'
bodyabuN=solarabuN
bodyabuNarr=solarabuNarr
if sunearthlable EQ 'YES' then goto, calratio
end
'solar': begin
solarN_keyR=[C2O, Mg2Si, Mantlecity, Fe2Ni]
solarN_keyRerr=[C2Oerr, Mg2Sierr, Mantlecityerr, Fe2Nierr]

ratiolable='earth'
bodyabuN=earthabuN
bodyabuNarr=earthabuNarr
goto, calratio
end
'earth': begin
earthN_keyR=[C2O, Mg2Si, Mantlecity, Fe2Ni]
earthN_keyRerr=[C2Oerr, Mg2Sierr, Mantlecityerr, Fe2Nierr]
end
;endratioN:
endcase
keyRlable=['C/O', 'Mg/Si', '(O-Mg-2Si)/Fe', 'Fe/Ni']
print, starlable+'starkeyR: '+keyRlable, starN_keyR
print, starlable+'starkeyRerr: '+keyRlable, starN_keyRerr
print, starlable+'planetkeyR: '+keyRlable, planetN_keyR
print, starlable+'planetkeyRerr: '+keyRlable, planetN_keyRerr
;
writecol, outpath+starlable+'NkeyR.txt',  keyRlable, planetN_keyR, planetN_keyRerr, starN_keyR, starN_keyRerr, fmt='(A, F,F,F,F)'


;;;>>>>>>>>RATIOS IN DEX<<<<<<<<<<<<<<<<<<<<<<
;;Key elemental ratios, e.g. [Mg/Si]_p = [Mg/Al]_p - [Si/Al]_p;; 
keyRlable='Mg2Si'
nelem1=nMg
nelem2=nSi
stardexc=stardex
if outputref2protosun eq 'YES' then stardexc= starabudex - protosundex
if starlable eq 'solar' then planetdexerr=ymodeldexsdmax

starMg2Si=stardexc[nelem1]-stardexc[nelem2]
starMg2Sierr=sqrt(stardexerr[nelem1]^2. + stardexerr[nelem2]^2.)

planetMg2Si=planetdex[nelem1]-planetdex[nelem2]
planetMg2Sierr=sqrt(planetdexerr[nelem1]^2. + planetdexerr[nelem2]^2.)
;

keyRlable='C2O'
nelem1=nC
nelem2=nO
starC2O=stardexc[nelem1]-stardexc[nelem2]
starC2Oerr=sqrt(stardexerr[nelem1]^2. + stardexerr[nelem2]^2.)
;

planetC2O=planetdex[nelem1]-planetdex[nelem2]
planetC2Oerr=sqrt(planetdexerr[nelem1]^2. + planetdexerr[nelem2]^2.)
;

keyRlable='Fe2Si'
nelem1=nFe
nelem2=nSi
starFe2Si=stardexc[nelem1] - stardexc[nelem2] ;- stardex[nO]
starFe2Sierr=sqrt(stardexerr[nelem1]^2. + stardexerr[nelem2]^2.) ; + stardexerr[nO]^2.)
;

planetFe2Si=(planetdex[nelem1] - planetdex[nelem2]) ;- planetdex[nFe] ; - planetdex[nO]
planetFe2Sierr=sqrt(planetdexerr[nelem1]^2. + planetdexerr[nelem2]^2.) ; + planetdexerr[nFe]^2.) ; + planetdexerr[nO]^2.)
;

if starlable eq 'solar' then planetdexarr = ymodeldexarr

IF histoforkeyRdexlable eq 'YES' THEN BEGIN
starMg2SiNarr=stardexarr[*,nMg] - stardexarr[*,nSi]
histogauss, starMg2SiNarr, A, /noplot
starMg2Sierr=A[2]
;
planetMg2SiNarr=planetdexarr[*,nMg] - planetdexarr[*,nSi]
if keyword_set(xplot) then begin
window, /free, title=starlable+ 'b '+keyRlable+' in DEX'
histogauss, planetMg2SiNarr, A
endif else histogauss, planetMg2SiNarr, A, /noplot
planetMg2Sierr=A[2]

starC2ONarr=stardexarr[*,nC] - stardexarr[*,nO]
histogauss, starC2ONarr, A, /noplot
starC2Oerr=A[2]
;
planetC2ONarr=planetdexarr[*,nC] - planetdexarr[*,nO]
if keyword_set(xplot) then begin
window, /free, title=starlable+ 'b '+keyRlable+' in DEX'
histogauss, planetC2ONarr, A
endif else histogauss, planetC2ONarr, A ,/noplot
planetC2Oerr=A[2]

starFe2SiNarr=stardexarr[*,nFe] - stardexarr[*,nSi]
histogauss, starFe2SiNarr, A, /noplot
starFe2Sierr=A[2]
;
planetFe2SiNarr=planetdexarr[*,nFe] - planetdexarr[*,nSi]
if keyword_set(xplot) then begin
window, /free, title=starlable+ 'b '+keyRlable+' in DEX'
histogauss, planetFe2SiNarr, A
endif else histogauss, planetFe2SiNarr, A,/noplot
planetFe2Sierr=A[2]
ENDIF

KeyRdexlable=['[Fe/Si]', '[Mg/Si]', '[C/O]']
planetkeyR=[planetFe2Si, planetMg2Si, planetC2O];,  planetM2O]
planetkeyRerr=[planetFe2Sierr, planetMg2Sierr, planetC2Oerr];,  planetM2Oerr]
starkeyR=[starFe2Si, starMg2Si, starC2O];, starM2O]
starkeyRerr=[starFe2Sierr, starMg2Sierr, starC2Oerr];, starM2Oerr]

writecol, outpath+starlable+'keyRdex.txt',  keyRdexlable, planetkeyR, planetkeyRerr, starkeyR, starkeyRerr, fmt='(A, F,F,F,F)'

;;;;Check the sufficiency of O per every group of compounds
;;;;;;;;;;;;;;;;

;;;;;>>>>set the situation with no depletion
val_id=where(finite(starabuN) eq 1)
IF Keyword_set(nodevostar) THEN BEGIN
planetabuN=starabuN
planetabuNarr=starabuNarr
print, 'No Devol used for stochiometric cal of compounds!'
ENDIF

Print, 'Simul Using the best Mean Values of Abundances...'

;Results_m=chemsysmodel_upgrade(planetabuN, /outinfull)
Results_m=chemsysmodel(planetabuN, /outinfull)
print, starlable+' mean results of mantle and core compositions and core mass fraction:'
print, 'mantle compounds name: ', Results_m.mantlecompsname
print, 'Number of groups of results:', Results_m.Nval
;;;
sizemantle=size(Results_m.mantlecompsmassfra)
Nmantlecomps=sizemantle[1] ;;14, as designed, changable
sizecore=size(Results_m.corecompsmassfra)
Ncorecomps=sizecore[1];; 3, as designed, changable

if Results_m.Nval le 1 then begin 
mantlecompsmassfra_mean=Results_m.mantlecompsmassfra
corecompsmassfra_mean=Results_m.corecompsmassfra
fcoremass_mean=Results_m.fcoremass
endif else begin 
mantlecompsmassfra_mean=make_array(Nmantlecomps, value=nan)
for i=0, Nmantlecomps-1 do mantlecompsmassfra_mean[i]=mean(Results_m.mantlecompsmassfra[i,*], /nan)
;
corecompsmassfra_mean=make_array(Ncorecomps, value=nan)
for i=0, Ncorecomps-1 do corecompsmassfra_mean[i]=mean(Results_m.corecompsmassfra[i,*], /nan)
;
fcoremass_mean=mean(Results_m.fcoremass, /nan)

endelse
print, 'Mean of BEST abundances of '+starlable
print, 'C_mantle: Mantle stoichimetric composition (BEST mean):', mantlecompsmassfra_mean
print, 'C_core: Core stoichimetric composition (BEST mean):', corecompsmassfra_mean
print, 'F_core: Core mass fraction (BEST mean):', fcoremass_mean
;;;>>>>>>>THE RESULTS ABOVE FOR VERIFICATION 
;STOP


print, 'Simul Using MC to Produce Error Bars...'
;;;>>>Generate randomn arrays for planet abundances (with uncertainties)
;Nsim=2e4; 2e4 ;2e4; 2e4
Ncut=Nsim/100 ;;less than 1/10th of the total simulation is not trusted. 

;planetabuNarrMC=planetabuNarr ;;Gaussian Distribution make_array(Nelems, Nmc, value=nan)
planetabuNarrMCup=make_array(Nsim/2., Nelems, value=nan)
planetabuNarrMCdn=make_array(Nsim/2., Nelems, value=nan)
for i=0, Nelems-1 do begin
planetabuNarrMCup[*,i]=planetabuN[i] + abs(randomn(seed, Nsim/2.)*planetabuNerrup[i])
planetabuNarrMCdn[*,i]=planetabuN[i] - abs(randomn(seed, Nsim/2.)*planetabuNerrdn[i])
endfor
planetabuNarrMC=[planetabuNarrMCup, planetabuNarrMCdn]
IF oxidesMClable eq 'U' THEN BEGIN
for i=0, Nelems-1 do begin
if finite(planetabuN[i]) eq 1 then planetabuNarrMC[*, i] = (planetabuN[i]-planetabuNerrdn[i]) + randomu(seed, Nmc)*(planetabuNerrup[i]+planetabuNerrdn[i])  ;;uniform distribution 
endfor
ENDIF

mantlecompsmassfraARR=[ ]
corecompsmassfraARR=[ ]
fcoremassARR=[ ]
;NvalARR=[ ]
Nvalsim=0 ;; valid simulation (i.e. generate at least one valid pair result for FeO, NiO, SO3 (and Fe, Ni, or S). 
Nvalarrsim=0 ;; valid simulation with more than 1 pair
Nvalarrtotal=0
for i=0, Nsim-1 do begin
print, starlable+' MC simul', i+1, ' of ', Nsim
;;;;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
;Results=chemsysmodel_upgrade(planetabuNarrMC[i, *], /outinfull) ;;;if multiple input and then concatenate, /outinfull must be used to keep the same dimensions
Results=chemsysmodel(planetabuNarrMC[i, *], /outinfull)
;;;;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mantlecompsmassfraARR=[[mantlecompsmassfraARR], [Results.mantlecompsmassfra]] ;;might be multi dimensional, then should have [ ]
corecompsmassfraARR=[[corecompsmassfraARR], [Results.corecompsmassfra]]
Ncoreval=n_elements(Results.corecompsmassfra)
fcoremassARR=[fcoremassARR, Results.fcoremass] ;; Results.fcoremass is 1 scalar or 1-dimension results, then concatenate by column
Nvalsim=Nvalsim + (Results.Nval ge 1)
Nvalarrsim=Nvalarrsim + (Results.Nval gt 1)
Nvalarrtotal=Nvalarrtotal+Results.Nval*(Results.Nval gt 1)
;NvarARR=[NvalARR, Results.Nval]
endfor
;
fullcompsname=Results.mantlecompsname
fullcorecompsname=Results.corecompsname
Nfullmantle=N_ELEMENTS(fullcompsname)
Nfullcore=N_ELEMENTS(fullcorecompsname)
nFeO=where(fullcompsname eq 'FeO')
nMetal=where(fullcompsname eq 'Metals')


;stop
;;Cal the mean and 1-sigma uncertainty of the results
mantlecompsmassfraERR=make_array(Nfullmantle, value=nan)
mantlecompsmassfraMAX=make_array(Nfullmantle, value=nan)
mantlecompsmassfraMIN=make_array(Nfullmantle, value=nan)
for i=0, Nfullmantle-1 do begin
if n_elements(where(finite(mantlecompsmassfraARR[i,*]) eq 1)) ge Ncut then begin
mantlecompsmassfraERR[i]=stddev(mantlecompsmassfraARR[i,*], /nan)
mantlecompsmassfraMAX[i]=max(mantlecompsmassfraARR[i,*], /nan)
mantlecompsmassfraMIN[i]=min(mantlecompsmassfraARR[i,*], /nan)
endif
endfor


corecompsmassfraERR=make_array(Nfullcore, value=nan)
corecompsmassfraMAX=make_array(Nfullcore, value=nan)
corecompsmassfraMIN=make_array(Nfullcore, value=nan)
for i=0, Nfullcore-1 do begin
if n_elements(where(finite(corecompsmassfraARR[i,*]) eq 1)) ge Ncut then begin
corecompsmassfraERR[i]=stddev(corecompsmassfraARR[i,*], /nan)
corecompsmassfraMAX[i]=max(corecompsmassfraARR[i,*], /nan)
corecompsmassfraMIN[i]=min(corecompsmassfraARR[i,*], /nan)
endif
endfor
;;normalize to 1

fcoremassERR=stddev(fcoremassARR, /nan)
fcoremassMAX=max(fcoremassARR, /nan)
fcoremassMIN=min(fcoremassARR, /nan)

;;;--OUTPUT RESULTS---
result={mantlecompsname:fullcompsname, mantlecompsmassfraBEST:mantlecompsmassfra_mean, mantlecompsmassfraERR:mantlecompsmassfraERR, corecompsname:fullcorecompsname, corecompsmassfraBEST:corecompsmassfra_mean, corecompsmassfraERR:corecompsmassfraERR, fcoremassBEST:fcoremass_mean, fcoremassERR:fcoremassERR, planetN_keyR:planetN_keyR, planetN_keyRerr:planetN_keyRerr, starN_keyR:starN_keyR, starN_keyRerr:starN_keyRerr, planetdex:planetdex, planetdexerr:planetdexerr, planetkeyR:planetkeyR, planetkeyRerr:planetkeyRerr, starkeyR:starkeyR, starkeyRerr:starkeyRerr}

print, 'output results to files...'
writecol, outpath+starlable+'_mantlecomp.txt', fullcompsname, mantlecompsmassfra_mean, mantlecompsmassfraERR, fmt='(A, F,F)'
writecol, outpath+starlable+'_corecomp.txt',  fullcorecompsname, corecompsmassfra_mean, corecompsmassfraERR, fmt='(A, F,F)' 
writecol, outpath+starlable+'_fcoremass.txt',  fcoremass_mean, fcoremassERR, fmt='(F,F)' ;;fcoremassmeanOK, fcoremassERROK,   

print, 'execution done for'+starlable
TOC
return, result
END
