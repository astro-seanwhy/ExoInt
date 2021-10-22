Function exomodel, stardex, stardexerr, stardexerrdn, starlabel=starlabel, quick=quick, xplot=xplot, withreferr=withreferr, F3sigerr=F3sig, F5sigerr=F5sig, upperlmt=upperlmt, lowerlmt=lowerlmt, residual_add=residual_add 
;;;;DO NOT MODIFY THIS CODE (EXCEPT CHOOSING DIFFERENT KEYWORDS), WITHOUT THE PERMISION OF THE CREATER
;;;;A Geo-&Cosmo-model [CODE-2] for estimating rocky exoplanetary
;;composition and structures: mantle and core composition, core mass fraction, core radius.
;;Created and maintained by Haiyang S. Wang
;;Correspondance: haiyang.wang@anu.edu.au 
;;Citation: Wang, H. S., Liu, F., Ireland, T., Brasser, R., Yong, D.,
;;and Lineweaver, C. H. 2019. Enhanced constraints on the interior composition and structure of terrestrial exoplanets. MNRAS 482:2222-2233. doi.org/10.1093/mnras/sty2749
;;SET DEFAULT VALUE FOR KEYWORDS
setdefaultvalue, gcelabel, 'NO'  ;;NO GCE SHOULD BE DONE FOR THE PURPOSE HERE: I.E. EXOPLANETARY COMPOSITION

TIC

if keyword_set(xplot) then begin
cgCleanup ;;This procedure cleans-up and/or destroys any open graphics or widget windows on the display.
set_plot, 'x'
endif

outpath='outputs/'
if gcelabel eq 'YES' then outpath='outputs_gce/'

if keyword_set(starlabel) then print, 'Run for '+starlabel
;othersign=''


;;;>>>>>>>KEYS SET UP

sunearthlable='YES' ;; if calculate the key elemental ratios (by number) of the sun and the earth.
histoforkeyRdexlable='NO' ;; histogauss cal for computing error of key R in dex 
symetriclable='NO' ;;set all lower error of inputs (model and stellar abu) equal to upper error bar 
oxidesMClable='G' ;; DISTRIBUTION OF MC FOR CALING OXIDES: G-GUASSIAN ; U-UNIFORM
;outputref2protosun='NO'
devolfactorplot='NO' ;; Interatively plot the scaled devol. factors or not
scalemethod=1    ;;0 - ymodel/fscaleR ; 1 - alphac=(alpha - fscaledex); betac=(beta + fscaledex*alog10(TD)); 2 - ymodel(nO)/(O/(Mg+2Si+Fe))
solarreflabel='proto' ;'A09', 'A20', or 'proto'
constrainlabel='YES'

Nelems=83
nan=!values.f_nan
Nsim=2e4 ;;; the final sim should be no less than 1e3 to make reasonable results
if keyword_set(quick) then Nsim=1e2 ; for testing/debugging or calculating the key elemental ratios only (i.e. without interior composition modeling). 

;;;>>>INPUT NECESSARY SUPPORTING DATA
readcol, 'data/atomwttc_new.txt', F='I,A,F,F', atomid, elemid, atomwt, elemtc
;readcol, 'data/protosunppmwhy.txt', F='x,d,d,d', solarppm, solarppmerrup, solarppmerrdn
readcol, 'data/PEppmwhy.txt', F='x,d,d,d', bulkearthppm, bulkearthppmerrup, bulkearthppmerrdn
readcol, 'data/Asplund09dex.txt', F='x,x,f,f', A09solardex, A09solardexerr, skipline=3
readcol, 'data/Asplund2020dex.txt', F='x,x,f,f', A20solardex, A20solardexerr, skipline=3
readcol, 'data/protosunwhy.txt', F='x,f,f,f', protosundex, protosundexerrup, protosundexerrdn
;if protosunMgwerr eq 'YES' then readcol, 'data/protosunwhy_mg_werr.txt', F='x,f,f,f', protosundex, protosundexerrup, protosundexerrdn
;readcol, 'data/devolmodel_lin_pub.txt', F='x,f,f,f', ymodelN, ymodelNsdup, ymodelNsddn
;readcol, 'data/devolmodel_log_pub.txt',  F='x,f,f,f', ymodellog, ymodeldexsdup, ymodeldexsddn
readcol, 'data/devolmodel_lin_pattern.txt', F='x,f,f,f,f,f', ymodelN, ymodelNsdup, ymodelNsddn, ymodelN3lmtup, ymodelN3lmtdn, ymodelN5lmtup, ymodelN5lmtdn
readcol, 'data/devolmodel_log_pattern.txt',  F='x,f,f,f,f,f', ymodellog, ymodeldexsdup, ymodeldexsddn, ymodeldex3lmtup, ymodeldex3lmtdn, ymodeldex5lmtup, ymodeldex5lmtdn
readcol, 'data/Earth_residual_toAl.txt', F = 'x, F', residual

;; ENDELSE
;stop
;;;process the datasets with -11 replaced by nan
A09solardex(where(A09solardex eq -11))=nan
A09solardexerr(where(A09solardexerr eq -11))=nan
A20solardex(where(A20solardex eq -11))=nan
A20solardexerr(where(A20solardexerr eq -11))=nan

stardex(where(stardex eq -11, /null))=nan ;;;Use /NULL, if no elements match then array will not be modified: OR Use Count to get the number of nonzero elements
stardexerr(where(stardexerr eq -11, /null))=nan ;; Only subscript the array if it is safe
stardexerrdn(where(stardexerrdn eq -11, /null))=nan


;;;;---reference solar abundances

refsolardex=A09solardex
refsolardexerr=A09solardexerr
refsolardexerrdn=A09solardexerr
if solarreflabel eq 'A20' then begin
refsolardex=A20solardex
refsolardexerr=A20solardexerr
refsolardexerrdn=A20solardexerr
endif
if solarreflabel eq 'proto' then begin
refsolardex=protosundex
refsolardexerr=protosundexerrup
refsolardexerrdn=protosundexerrdn
endif

;;>>when input star is the sun or proto-sun, then fix the reference star as it is.
if starlabel eq 'protosun' then begin
   refsolardex=protosundex
   refsolardexerr=protosundexerrup
   refsolardexerrdn=protosundexerrdn
   solarreflabel='proto'
endif
if starlabel eq 'solar09' then begin
   refsolardex=A09solardex
   refsolardexerr=A09solardexerr
   refsolardexerrdn=A09solardexerr
   solarreflabel='A09'
endif
if starlabel eq 'solar20' then begin
   refsolardex=A20solardex
   refsolardexerr=A20solardexerr
   refsolardexerrdn=A20solardexerr
   solarreflabel='A20'
endif

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
 if elemid[i] eq 'Eu' then nEu=i
endfor

if finite(stardex(nO)) eq 0 OR finite(stardex(nMg)) eq 0 OR finite(stardex(nSi)) eq 0 OR finite(stardex(nFe)) eq 0 then begin 
print, 'ERROR: You must have measurements for O, Mg, Si, and Fe, AT LEAST, for estimating a reasonable rocky planetary composition !'
stop
endif

IF finite(stardex(nC)) eq 0 THEN print, 'WARNING: No C/O can be outputed when missing stellar C abundance'
stardex_C2O= stardex[nC] -stardex[nO]
stardex_Mg2Si= stardex[nMg] - stardex[nSi]


 ;;;>>>>>>>>>>>calculate the atomic abundance of the available elements, and thus
;;;the mass fractions of the probably formed compounds. 
;;first, convert solar dex to atomic abundance [N(H)=10^12]
;;SOLAR
Nmc=1e5 ;;fixed
nref=nAl
solarabuN=10^(refsolardex - refsolardex[nref]) * 100 ;100. ;to Al=100 ;1e6  ;; to Al=10^6
solarabuNerrup=(10^refsolardexerr-1)*solarabuN ;;conservative error bars
solarabuNerrdn=(1-10^(-refsolardexerrdn))*solarabuN

solarabuNerr=solarabuNerrup
errdiff = solarabuNerrup - solarabuNerrdn
solarabuNerr(where(errdiff lt 0)) = solarabuNerrdn(where(errdiff lt 0)) ;;conservative error bars


solarabuNarrup=fltarr(Nmc/2., Nelems)
solarabuNarrdn=fltarr(Nmc/2., Nelems)
solarabuNarr=fltarr(Nmc, Nelems)
for i=0, Nelems-1 do begin 
solarabuNarrup[*,i]=solarabuN[i] + abs(randomn(seed, Nmc/2.)*solarabuNerrup[i])
solarabuNarrdn[*,i]=solarabuN[i] - abs(randomn(seed, Nmc/2.)*solarabuNerrdn[i])
solarabuNarr[*,i]=shuffle([solarabuNarrup[*,i], solarabuNarrdn[*,i]])
endfor

;;;>>>if use symetric, conservative error bars 
if symetriclable eq 'YES' then for i=0, Nelems-1 do solarabuNarr[*,i] = solarabuN[i] + randomn(see, Nmc)*solarabuNerr[i]

;stop

;;Earth
earthabuN=(bulkearthppm/atomwt)/(bulkearthppm[nref]/atomwt[nref]) * 100. ;1e6
;bulkearthppmerr=bulkearthppmerrup ;; for C, O, Mg, Si, Fe, Ca, Al, Na, and S, error bars of Wang+2018 are symetric
earthabuNerrup=(bulkearthppmerrup/bulkearthppm)*earthabuN ;
earthabuNerrdn=(bulkearthppmerrdn/bulkearthppm)*earthabuN 

earthabuNarrup=fltarr(Nmc/2., Nelems)
earthabuNarrdn=fltarr(Nmc/2., Nelems)
earthabuNarr=fltarr(Nmc, Nelems)
for i=0, Nelems-1 do begin
earthabuNarrup[*,i]=earthabuN[i] + abs(randomn(seed, Nmc/2.)*earthabuNerrup[i])
earthabuNarrdn[*,i]=earthabuN[i] - abs(randomn(seed, Nmc/2.)*earthabuNerrdn[i])
earthabuNarr[*,i]=shuffle([earthabuNarrup[*,i], earthabuNarrdn[*,i]])
endfor


earthabuNerr = earthabuNerrup
errdiff = earthabuNerrup - earthabuNerrdn
earthabuNerr(where(errdiff lt 0)) = earthabuNerrdn(where(errdiff lt 0))
if symetriclable eq 'YES' then for i=0, Nelems-1 do earthabuNarr[*,i] = earthabuN[i] + randomn(see, Nmc)*earthabuNerr[i]

;;;>>write out protosolar and earth abundances (to Al=1e6) for other use
;fname=outpath+'earthabuN.txt'
;writecol, fname, elemid, earthabuN, earthabuNerrup, earthabuNerrdn, fmt='(A, F,F,F)'  

if starlabel eq 'earth' then begin
   planetabuN=earthabuN 
   planetabuNerr=earthabuNerr
;   planetabuNerrdn=earthabuNerrdn
   goto, sim
endif


starabudex=stardex+refsolardex  ;;;differential dex;; irrelevant to solar abundances

if keyword_set(withreferr) then begin
   starabudexerrup=sqrt(stardexerr^2. + refsolardexerr^2.)
   starabudexerrdn=sqrt(stardexerrdn^2. + refsolardexerrdn^2.)
   if starlabel eq 'solar' or starlabel eq 'protosun' then begin
      starabudexerrup=refsolardexerr
      starabudexerrdn=refsolardexerrdn
   endif
endif else begin
   starabudexerrup=stardexerr
   starabudexerrdn=stardexerrdn
;   if starlabel eq 'protosun' then starabudexerrdn=protosundexerrdn
endelse


;nref=nAl ;Fe ;nAl
if finite(stardex[nAl]) eq 0 AND finite(stardex[nCa]) eq 1 then nref=nCa
if finite(stardex[nAl]) eq 0 AND finite(stardex[nCa]) eq 0 AND finite(stardex[nFe]) eq 1 then nref=nFe
starabuN=10^(starabudex - starabudex[nref]) * 100 ;1e6 
starabuNerrup=(10^starabudexerrup - 1)*starabuN
starabuNerrdn=(1-10^(-starabudexerrdn))*starabuN
;;---as an option for stellar key ratios only---
starNerrup=(10^stardexerr - 1)*starabuN
starNerrdn=(1-10^(-stardexerrdn))*starabuN
;;-------------------------------------------

;
starabuNarrup=fltarr(Nmc/2., Nelems)
starabuNarrdn=fltarr(Nmc/2., Nelems)
starabuNarr=fltarr(Nmc, Nelems)
for i=0, Nelems-1 do begin
starabuNarrup[*,i]=starabuN[i] + abs(randomn(seed, Nmc/2.)*starabuNerrup[i])
starabuNarrdn[*,i]=starabuN[i] - abs(randomn(seed, Nmc/2.)*starabuNerrdn[i])
starabuNarr[*,i]=shuffle([starabuNarrup[*,i], starabuNarrdn[*,i]])
endfor

starabuNerr=starabuNerrup
errdiff = starabuNerrup - starabuNerrdn
starabuNerr(where(errdiff lt 0)) = starabuNerrdn(where(errdiff lt 0))
if symetriclable eq 'YES' then for i=0, Nelems-1 do starabuNarr[*,i] = starabuN[i] + randomn(see, Nmc)*starabuNerr[i]


;;;>>>INDEPENDENT CALCULATIONS BY ERROR PROPAGATION 
;;STAR's planet
;;>>3sigma
if keyword_set(F3sig) then begin
ymodelNsdup=ymodelN3lmtup-ymodelN
ymodelNsddn=ymodelN-ymodelN3lmtdn
ymodeldexsdup=ymodeldex3lmtup - ymodellog
ymodeldexsddn=ymodellog - ymodeldex3lmtdn

if keyword_set(upperlmt) then begin
ymodelN=ymodelN3lmtup
ymodelNsdup=fltarr(Nelems)
ymodelNsddn=fltarr(Nelems)
ymodellog=ymodeldex3lmtup
ymodeldexsdup=fltarr(Nelems)
ymodeldexsddn=fltarr(Nelems)
endif
if keyword_set(lowerlmt) then begin
ymodelN=ymodelN3lmtdn
ymodelNsdup=fltarr(Nelems)
ymodelNsddn=fltarr(Nelems)
ymodellog=ymodeldex3lmtdn
ymodeldexsdup=fltarr(Nelems)
ymodeldexsddn=fltarr(Nelems)
endif
if keyword_set(upperlmt) and keyword_set(lowerlmt) then print, 'WARNING: Both upper and lower limits are set; lower limit is applied...'
endif else begin
;;>>F5sig?
if keyword_set(F5sig) then begin
ymodelNsdup=ymodelN5lmtup-ymodelN
ymodelNsddn=ymodelN-ymodelN5lmtdn
ymodeldexsdup=ymodeldex5lmtup - ymodellog
ymodeldexsddn=ymodellog - ymodeldex5lmtdn

if keyword_set(upperlmt) then begin
ymodelN=ymodelN5lmtup
ymodelNsdup=fltarr(Nelems)
ymodelNsddn=fltarr(Nelems)
ymodellog=ymodeldex5lmtup
ymodeldexsdup=fltarr(Nelems)
ymodeldexsddn=fltarr(Nelems)
endif
if keyword_set(lowerlmt) then begin
ymodelN=ymodelN5lmtdn
ymodelNsdup=fltarr(Nelems)
ymodelNsddn=fltarr(Nelems)
ymodellog=ymodeldex5lmtdn
ymodeldexsdup=fltarr(Nelems)
ymodeldexsddn=fltarr(Nelems)
endif
if keyword_set(upperlmt) and keyword_set(lowerlmt) then print, 'WARNING: Both upper and lower limits are set; lower limit is applied...'
endif else begin
if keyword_set(upperlmt) or keyword_set(lowerlmt) then begin
print, 'Please indicate the sigma level (F3sig or F5sig) if upper or lower limit is set'
stop
endif
endelse

;;>>if neither F3sig or F5sig set 
if keyword_set(upperlmt) or keyword_set(lowerlmt) then begin
print, 'Please indicate the sigma level (F3sig or F5sig) if upper or lower limit is set'
stop
endif
;;otherwise, run normally with 1sigma uncertainties
print, 'Apply normally the best-fit devol pattern (with 1sigma uncertainties)...'
endelse

if keyword_set(F3sig) and keyword_set(F5sig) then print, 'WARNING: Both 3- and 5-sigma are set; 5-sigma is applied...'

;;>>devol pattern applied
planetabuN=starabuN*ymodelN  ;;;starabuN and ymodelN are totally independent sources
planetabuNerrup=sqrt(starabuNerrup^2.*ymodelN^2. + ymodelNsdup^2.*starabuN^2.) ;;up-up
planetabuNerrdn=sqrt(starabuNerrdn^2.*ymodelN^2. + ymodelNsddn^2.*starabuN^2.) ;;dn-dn

;stop

if keyword_set(residual_add) then begin
planetabuN0=planetabuN
planetabuN=planetabuN0 + residual
endif
;stop

planetabuN_toFe=planetabuN/planetabuN[nFe]*1e2
planetabuN_toFeerrup=planetabuN_toFe*(planetabuNerrup/planetabuN)
planetabuN_toFeerrdn=planetabuN_toFe*(planetabuNerrdn/planetabuN)

planetabuNarrup=make_array(Nmc/2., Nelems, value=nan)
planetabuNarrdn=make_array(Nmc/2., Nelems, value=nan)
planetabuNarr=fltarr(Nmc, Nelems)
for i=0, Nelems-1 do begin
planetabuNarrup[*,i]=planetabuN[i] + abs(randomn(seed, Nmc/2.)*planetabuNerrup[i])
planetabuNarrdn[*,i]=planetabuN[i] - abs(randomn(seed, Nmc/2.)*planetabuNerrdn[i])
planetabuNarr[*,i]=shuffle([planetabuNarrup[*,i], planetabuNarrdn[*,i]])
endfor

planetabuNerr = planetabuNerrup
errdiff = planetabuNerrup - planetabuNerrdn
planetabuNerr(where(errdiff lt 0)) = planetabuNerrdn(where(errdiff lt 0))
if symetriclable eq 'YES' then for i=0, Nelems-1 do planetabuNarr[*,i] = planetabuN[i] + randomn(see, Nmc)*planetabuNerr[i]


if finite(starabudex[nAl]) eq 0 and finite(starabudex[nCa]) eq 1 then nref=nCa
if finite(starabudex[nAl]) eq 0 and finite(starabudex[nCa]) eq 0 AND finite(starabudex[nFe]) eq 1 then nref=nFe
starabudex2Al=starabudex - starabudex[nref]

ymodeldex=ymodellog
planetabudex2Al=ymodeldex + starabudex2Al
;planetabudex2Alerr0=sqrt(ymodeldexsdmax^2. + starabudex2Alerr^2.)
;stop

ymodeldexarrup=fltarr(Nmc/2., Nelems)
ymodeldexarrdn=fltarr(Nmc/2., Nelems)
ymodeldexarr = fltarr(Nmc, Nelems)
for i=0, Nelems-1 do begin 
ymodeldexarrup[*,i]=ymodellog[i] + abs(randomn(seed, Nmc/2.)*ymodeldexsdup[i])
ymodeldexarrdn[*,i]=ymodellog[i] - abs(randomn(seed, Nmc/2.)*ymodeldexsddn[i]) 
ymodeldexarr[*,i]=shuffle([ymodeldexarrup[*,i], ymodeldexarrdn[*,i]])
endfor

ymodeldexsd = ymodeldexsdup
errdiff = ymodeldexsdup - ymodeldexsddn
ymodeldexsd(where(errdiff lt 0)) = ymodeldexsddn(where(errdiff lt 0))
if symetriclable eq 'YES' then for i=0, Nelems-1 do ymodeldexarr[*,i] = ymodellog[i] + randomn(seed, Nmc)*ymodeldexsd[i]


stardexarrup=fltarr(Nmc/2., Nelems)
stardexarrdn=fltarr(Nmc/2., Nelems)
stardexarr=fltarr(Nmc, Nelems)
for i=0, Nelems-1 do begin
   stardexarrup[*,i]=starabudex2Al[i] + abs(randomn(seed, Nmc/2.)*starabudexerrup[i])
   stardexarrdn[*,i]=starabudex2Al[i] - abs(randomn(seed, Nmc/2.)*starabudexerrdn[i])   
   stardexarr[*,i]=shuffle([stardexarrup[*,i], stardexarrdn[*,i]])
endfor

starabudexerr = starabudexerrup
errdiff = starabudexerrup - starabudexerrdn
starabudexerr(where(errdiff lt 0)) = starabudexerrdn(where(errdiff lt 0))
if symetriclable eq 'YES' then for i=0, Nelems-1 do stardexarr[*,i] = starabudex2Al[i] + randomn(see, Nmc)*starabudexerr[i]


planetdexarr=stardexarr + ymodeldexarr ;ymodellogarrc


;if keyword_set(xplot) then window, /free, title='Histograms of Planet Abundance in Dex'
planetabudex2Alerr=make_array(Nelems, value=nan)
for i=0, Nelems-1 do begin
if  where(DPelemshow eq elemid[i]) NE -1 AND finite(starabudex2Al[i]) eq 1 then begin
histogauss, planetdexarr[*,i], A, /noplot
planetabudex2Alerr[i]=A[2] ;; A[0]-height, A[1]-mean, A[2]-weight of histo (standard deviation)
;stop
endif
endfor


;;;>>cal the differential abundances for planets, for further calculating key elemental ratios in dex
planetdex=planetabudex2Al - (refsolardex - refsolardex[nref])  ;;2Al2sun, namely [X/Ca]
planetdexerr=planetabudex2Alerr

;print, 'output planet abundance in dex to file'
;writecol, outpath+starlabel+'exopldex.txt', atomid, elemid, planetdex, planetdexerr, fmt='(i2, A, F,F)'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print, 'Cal key elemental ratios in both linear and logarithm...'
;;;;;;
planetC2O=planetabuN[nC] / planetabuN[nO]
print, starlabel+' planetary C/O: ', planetC2O  ;;According to Unterborn+2014, C would be oxidized after Fe, and thus it is reansonable to consider C as the last. ; But this must be the case when C/O < 0.8; otherwise, diamond or carbide (e.g. FeC, CaC) mantle would form.  
Mantleicity0=(planetabuN[nMg] + 2.*(planetabuN[nSi])) / (planetabuN[nO])
print, starlabel+' planetary mean: (Mg + 2Si) / O: ', Mantleicity0
Mantleicity1=total([planetabuN[nMg], 2.*(planetabuN[nSi]),  planetabuN[nCa], 3/2.*planetabuN[nAl]], /nan) / (planetabuN[nO])
print, starlabel+' planetary mean: (Mg + 2Si + Ca + 3/2Al) / O: ', Mantleicity1


;stop
;;;>>>>cal the abundance ratios of the bodies themselves
;;solar
;set_plot, 'x'
ratiolable='star'
bodyabuN=starabuN
bodyabuNarr=starabuNarr
bodyabuNerrup=starabuNerrup ;starabuNerrup
bodyabuNerrdn=starabuNerrdn ;starabuNerrdn

calratio:
;Nmc=1e5
C2O=bodyabuN[nC]/bodyabuN[nO]
C2ONarr=bodyabuNarr[*,nC]/bodyabuNarr[*,nO]
C2Oerr_up=sqrt(bodyabuNerrup[nC]^2./bodyabuN[nO]^2. + bodyabuNerrdn[nO]^2.*(bodyabuN[nC]/bodyabuN[nO]^2.)^2.)
C2Oerr_dn=sqrt(bodyabuNerrdn[nC]^2./bodyabuN[nO]^2. + bodyabuNerrup[nO]^2.*(bodyabuN[nC]/bodyabuN[nO]^2.)^2.)

if keyword_set(xplot) then begin
window, 1, title=starlabel+ratiolable+' C/O'
histogauss, C2ONarr, A
rectangular, A[1]-A[2], 0, A[1]+A[2], A[0]
rectangular, C2O-C2Oerr_dn, 0, C2O+C2Oerr_up, A[0], linestyle=2
endif else histogauss, C2ONarr, A, /noplot
C2Oerr=A[2]
;stop

Mg2Si=bodyabuN[nMg]/bodyabuN[nSi]
Mg2SiNarr=bodyabuNarr[*,nMg]/bodyabuNarr[*,nSi]
Mg2Sierr_up=sqrt(bodyabuNerrup[nMg]^2./bodyabuN[nSi]^2. + bodyabuNerrdn[nSi]^2.*(bodyabuN[nMg]/bodyabuN[nSi]^2.)^2.)
Mg2Sierr_dn=sqrt(bodyabuNerrdn[nMg]^2./bodyabuN[nSi]^2. + bodyabuNerrup[nSi]^2.*(bodyabuN[nMg]/bodyabuN[nSi]^2.)^2.)

if keyword_set(xplot) then begin
window, 2, title=starlabel+ratiolable+' Mg/Si'
histogauss, Mg2SiNarr, A 
rectangular, A[1]-A[2], 0, A[1]+A[2], A[0]
rectangular, Mg2Si-Mg2Sierr_dn, 0, Mg2Si+Mg2Sierr_up, A[0], linestyle=2
endif else histogauss, Mg2SiNarr, A ,/noplot
Mg2Sierr=A[2]

maxFeO= bodyabuN[nO] - (bodyabuN[nMg] + bodyabuN[nSi]*2.)
mantlecity=maxFeO/bodyabuN[nFe]
maxFeONarr=bodyabuNarr[*,nO] - (bodyabuNarr[*,nMg] + bodyabuNarr[*,nSi]*2.)
mantlecityNarr=maxFeONarr/bodyabuNarr[*,nFe]
maxFeOerr_up=sqrt(bodyabuNerrup[nO]^2. + bodyabuNerrdn[nMg]^2. + (bodyabuNerrdn[nSi]*2.)^2.)
maxFeOerr_dn=sqrt(bodyabuNerrdn[nO]^2. + bodyabuNerrup[nMg]^2. + (bodyabuNerrup[nSi]*2.)^2.)
mantlecityerr_up=sqrt(maxFeOerr_up^2./bodyabuN[nFe]^2. + bodyabuNerrdn[nFe]^2.*(maxFeO/bodyabuN[nFe]^2.)^2.)
mantlecityerr_dn=sqrt(maxFeOerr_dn^2./bodyabuN[nFe]^2. + bodyabuNerrup[nFe]^2.*(maxFeO/bodyabuN[nFe]^2.)^2.)

if keyword_set(xplot) then begin 
window, 3, title=starlabel+ratiolable+' (O-Mg-2Si)/Fe'
histogauss, mantlecityNarr, A 
rectangular, A[1]-A[2], 0, A[1]+A[2], A[0]
rectangular, mantlecity-mantlecityerr_dn, 0, mantlecity+mantlecityerr_up, A[0], linestyle=2
endif else histogauss, mantlecityNarr, A ,/noplot
mantlecityerr=A[2]
;STOP



Fe2Mg=bodyabuN[nFe]/bodyabuN[nMg]
Fe2MgNarr=bodyabuNarr[*,nFe]/bodyabuNarr[*,nMg]
Fe2Mgerr_up=sqrt(bodyabuNerrup[nFe]^2./bodyabuN[nMg]^2. + bodyabuNerrdn[nMg]^2.*(bodyabuN[nFe]/bodyabuN[nMg]^2.)^2.)
Fe2Mgerr_dn=sqrt(bodyabuNerrdn[nFe]^2./bodyabuN[nMg]^2. + bodyabuNerrup[nMg]^2.*(bodyabuN[nFe]/bodyabuN[nMg]^2.)^2.)

;if finite(bodyabuN[nNi]) eq 1 then begin
if keyword_set(xplot) then begin
window, 3, title=starlabel+ratiolable+' Fe/Mg'
histogauss, Fe2MgNarr, A
rectangular, A[1]-A[2], 0, A[1]+A[2], A[0]
rectangular, Fe2Mg-Fe2Mgerr_dn, 0, Fe2Mg+Fe2Mgerr_up, A[0], linestyle=2
endif else histogauss, Fe2MgNarr, A,/noplot
Fe2Mgerr=A[2]
;endif

;;>>(Mg + Si)/Fe
MgSi2Fe=(bodyabuN[nSi]+bodyabuN[nMg])/bodyabuN[nFe]
MgSi2FeNarr=(bodyabuNarr[*,nSi]+bodyabuNarr[*,nMg])/bodyabuNarr[*,nFe]
MgSi2Feerr_up=sqrt((bodyabuNerrup[nSi]^2. + bodyabuNerrup[nMg]^2.)/bodyabuN[nFe]^2. + bodyabuNerrdn[nFe]^2.*((bodyabuN[nMg]+bodyabuN[nSi])/bodyabuN[nFe]^2.)^2.)
MgSi2Feerr_dn=sqrt((bodyabuNerrdn[nSi]^2. + bodyabuNerrdn[nMg]^2.)/bodyabuN[nFe]^2. + bodyabuNerrup[nFe]^2.*((bodyabuN[nMg]+bodyabuN[nSi])/bodyabuN[nFe]^2.)^2.)

;if finite(bodyabuN[nNi]) eq 1 then begin
if keyword_set(xplot) then begin
window, 4, title=starlabel+ratiolable+' (Mg+Si)/Fe'
histogauss, MgSi2FeNarr, A
rectangular, A[1]-A[2], 0, A[1]+A[2], A[0]
rectangular, MgSi2Fe-MgSi2Feerr_dn, 0, MgSi2Fe+MgSi2Feerr_up, A[0], linestyle=2
endif else histogauss, MgSi2FeNarr, A,/noplot
MgSi2Feerr=A[2]


;;>>Eu/Mg
Eu2Mg=bodyabuN[nEu]*1e8/bodyabuN[nMg] ;;;1e8 is artificially chosen to enlarge Eu value to make the ratio to be float
Eu2MgNarr=bodyabuNarr[*,nEu]*1e8/bodyabuNarr[*,nMg]
Eu2Mgerr_up=sqrt((bodyabuNerrup[nEu]*1e8)^2./bodyabuN[nMg]^2. + bodyabuNerrdn[nMg]^2.*(bodyabuN[nEu]*1e8/bodyabuN[nMg]^2.)^2.)
Eu2Mgerr_dn=sqrt((bodyabuNerrdn[nEu]*1e8)^2./bodyabuN[nMg]^2. + bodyabuNerrup[nMg]^2.*(bodyabuN[nEu]*1e8/bodyabuN[nMg]^2.)^2.)

if keyword_set(xplot) then begin
window, 5, title=starlabel+ratiolable+' Eu/Mg'
histogauss, Eu2MgNarr, A
rectangular, A[1]-A[2], 0, A[1]+A[2], A[0]
rectangular, Eu2Mg-Eu2Mgerr_dn, 0, Eu2Mg+Eu2Mgerr_up, A[0], linestyle=2
endif else histogauss, Eu2MgNarr, A,/noplot
Eu2Mgerr=A[2]


;;>>Eu/(Mg+Si)
Eu2MgSi=bodyabuN[nEu]*1e8/(bodyabuN[nMg]+bodyabuN[nSi]) ;;;1e8 is artificially chosen to enlarge Eu value to make the ratio to be float
Eu2MgSiNarr=bodyabuNarr[*,nEu]*1e8/(bodyabuNarr[*,nMg] + bodyabuNarr[*,nSi])
Eu2MgSierr_up=sqrt((bodyabuNerrup[nEu]*1e8)^2./(bodyabuN[nMg]+bodyabuN[nSi])^2. + (bodyabuNerrdn[nMg]^2. + bodyabuNerrdn[nSi]^2.)*(bodyabuN[nEu]*1e8/(bodyabuN[nMg] + bodyabuN[nSi])^2.)^2.)
Eu2MgSierr_dn=sqrt((bodyabuNerrdn[nEu]*1e8)^2./(bodyabuN[nMg]+bodyabuN[nSi])^2. + (bodyabuNerrup[nMg]^2. + bodyabuNerrup[nSi]^2.)*(bodyabuN[nEu]*1e8/(bodyabuN[nMg] + bodyabuN[nSi])^2.)^2.)

if keyword_set(xplot) then begin
window, 6, title=starlabel+ratiolable+' Eu/(Mg+Si)'
histogauss, Eu2MgSiNarr, A
rectangular, A[1]-A[2], 0, A[1]+A[2], A[0]
rectangular, Eu2MgSi-Eu2MgSierr_dn, 0, Eu2MgSi+Eu2MgSierr_up, A[0], linestyle=2
endif else histogauss, Eu2MgSiNarr, A,/noplot
Eu2MgSierr=A[2]


;;>>fo2_proxy: O/(Mg + 2Si)
fo2_p=bodyabuN[nO] / (bodyabuN[nMg] + bodyabuN[nSi]*2.)
fo2_parr=bodyabuNarr[*, nO] / (bodyabuNarr[*, nMg] + bodyabuNarr[*, nSi]*2.)
fo2_perrup = sqrt(bodyabuNerrup[nO]^2./(bodyabuN[nMg] + bodyabuN[nSi]*2.)^2. + (bodyabuNerrdn[nMg]^2. + (bodyabuNerrdn[nSi]*2)^2.)*(bodyabuN[nO] / (bodyabuN[nMg] + bodyabuN[nSi]*2.)^2.)^2.)
fo2_perrdn = sqrt(bodyabuNerrdn[nO]^2./(bodyabuN[nMg] + bodyabuN[nSi]*2.)^2. + (bodyabuNerrup[nMg]^2. + (bodyabuNerrup[nSi]*2)^2.)*(bodyabuN[nO] / (bodyabuN[nMg] + bodyabuN[nSi]*2.)^2.)^2.)
if keyword_set(xplot) then begin
window, 7, title=starlabel+ratiolable+ ' O/(Mg+2Si)'
histogauss, fo2_parr, A
rectangular, A[1]-A[2], 0, A[1]+A[2], A[0]
rectangular, fo2_p - fo2_perrdn, 0, fo2_p + fo2_perrup, A[0], linestyle=2
endif else histogauss, fo2_parr, A,/noplot
fo2_perr=A[2]


case ratiolable of
'star': begin
starN_keyR=[C2O, Mg2Si, Fe2Mg, Eu2Mg, mantlecity, fo2_p, MgSi2Fe, Eu2MgSi]
starN_keyRerr=[C2Oerr, Mg2Sierr, Fe2Mgerr, Eu2Mgerr, mantlecityerr, fo2_perr, MgSi2Feerr, Eu2MgSierr]
starN_keyRerr_up=[C2Oerr_up, Mg2Sierr_up, Fe2Mgerr_up, Eu2Mgerr_up, mantlecityerr_up, fo2_perrup, MgSi2Feerr_up, Eu2MgSierr_up]
starN_keyRerr_dn=[C2Oerr_dn, Mg2Sierr_dn, Fe2Mgerr_dn, Eu2Mgerr_dn, mantlecityerr_dn, fo2_perrdn,  MgSi2Feerr_dn, Eu2MgSierr_dn]
;stop
ratiolable='planet'
bodyabuN=planetabuN
bodyabuNarr=planetabuNarr
bodyabuNerrup=planetabuNerrup
bodyabuNerrdn=planetabuNerrdn
;stop
goto, calratio
end
'planet': begin
planetN_keyR=[C2O, Mg2Si, Fe2Mg, Eu2Mg, mantlecity, fo2_p, MgSi2Fe, Eu2MgSi]
planetN_keyRerr=[C2Oerr, Mg2Sierr,  Fe2Mgerr, Eu2Mgerr,mantlecityerr,  fo2_perr, MgSi2Feerr, Eu2MgSierr]
planetN_keyRerr_up=[C2Oerr_up, Mg2Sierr_up, Fe2Mgerr_up, Eu2Mgerr_up, mantlecityerr_up, fo2_perrup, MgSi2Feerr_up, Eu2MgSierr_up]
planetN_keyRerr_dn=[C2Oerr_dn, Mg2Sierr_dn, Fe2Mgerr_dn, Eu2Mgerr_dn, mantlecityerr_dn, fo2_perrdn, MgSi2Feerr_dn, Eu2MgSierr_dn]
;stop
ratiolable='solar'
bodyabuN=solarabuN
bodyabuNarr=solarabuNarr
bodyabuNerrup=solarabuNerrup
bodyabuNerrdn=solarabuNerrdn
;if sunearthlable EQ 'YES' then 
goto, calratio
end
'solar': begin
solarN_keyR=[C2O, Mg2Si, Fe2Mg, Eu2Mg, mantlecity, fo2_p, MgSi2Fe, Eu2MgSi]
solarN_keyRerr=[C2Oerr, Mg2Sierr, Fe2Mgerr, Eu2Mgerr, mantlecityerr, fo2_perr, MgSi2Feerr, Eu2MgSierr]
solarN_keyRerr_up=[C2Oerr_up, Mg2Sierr_up, Fe2Mgerr_up, Eu2Mgerr_up, mantlecityerr_up, fo2_perrup, MgSi2Feerr_up, Eu2MgSierr_up]
solarN_keyRerr_dn=[C2Oerr_dn, Mg2Sierr_dn, Fe2Mgerr_dn, Eu2Mgerr_dn, mantlecityerr_dn, fo2_perrdn, MgSi2Feerr_dn, Eu2MgSierr_dn]
;stop
ratiolable='earth'
bodyabuN=earthabuN
bodyabuNarr=earthabuNarr
bodyabuNerrup=earthabuNerrup
bodyabuNerrdn=earthabuNerrdn
;if sunearthlable eq 'YES' then 
goto, calratio
end
'earth': begin
earthN_keyR=[C2O, Mg2Si, Fe2Mg, Eu2Mg, mantlecity, fo2_p, MgSi2Fe, Eu2MgSi]
earthN_keyRerr=[C2Oerr, Mg2Sierr,  Fe2Mgerr, Eu2Mgerr, mantlecityerr, fo2_perr, MgSi2Feerr, Eu2MgSierr]
earthN_keyRerr_up=[C2Oerr_up, Mg2Sierr_up, Fe2Mgerr_up, Eu2Mgerr_up, mantlecityerr_up, fo2_perrup, MgSi2Feerr_up, Eu2MgSierr_up]
earthN_keyRerr_dn=[C2Oerr_dn, Mg2Sierr_dn, Fe2Mgerr_dn, Eu2Mgerr_dn, mantlecityerr_dn, fo2_perrdn, MgSi2Feerr_dn, Eu2MgSierr_dn]
end
;endratioN:
endcase
keyRlable=['C/O', 'Mg/Si', 'Fe/Mg', 'Eu/Mg', '(O-Mg-2Si)/Fe', 'O/(Mg+2Si)', '(Mg+Si)/Fe', 'Eu/(Mg+Si)']
print, starlabel+'starkeyR: '+keyRlable, starN_keyR
print, starlabel+'starkeyRerr: '+keyRlable, starN_keyRerr
print, starlabel+'planetkeyR: '+keyRlable, planetN_keyR
print, starlabel+'planetkeyRerr: '+keyRlable, planetN_keyRerr

;stop
writecol, outpath+starlabel+'_starabuN.txt', elemid, starabuN, starabuNerrup, starabuNerrdn, fmt='(A, F, F, F)'
writecol, outpath+starlabel+'_PlanetabuN.txt', elemid, PlanetabuN, PlanetabuNerrup, PlanetabuNerrdn, fmt='(A, F, F, F)' 
writecol, outpath+starlabel+'_PlanetabuN_toFe.txt', elemid, PlanetabuN_toFe, PlanetabuN_toFeerrup, PlanetabuN_toFeerrdn, fmt='(A, F, F, F)'
writecol, outpath+starlabel+'_NkeyR_'+solarreflabel+'.txt',  keyRlable, starN_keyR, starN_keyRerr, starN_keyRerr_up, starN_keyRerr_dn, planetN_keyR, planetN_keyRerr, planetN_keyRerr_up, planetN_keyRerr_dn,fmt='(A, F,F,F,F,F,F,F,F)'
;if sunearthlable eq 'YES' then 
writecol, outpath+'RefSunEarth_NkeyR_'+solarreflabel+'.txt',  keyRlable, solarN_keyR, solarN_keyRerr, solarN_keyRerr_up, solarN_keyRerr_dn, earthN_keyR, earthN_keyRerr, earthN_keyRerr_up, earthN_keyRerr_dn, fmt='(A, F, F,F,F,F,F,F,F)'

;stop

;;;>>>>>>>>RATIOS IN DEX<<<<<<<<<<<<<<<<<<<<<<
;;Key elemental ratios, e.g. [Mg/Si]_p = [Mg/Al]_p - [Si/Al]_p;; 
keyRlable='Mg2Si'
nelem1=nMg
nelem2=nSi
stardexc=stardex
;if outputref2protosun eq 'YES' then stardexc= starabudex - protosundex
;if starlabel eq 'solar' then planetdexerr=ymodeldexsdmax

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

keyRlable='Fe2Mg'
nelem1=nFe
nelem2=nMg
starFe2Mg=stardexc[nelem1] - stardexc[nelem2] ;- stardex[nO]
starFe2Mgerr=sqrt(stardexerr[nelem1]^2. + stardexerr[nelem2]^2.) ; + stardexerr[nO]^2.)
;

planetFe2Mg=(planetdex[nelem1] - planetdex[nelem2]) ;- planetdex[nFe] ; - planetdex[nO]
planetFe2Mgerr=sqrt(planetdexerr[nelem1]^2. + planetdexerr[nelem2]^2.)


;if starlabel eq 'solar' then planetdexarr = ymodeldexarr

IF histoforkeyRdexlable eq 'YES' THEN BEGIN
starMg2SiNarr=stardexarr[*,nMg] - stardexarr[*,nSi]
histogauss, starMg2SiNarr, A, /noplot
starMg2Sierr=A[2]
;
planetMg2SiNarr=planetdexarr[*,nMg] - planetdexarr[*,nSi]
if keyword_set(xplot) then begin
window, /free, title=starlabel+ 'b '+keyRlable+' in DEX'
histogauss, planetMg2SiNarr, A
endif else histogauss, planetMg2SiNarr, A, /noplot
planetMg2Sierr=A[2]

starC2ONarr=stardexarr[*,nC] - stardexarr[*,nO]
histogauss, starC2ONarr, A, /noplot
starC2Oerr=A[2]
;
planetC2ONarr=planetdexarr[*,nC] - planetdexarr[*,nO]
if keyword_set(xplot) then begin
window, /free, title=starlabel+ 'b '+keyRlable+' in DEX'
histogauss, planetC2ONarr, A
endif else histogauss, planetC2ONarr, A ,/noplot
planetC2Oerr=A[2]

starFe2SiNarr=stardexarr[*,nFe] - stardexarr[*,nSi]
histogauss, starFe2SiNarr, A, /noplot
starFe2Sierr=A[2]
;
planetFe2SiNarr=planetdexarr[*,nFe] - planetdexarr[*,nSi]
if keyword_set(xplot) then begin
window, /free, title=starlabel+ 'b '+keyRlable+' in DEX'
histogauss, planetFe2SiNarr, A
endif else histogauss, planetFe2SiNarr, A,/noplot
planetFe2Sierr=A[2]

starFe2MgNarr=stardexarr[*,nFe] - stardexarr[*,nMg]
histogauss, starFe2MgNarr, A, /noplot
starFe2Mgerr=A[2]
;
planetFe2MgNarr=planetdexarr[*,nFe] - planetdexarr[*,nMg]
if keyword_set(xplot) then begin
window, /free, title=starlabel+ 'b '+keyRlable+' in DEX'
histogauss, planetFe2MgNarr, A
endif else histogauss, planetFe2MgNarr, A,/noplot
planetFe2Mgerr=A[2]

ENDIF

KeyRdexlable=['[Fe/Si]', '[Fe/Mg]', '[Mg/Si]', '[C/O]']
planetkeyR=[planetFe2Si, planetFe2Mg, planetMg2Si, planetC2O];,  planetM2O]
planetkeyRerr=[planetFe2Sierr, planetFe2Mgerr, planetMg2Sierr, planetC2Oerr];,  planetM2Oerr] >>>these errors have not included the errors on the reference solar abundances, so its error is smaller than the actual values as determined in the linear scale above. 
starkeyR=[starFe2Si, starFe2Mg, starMg2Si, starC2O];, starM2O]
starkeyRerr=[starFe2Sierr, starFe2Mgerr, starMg2Sierr, starC2Oerr];, starM2Oerr]

;stop

;writecol, outpath+starlabel+'keyRdex.txt',  keyRdexlable, planetkeyR, planetkeyRerr, starkeyR, starkeyRerr, fmt='(A, F,F,F,F)'

;;;;Check the sufficiency of O per every group of compounds
;;;;;;;;;;;;;;;;

;;;;;>>>>set the situation with no depletion
val_id=where(finite(starabuN) eq 1)
IF Keyword_set(nodevostar) THEN BEGIN
planetabuN=starabuN
planetabuNarr=starabuNarr
print, 'No Devol used for stochiometric cal of compounds!'
ENDIF

;stop

sim:
Print, 'Simul Using the best Mean Values of Abundances...'
;stop
;Results_m=chemsysmodel_upgrade(planetabuN, /outinfull)
chemsyslabel='new'
if keyword_set(quick) then Results_m=chemsysmodel(planetabuN, /outinfull, /quick) else  Results_m=chemsysmodel(planetabuN, /outinfull)

print, 'Overall valid results from uniform MC or not:', Results_m.Nval
print, 'Normalised valid runs for each planet abundance input: '
print, 'Mantle: ', Results_m.mantle_Nval
print, 'Core: ', Results_m.core_Nval
print, 'CMF: ', Results_m.cmf_Nval

meanprocess:
sizemantle=size(Results_m.mantlecompsmassfra)
Nmantlecomps=sizemantle[1] ;;14, as designed, changable
sizecore=size(Results_m.corecompsmassfra)
Ncorecomps=sizecore[1];; 3, as designed, changable
perc1=[0.1587, 0.50, 0.8413] ;; the chosen 1-sigma percentiles, with 0th 

if Results_m.Nval le 1 then begin 
mantlecompsmassfraBEST=Results_m.mantlecompsmassfra
corecompsmassfraBEST=Results_m.corecompsmassfra
fcoremassBEST=Results_m.fcoremass

endif else begin 
   mantlecompsmassfraBEST=make_array(Nmantlecomps, value=nan)
   mantlecompsmassfraBESTperc=make_array(Nmantlecomps, value=nan)
   for i=0, Nmantlecomps-1 do  begin
   val_id=where(finite(Results_m.mantlecompsmassfra[i,*]) eq 1, count, /null) 
   mantlecompsmassfraBEST[i]=mean(Results_m.mantlecompsmassfra[i,*], /nan)
if count ge 1 then mantlecompsmassfraBESTperc[i, *]=cgpercentiles(Results_m.mantlecompsmassfra[i,val_id], percentiles=[0.50])
endfor

;
corecompsmassfraBEST=make_array(Ncorecomps, value=nan)
corecompsmassfraBESTperc=make_array(Ncorecomps, value=nan)
for i=0, Ncorecomps-1 do begin
   val_id=where(finite(Results_m.corecompsmassfra[i,*]) eq 1, count, /null)
   corecompsmassfraBEST[i]=mean(Results_m.corecompsmassfra[i,*], /nan)
if count ge 1 then   corecompsmassfraBESTperc[i, *]=cgpercentiles(Results_m.corecompsmassfra[i,val_id], percentiles=[0.50])   
endfor
;
fcoremassBEST=mean(Results_m.fcoremass, /nan)
fcoremassBESTperc=nan; make_array(3, value=nan)
val_id=where(finite(Results_m.fcoremass) eq 1, count, /null)
if count ge 1 then fcoremassBESTperc=cgpercentiles(Results_m.fcoremass[val_id], percentiles=[0.50])


percusedforbest='YES'
if percusedforbest eq 'YES' then begin
   mantlecompsmassfraBEST=mantlecompsmassfraBESTperc
   corecompsmassfraBEST=corecompsmassfraBESTperc
   fcoremassBEST=fcoremassBESTperc
endif
endelse




print, 'Mean of BEST abundances of '+starlabel
print, 'C_mantle: Mantle stoichimetric composition (BEST mean):', Results_m.mantlecompsname, mantlecompsmassfraBEST
print, 'C_core: Core stoichimetric composition (BEST mean):', Results_m.corecompsname, corecompsmassfraBEST
print, 'F_core: Core mass fraction (BEST mean):', fcoremassBEST

print, 'Simul Using MC to Produce Error Bars...'
;;;>>>Generate randomn arrays for planet abundances (with uncertainties)
Ncut=Nsim*0.01 ;;less than 1/100 of the total simulation is not trusted. 
if Ncut lt 1. then Ncut=1.


planetabuNarrMCup=make_array(Nsim/2., Nelems, value=nan)
planetabuNarrMCdn=make_array(Nsim/2., Nelems, value=nan)
planetabuNarrMC=make_array(Nsim, Nelems, value=nan)
for i=0, Nelems-1 do begin
planetabuNarrMCup[*,i]=planetabuN[i] + abs(randomn(seed, Nsim/2.)*planetabuNerrup[i])
planetabuNarrMCdn[*,i]=planetabuN[i] - abs(randomn(seed, Nsim/2.)*planetabuNerrdn[i])
planetabuNarrMC[*,i]=shuffle([planetabuNarrMCup[*,i], planetabuNarrMCdn[*,i]])
endfor



if symetriclable eq 'YES' then for i=0, Nelems-1 do planetabuNarrMC[*,i] = planetabuN[i] + randomn(seed, Nsim)*planetabuNerr[i]


IF oxidesMClable eq 'U' THEN BEGIN
for i=0, Nelems-1 do begin
if finite(planetabuN[i]) eq 1 then planetabuNarrMC[*, i] = (planetabuN[i]-planetabuNerrdn[i]) + randomu(seed, Nmc)*(planetabuNerrup[i]+planetabuNerrdn[i])  ;;uniform distribution 
endfor
ENDIF



mantlecompsmassfraARR=[ ]
corecompsmassfraARR=[ ]
fcoremassARR=[ ]
mantlestoicmassARR=[]
corestoicmassARR=[]
mantle_sim_Nval=fltarr(Nmantlecomps)
core_sim_Nval=fltarr(Ncorecomps)
cmf_sim_Nval=0.


IF keyword_set(quick) THEN BEGIN
for i=0, Nsim-1 do begin
print, starlabel+' MC simul', i+1, ' of ', Nsim
Results=chemsysmodel(planetabuNarrMC[i, *], /outinfull, /quick)


if finite(Results.Nval) eq 0 then Ninval=Ninval+1. 
;;;;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mantlecompsmassfraARR=[[mantlecompsmassfraARR], [Results.mantlecompsmassfra]] ;;might be multi dimensional, then should have [ ]
corecompsmassfraARR=[[corecompsmassfraARR], [Results.corecompsmassfra]]
;Ncoreval=n_elements(Results.corecompsmassfra)
fcoremassARR=[fcoremassARR, Results.fcoremass] ;; Results.fcoremass is 1 scalar or 1-dimension results, then concatenate by column
mantlestoicmassARR=[[mantlestoicmassARR], [Results.mantlestoicmass]]
corestoicmassARR=[[corestoicmassARR], [Results.corestoicmass]]

;;>>count the valid outputs for each compound/quantity of each
;;component
mantle_sim_Nval=mantle_sim_Nval + Results.mantle_Nval
core_sim_Nval=core_sim_Nval + Results.core_Nval
cmf_sim_Nval=cmf_sim_Nval + Results.cmf_Nval
print, 'mantle_sim_nval: ', mantle_sim_nval
;stop
endfor
ENDIF ELSE BEGIN
   for i=0, Nsim-1 do begin
print, starlabel+' MC simul', i+1, ' of ', Nsim
;;;;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Results=chemsysmodel(planetabuNarrMC[i, *], /outinfull)


if finite(Results.Nval) eq 0 then Ninval=Ninval+1. 
;;;;>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
mantlecompsmassfraARR=[[mantlecompsmassfraARR], [Results.mantlecompsmassfra]] ;;might be multi dimensional, then should have [ ]
corecompsmassfraARR=[[corecompsmassfraARR], [Results.corecompsmassfra]]
;Ncoreval=n_elements(Results.corecompsmassfra)
fcoremassARR=[fcoremassARR, Results.fcoremass] ;; Results.fcoremass is 1 scalar or 1-dimension results, then concatenate by column
mantlestoicmassARR=[[mantlestoicmassARR], [Results.mantlestoicmass]]
corestoicmassARR=[[corestoicmassARR], [Results.corestoicmass]]


;;>>count the valid outputs for each compound/quantity of each
;;component
mantle_sim_Nval=mantle_sim_Nval + Results.mantle_Nval
core_sim_Nval=core_sim_Nval + Results.core_Nval
cmf_sim_Nval=cmf_sim_Nval + Results.cmf_Nval
print, 'mantle_sim_nval: ', mantle_sim_nval
;stop
endfor
ENDELSE




print, 'valid sims among ', Nsim, 'of MC sims'
print, 'mantle_sim_Nval: ', mantle_sim_Nval
print, 'core_sim_Nval: ', core_sim_Nval
print, 'cmf_sim_Nval: ', cmf_sim_Nval
print, 'Valid run ratio for each compound/quantity of each component under ', Nsim, ' planet abu MC inputs'
R_Nsim_m=mantle_sim_Nval/Nsim
R_Nsim_c=core_sim_Nval/Nsim
R_Nsim_cmf=cmf_sim_Nval/Nsim
print, 'Mantle: ', Results.mantlecompsname, R_Nsim_m
print, 'Core: ', Results.corecompsname, R_Nsim_c
print, 'CMF: ', R_Nsim_cmf

;;CONCONATE THE BEST DATASET
mantlecompsmassfraARR=[[mantlecompsmassfraARR], [Results_m.mantlecompsmassfra]] ;;might be multi dimensional, then should have [ ]
corecompsmassfraARR=[[corecompsmassfraARR], [Results_m.corecompsmassfra]]
fcoremassARR=[fcoremassARR, Results_m.fcoremass]
mantlestoicmassARR=[[mantlestoicmassARR], [Results_m.mantlestoicmass]]
corestoicmassARR=[[corestoicmassARR], [Results_m.corestoicmass]]


;stop

fullcompsname=Results.mantlecompsname
fullcorecompsname=Results.corecompsname
Nfullmantle=N_ELEMENTS(fullcompsname)
Nfullcore=N_ELEMENTS(fullcorecompsname)
nFeO=where(fullcompsname eq 'FeO')
nMetal=where(fullcompsname eq 'Metals')
;stop
;;Cal the mean, standard deviations, and [0.1587, 0.5, 0.8413]
;;percentiles  of the results
mantlecompsmassframean=make_array(Nfullmantle, value=nan)
mantlecompsmassfraERR=make_array(Nfullmantle, value=nan)
mantlecompsmassfraMAX=make_array(Nfullmantle, value=nan)
mantlecompsmassfraMIN=make_array(Nfullmantle, value=nan)
mantlecompsmassfra_perc=make_array(Nfullmantle, 3, value=nan)
mantlecompsmassfra_boot=make_array(Nfullmantle, 3, value=nan)
mantlecompsmassfrac_perc=make_array(Nfullmantle, 3, value=nan) ;;constrained results
mantlestoicmassARRc = mantlestoicmassARR
mantlecompsmassfraARRc = mantlecompsmassfraARR

mantle_val_count=fltarr(Nfullmantle)
for i=0, Nfullmantle-1 do begin
val_id=where(finite(mantlecompsmassfraARR[i,*]) eq 1, count, /null)
mantlestoicmassARRc[i,*] = mantlestoicmassARR[i,*]*R_Nsim_m[i] 
mantlecompsmassfraARRc[i,*] = mantlecompsmassfraARR[i,*]*R_Nsim_m[i]
;mantle_val_count[i]=count  
   if count ge Ncut then begin
      mantlecompsmassframean[i]=mean(mantlecompsmassfraARR[i,*], /nan)
mantlecompsmassfraERR[i]=stddev(mantlecompsmassfraARR[i,*], /nan)
mantlecompsmassfraMAX[i]=max(mantlecompsmassfraARR[i,*], /nan)
mantlecompsmassfraMIN[i]=min(mantlecompsmassfraARR[i,*], /nan)
mantlecompsmassfra_perc[i,*]=cgpercentiles(mantlecompsmassfraARR[i, val_id], percentiles=[0.1587, 0.50, 0.8413]) ;; cgpercentiles would return different results if nan is included in an array. 
;mantlecompsmassfra_boot[i,*] = BOOTSTRAP_MEDIAN(mantlecompsmassfraARR[i, val_id]) ;;default number of boostrap resampling: 1000; default confidence limit: 0.68
mantlecompsmassfrac_perc[i,*]=cgpercentiles(mantlecompsmassfraARRc[i, val_id], percentiles=[0.1587, 0.50, 0.8413]) 
endif
endfor
corecompsmassframean=make_array(Nfullcore, value=nan)
corecompsmassfraERR=make_array(Nfullcore, value=nan)
corecompsmassfraMAX=make_array(Nfullcore, value=nan)
corecompsmassfraMIN=make_array(Nfullcore, value=nan)
corecompsmassfra_perc=make_array(Nfullcore, 3, value=nan)
corecompsmassfra_boot=make_array(Nfullcore, 3, value=nan)
corecompsmassfrac_perc=make_array(Nfullcore, 3, value=nan)
corestoicmassARRc = corestoicmassARR
corecompsmassfraARRc = corecompsmassfraARR

core_val_count=fltarr(Nfullcore)
for i=0, Nfullcore-1 do begin
   val_id=where(finite(corecompsmassfraARR[i,*]) eq 1, count, /null)
corestoicmassARRc[i,*] = corestoicmassARR[i,*]*R_Nsim_c[i]
corecompsmassfraARRc[i,*] = corecompsmassfraARR[i,*]*R_Nsim_c[i]
;core_val_count[i]=count
if count ge Ncut then begin
corecompsmassframean[i]=mean(corecompsmassfraARR[i,*], /nan)
   corecompsmassfraERR[i]=stddev(corecompsmassfraARR[i,*], /nan)
corecompsmassfraMAX[i]=max(corecompsmassfraARR[i,*], /nan)
corecompsmassfraMIN[i]=min(corecompsmassfraARR[i,*], /nan)
corecompsmassfra_perc[i,*]=cgpercentiles(corecompsmassfraARR[i, val_id], percentiles=[0.1587, 0.50, 0.8413])
;corecompsmassfra_boot[i,*] = BOOTSTRAP_MEDIAN(corecompsmassfraARR[i, val_id])
corecompsmassfrac_perc[i,*]=cgpercentiles(corecompsmassfraARRc[i, val_id], percentiles=[0.1587, 0.50, 0.8413])
endif
endfor
fcoremassmean=mean(fcoremassARR, /nan)
fcoremassERR=stddev(fcoremassARR, /nan)
fcoremassMAX=max(fcoremassARR, /nan)
fcoremassMIN=min(fcoremassARR, /nan)
cmf_perc=make_array(3, value=nan)
cmf_boot=make_array(3, value=nan)
val_id=where(finite(fcoremassARR) eq 1, cmf_count, /null)

totalmantlestoicmassc = total(mantlestoicmassARRc, 1, /nan)
totalcorestoicmassc = total(corestoicmassARRc, 1, /nan)
fcoremassARRc = totalcorestoicmassc / (totalmantlestoicmassc + totalcorestoicmassc)
cmfc_perc = make_array(3, value=nan)
if cmf_count ge Ncut then begin
cmf_perc=cgpercentiles(fcoremassARR(val_id), percentiles=[0.1587, 0.5, 0.8413])
cmfc_perc = cgpercentiles(fcoremassARRc(val_id), percentiles=[0.1587, 0.5, 0.8413])
endif

;;;>>>>>>>>output the histograms for check and verification

;;>>>Info of Earth as reference
mantlefra_control=make_array(Nfullmantle, value=nan)
mantleE=[0.0036, 0.0355, 0.378, 0.0445, 0.45, 0.0805, 0.0025, NAN, NAN, NAN, NAN, NAN]
mantlefra_control=mantleE


corefraE=[0.8671, 0.0527, 0.0193]
if Nfullcore eq 4 then corefraE=[corefraE , 0.0609]
corefra_control=corefraE/total(corefraE)

cmf_control=0.325 ; from Wang et al. 2018, which is based on geophysical constraints. 
cmferr_control=0.003



;;;>>>>write out the ARRAY FILES
;WRITE_CSV, outpath+'mantleARR_'+starlabel+'.csv', mantlecompsmassfraARR, header=fullcompsname
;WRITE_CSV, outpath+'coreARR_'+starlabel+'.csv', corecompsmassfraARR, header=fullcorecompsname
;WRITE_CSV, outpath+'cmfARR_'+starlabel+'.csv', fcoremassARR
;paranames1=['Nfullmantle', 'Nfullcore', 'Ncut']
;paravalues1=[Nfullmantle,  Nfullcore,  Ncut]
;WRITE_CSV, outpath+'parameters_1_'+starlabel+'.csv', Nfullmantle,  Nfullcore,  Ncut, header=paranames1
;paranames2=['fullcompsname', 'mantle_sim_Nval', 'R_Nsim_m']
;WRITE_CSV, outpath+'parameters_2_'+starlabel+'.csv', fullcompsname, mantle_sim_Nval, R_Nsim_m, header=paranames2
;paranames3=['fullcorecompsname','core_sim_Nval', 'R_Nsim_c']
;WRITE_CSV, outpath+'parameters_3_'+starlabel+'.csv', fullcorecompsname, core_sim_Nval, R_Nsim_c, header=paranames3
;paranames4=['cmf_sim_Nval', 'R_Nsim_cmf']
;WRITE_CSV, outpath+'parameters_4_'+starlabel+'.csv', cmf_sim_Nval, R_Nsim_cmf, header=paranames4
;percnames1=['fullcompsname','mantlecompsmassfra_perc16', 'mantlecompsmassfra_perc50', 'mantlecompsmassfra_perc84']
;WRITE_CSV, outpath+'perc_1_'+starlabel+'.csv', fullcompsname, mantlecompsmassfra_perc[*,0],mantlecompsmassfra_perc[*,1], mantlecompsmassfra_perc[*,2], header=percnames1
;percnames2=['fullcorecompsname','corecompsmassfra_perc16', 'corecompsmassfra_perc50', 'corecompsmassfra_perc84']
;WRITE_CSV, outpath+'perc_2_'+starlabel+'.csv',  fullcorecompsname, corecompsmassfra_perc[*,0], corecompsmassfra_perc[*,1], corecompsmassfra_perc[*,2], header=percnames2
;percnames3=['cmf_perc16', 'cmf_perc50', 'cmf_perc84']
;WRITE_CSV, outpath+'perc_3_'+starlabel+'.csv', cmf_perc[0], cmf_perc[1], cmf_perc[2], header=percnames3

;stop
;;;--OUTPUT RESULTS---
if starlabel eq 'earth' then result={mantlecompsname:fullcompsname, mantlecompsmassfraBEST:mantlecompsmassfraBEST,  mantle_perc:mantlecompsmassfra_perc, corecompsname:fullcorecompsname, corecompsmassfraBEST:corecompsmassfraBEST, core_perc:corecompsmassfra_perc, fcoremassBEST:fcoremassBEST, cmf_perc:cmf_perc, Nsim:Nsim, R_Nsim_m:R_Nsim_m, R_Nsim_c:R_Nsim_c, R_Nsim_cmf:R_Nsim_cmf} else $

Result={mantlecompsname:fullcompsname, mantlecompsmassfraBEST:mantlecompsmassfraBEST,  mantle_perc:mantlecompsmassfra_perc, mantlec_perc:mantlecompsmassfrac_perc,  corecompsname:fullcorecompsname, corecompsmassfraBEST:corecompsmassfraBEST, core_perc:corecompsmassfra_perc, corec_perc:corecompsmassfrac_perc,  fcoremassBEST:fcoremassBEST, cmf_perc:cmf_perc, cmfc_perc:cmfc_perc, planetdex:planetdex, planetdexerr:planetdexerr, planetkeyR:planetkeyR, planetkeyRerr:planetkeyRerr, starkeyR:starkeyR, starkeyRerr:starkeyRerr, KeyRdexlable:KeyRdexlable, Nsim:Nsim, R_Nsim_m:R_Nsim_m, R_Nsim_c:R_Nsim_c, R_Nsim_cmf:R_Nsim_cmf}

TOC
return, result
END
