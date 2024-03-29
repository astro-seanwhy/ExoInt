Function chemsysmodel, planetabuN, outinfull=outinfull, quick=quick
TIC
;;;;;;DO NOT MODIFY THIS CODE (EXCEPT ADJUSTING KEYWORDS IF NECESSARY), WITHOUT THE PERMISSION OF THE CREATER
;;;;;;A Geo-&Cosmo-model [CODE-3] for estimating rocky exoplanetary
;;composition and structures: mantle and core composition, core mass fraction, core radius.
;;Created and maintained by Haiyang S. Wang
;;Citation: Wang, H. S., Liu, F., Ireland, T., Brasser, R., Yong, D.,
;;and Lineweaver, C. H. 2019. Enhanced constraints on the interior
;;composition and structure of terrestrial exoplanets. MNRAS 482:2222-2233. doi.org/10.1093/mnras/sty2749

;;History of updates
;;--------7 Feb 2020 ----------
;;1. the notion 'mantleSO2abu' has been changed to 'mantleSO3abu';
;;just a change to the notion, rather than to the variable, since the
;;variable had already been calculauted under the asummption of SO3
;;2. correct the line for calculating the remaining oxygen after
;;oxdizing C to CO2, from
;;;    planetOleft_res2=planetOleft_res - planetabuN[nC]*3.
;;;to
;;;    planetOleft_res2=planetOleft_res - planetabuN[nC]*2.
;;;And the following line
;;;    compsabu=planetabuN[nC] + planetOleft_res2/3.
;;;to 
;;;    compsabu=planetabuN[nC] + planetOleft_res2/2.
;;;As a result, subtle changes to the amounts of CO2 and eventually remaining
;;oxygen atoms would be expected. But the changes will have no effect
;;to the results of planet host cases in Wang et al. 2019 MNRAS, since
;;no case has enough oxygen to oxidze C to CO2 and thus has oxygen remained. 
;;-------------------------------
;;--------10 Apr 2020--------------
;;3. The uniform distribution for dividing Fe, Ni, and S is employed
;;to calculate the mean of the division and output as the result of
;;one iteration with a set of planetary abundances drawn from MC
;;simulation
;;4. SiO2 was put backward after Al2O3 in the oxidation sequence to
;;allow that Si may be in the core when oxygen is not enough. 
;;5. Na2O should be firstly oxidsed than CaO, as according to Johnson
;;& Warr 1999. This is a correction to what was mistaken in Wang et
;;al. 2019, but there is no affect to the modelling results, since
;;both Na and Ca would also be fully oxidised in all cases. 
;;-------2 May 2020 --------------
;;6. Making Si in the core as a choice, rather than being
;;determined. ;; according to Takashi & McDonough 2020, Si is not
;;considered in the Mars' core, with respect to the P-T-fO2
;;condition of the Martian interior (also refer to Wade & Wood, 2005;
;;Corgne et al. 2008).
;;----------------------------------------
;;-----22 March 2021 ----------------
;;7. updated the treatment of molar mass fractions, but this update is more for 
;;style change and basically has no effect on final results.  

Nelems=83
nan=!values.f_nan
secondaryOxidizeslable='NO' ;; if do secondary oxidizes: TiO2, Cr2O3, MnO, K2O, and P2O5. 
 ;; YES OR NO to consider metals and graphite/diamond in the mantle and un-oxidized Si into the core. 
;SiO2firstlable='YES'
SiIncorelable='YES' ;;otherwsie in the mantle as Si, which may form other non-oidies, e.g. SiC by combining with C
Oincorelable='NO' ;; otherwise in the mantle to be ready for other minor oxides
MantleMetallable='YES'
unimeanlable='YES' ;;use the mean from the uniform draws of Fe, Ni, S or Si as the output
meanmolarmass='YES'
molarfralabel='NO'

Nmc=1e5
;if keyword_set(quick) then Nmc=1e3

readcol, 'data/atomwttc_new.txt', F='I,A,F,F', atomid, elemid, atomwt, elemtc
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

;;process the planetabuN, in which the negative abundance
;;(e.g. generated by MC occasionally) may be set as N/A
planetabuN(where(planetabuN LT 0.))=nan  ;;the negative abundance randomly generated by Gaussian distribution is unrealistic, and thus set to be nan then. 
Nvalc=1 ;; initially set number of groups of results, to be updated upon the random distribution of Fe, Ni, and S in the mantle and the core
;;;;set the limits of NiO and SO4 firstly
coreFe2NiUL=18 + 4. ;; Assume 17.6 +/-4 is drawn from the distribution of Fe/Ni in the more than 5000 stars of Hypatia. Assume 17 +/- 1 as the range of ratio of Fe to Ni in the core, based on McDonough 2017 and personal communication with McDonough. 
coreFe2NiLL=18 - 4. 

;;preset the C (graphite/diamond) in mantle, metals (could be All, but
;;Si and S, of elements considered here) in mantle
CmantleN=nan ;; Graphite/diamond in mantle
MetalN=nan   ;; metals will only be in the form of Fe and/or Ni when oxygen is not enough and also out of the Fe/Ni constraints
MetalName=nan
Metalwt=nan
SinativeN=nan
OnativeN=nan
ExtraO=nan ;;extra oxygen will be put into the core instead. 

;;;>>>>SET UP THE FRAME
compsname=['Na2O', 'CaO', 'MgO', 'Al2O3', 'SiO2', 'FeO', 'NiO', 'SO3', 'CO2', 'C'] ;;['SiO2', 'CaO', 'Na2O', 'MgO', 'Al2O3', 'FeO', 'NiO', 'SO4', 'CO2', 'C']
compswt=[atomwt[nNa]*2.+atomwt[nO], atomwt[nCa]+atomwt[nO], atomwt[nMg]+atomwt[nO], atomwt[nAl]*2+atomwt[nO]*3, $
         atomwt[nSi]+atomwt[nO]*2., $
         atomwt[nFe]+atomwt[nO], atomwt[nNi]+atomwt[nO], atomwt[nS]+atomwt[nO]*3., $
         atomwt[nC]+atomwt[nO]*2., atomwt[nC]]
if MantleMetallable eq 'YES' then begin
   compsname=[compsname, 'Metals'];; the metals are mainly Fe and Ni, but plus any element which is not oxidised except Fe, Ni and S. 
   compswt=[compswt, nan] ;; the metal atomwt TBD later. 
   mantlemetals=[]
   metalabu=[]
   metalwt=[]
endif


corecompsname=['Fe', 'Ni', 'S'];, 'Si']
corecompswt=[atomwt[nFe], atomwt[nNi], atomwt[nS]];, atomwt[nSi]]
if Siincorelable eq 'YES' then begin
   corecompsname=[corecompsname, 'Si']
   corecompswt=[corecompswt, atomwt[nSi]]
endif 

if Oincorelable ne 'YES' then begin
   compsname=[compsname, 'ExtraO'];; extra oxygen can never be put in the core, since that is the case that both Fe and Ni have been fully oxidised and trapped in the mantle. 
   compswt=[compswt, atomwt[nO]]
endif

Ncomps=n_elements(compsname)

Ncorecomps=n_elements(corecompsname)

print, 'mantle compounds: ', compsname
print, 'core constituents: ', corecompsname
compsabu=replicate(nan, Ncomps)
corecompsabu=replicate(nan, Ncorecomps)
mantlecompsmassfra=replicate(nan, Ncomps) ;;9; 14 is the total mantle compounds number as designed, changable 
corecompsmassfra=replicate(nan, Ncorecomps) ;; 3 is the total core compounds assumed, changable
fcoremass=nan
mantlemolarmeanmass=nan
coremolarmeanmass=nan
Rmolarmeanmass_m2c=nan
compsmolarmassfra=replicate(nan, Ncomps)
corecompsmolarmassfra=replicate(nan, Ncorecomps)
mantlemolarmass = replicate(nan, Ncomps)
coremolarmass = replicate(nan, Ncorecomps)

Nval=0.
mantle_Nval=fltarr(Ncomps)
core_Nval=fltarr(Ncorecomps)
cmf_Nval=0.
;;;;;<<<<<<<<
;stop

;;>>>Filter 1--Elements that must have
avaielemid=elemid(where(finite(planetabuN) eq 1, /null))
print, 'elements with valid abundances by MC draws: ', avaielemid
IF where(avaielemid eq 'O') eq -1 OR where(avaielemid eq 'Mg') eq -1 OR where(avaielemid eq 'Si') eq -1 OR where(avaielemid eq 'Fe') eq -1 THEN BEGIN
print, 'WARNING: You must have measurements for O, Mg, Si, and Fe, AT LEAST, for estimating a reasonable planetary composition !'
print, 'Skip the sim after waiting for 2s!'
Wait, 2
;stop
goto, outputresults
ENDIF

;;>>>>negative values
negid=where(planetabuN lt 0, negcount, /null)
if negcount ne 0 then begin
print, 'the following elements have NEGATIVE abundances -- unvalid draw. RETURN.'
print, elemid(negid)
print, 'Skip the sim after waiting for 2s!'
Wait, 2
stop
goto, outputresults
endif


;;>>>> Filter2 --THIS IS ONLY APPLICABLE TO STAR ABUNDANCES WITH NO
;;DEVO APPLICATION YET. 
planetN_C2O=planetabuN[nC] / planetabuN[nO]
IF planetN_C2O ge 0.8 THEN BEGIN
print, 'Carbide Planet !!'
print, 'It may form CO, SiC, Mg2C, Fe3C, and CaC2... beyond the scope of a silicate planet!'
print, 'Calculation is DISMISSED; N/A Results for all are returned'
goto, outputresults
ENDIF


;;;;;>>>pre-set the arrays for two cases (i.e. when S is available but
;;;;;assumed all in mantle or all in core) 

;;OTHERWISE: SILICATE PLANET
;;do the order of ease of oxidation of metals, according to the
;;lithophibility of elements and the ease of oxidation as in Johnson
;;and Warr 1999. The order should be the following:
;;Level 1: Na, Ca, Mg, and Al -- sorted in oder of ease of oxidation and they
;;                      are also all lithophiles and supposed to be
;;fully oxidised in mantle  (it may be unreasonable to have no
;;sufficient oxygen to oxidize these four main pure lithophiles -- and
;;                                                                 thus
;;such draw of planet abundances with MC will be excluded)

;;Level 2: Si  -- the ease of oxidation of Si is not listed in Johnson
;;                & Warr 1999, but since Si may appear in the core as
;;in McDonough 2003 (Si-bearing core model) and the meta-analysis in
;;Wang et al. 2018, it is reasonable to put it follow the previous
;;lithophiles. If oxygen atoms are not sufficient to fully oxidze Si
;;atoms, then the remaining Si metal will be put into the core.
;;Level 3: Fe, Ni, and S ;; Fe has been listed less to be oxidized
;;than Al, so it is reasonable to put it here (but after Si -- a more
;;                                                             lithophile
;;than Fe).
;; Level 3: C --> CO2 or graphite in the mantle
;; Level 4: O --> ExtraO in the mantle or otherwise O in the core

;;Calculate the second cretiria -- Oxygen to (Na/2 + Ca + Mg + Al*3/2)

;;>> LEVEL-1 OXIDES: Na2O, CaO, MgO, Al2O3
OtoLithR=planetabuN[nO]/total([planetabuN[nNa]/2. + planetabuN[nCa] + planetabuN[nMg] + planetabuN[nAl]*3/2], /nan)
IF OtoLithR LT 1 THEN BEGIN
   print, 'O/(Na/2+Ca+Mg+Al*3/2) lt 1 -- An invalid MC draw of planet abundances', OtoLithR
   print, 'Re-draw abudances in process...'
   goto, outputresults
ENDIF ELSE BEGIN
;   compsname=['Na2O', 'CaO', 'MgO', 'Al2O3']
;   compswt=[atomwt[nNa]*2.+atomwt[nO], atomwt[nCa]+atomwt[nO], atomwt[nMg]+atomwt[nO], atomwt[nAl]*2+atomwt[nO]*3]
   compsabu[0:3]=[planetabuN[nNa]/2., planetabuN[nCa], planetabuN[nMg], planetabuN[nAl]/2.]

   planetOleft=planetabuN[nO] - total([planetabuN[nNa]/2. + planetabuN[nCa] + planetabuN[nMg] + planetabuN[nAl]*3./2.], /nan) ;;3/2=1; 3/2.=1.5!
ENDELSE

;;>> LEVEL-2 OXIDES: SiO2

if planetOleft/2. lt planetabuN[nSi] then begin
   print, 'O has been used up after forming SiO2'
   compsabu[where(compsname eq 'SiO2')]=planetOleft/2.

   case Siincorelable of
      'YES': corecompsabu[where(corecompsname eq 'Si')]=planetabuN[nSi]-planetOleft/2.
      'NO':  metalSiabu=planetabuN[nSi]-planetOleft/2.
   endcase
      goto, modelcore      
endif else begin
   compsabu[where(compsname eq 'SiO2')]=planetabuN[nSi]
   planetOleft=planetOleft - planetabuN[nSi]*2.
endelse


;; LEVEL-3 OXIDES: FeO, (NiO, and SO3)

;;;>>>ASSUME ALL O LEFT WOULD BE CONSUMED UP BY Fe, Ni, and/or S. 
Nmc=Nmc ;; sampling numbers of MC
NiOmax=planetabuN[nNi] ;; it is fair considering >90% Ni into the core, in the Earth, McDonough 2017
NiOmin=0. 
SO3max=planetabuN[nS] ;; it is fair considering >95% S into the core, in the Earth, McDonough 2017
SO3min=0.

;;;surfur is processed together with Fe and Ni (the two should be always bounded together, as their cosmochemical ratio of 17 in various meteorites) 
CASE finite(planetabuN[nS]) OF 
0: BEGIN  ;; NO S 
print, 'Planet has NO S...'

IF planetOleft ge total([planetabuN[nFe], planetabuN[nNi]], /nan) THEN BEGIN ;;ASSUME METALS TO BE OXIDIZED IF THERE IS ENGOUGH OXYGEN. 
print, 'A Coreless planet is formed!'
mantleFeOabu=planetabuN[nFe]
mantleNiOabu=planetabuN[nNi]
compsabu[where(compsname eq 'FeO')]=mantleFeOabu
compsabu[where(compsname eq 'NiO')]=mantleNiOabu

planetOleft=planetOleft-total([planetabuN[nFe], planetabuN[nNi]], /nan) 

;;;>>>>>C as either CO2 or otherwise graphite/kerogen/diamond
if planetOleft/2. le planetabuN[nC] then begin
   compsabu[where(compsname eq 'CO2')]=planetOleft/2.
   compsabu[where(compsname eq 'C')]=planetabuN[nC] - planetOleft/2.
   print, 'O atoms have been exhusted after forming CO2!'
endif else begin
   compsabu[where(compsname eq 'CO2')]=planetabuN[nC]
   planetOleft_res=planetOleft - planetabuN[nC]*2.

   compsabu(where(compsname eq 'ExtraO'))=planetOleft_res
   print, 'There are extra O atoms of: ', planetOleft_res
   print, 'And they may be used for forming other minor oxides that are not considered yet - e.g. MnO, P2O5, K2O, etc.'

;   if Oincorelable eq 'YES' then begin  ;;; FOR A CORELESS PLANET, IT
;   IS IRATIONAL TO HAVE O IN THE CORE
;      print, 'Opt to put the extra O into the core...'
;      corecompsabu[(where(corecompsname eq 'O')]=planetOleft_res
;   endif
endelse
goto, chemsysdone
ENDIF ELSE BEGIN
mantleNiOabuarr=NiOmin + (NiOmax-NiOmin)*randomu(seed, Nmc)  ;;;legitimate to assume uniform distribution, as the values between the max and min could be equally possible. 
mantleFeOabuarr=planetOleft - mantleNiOabuarr

;;CORE
coreNiabuarr=planetabuN[nNi] - mantleNioabuarr
coreFeabuarr=planetabuN[nFe] - mantleFeOabuarr


;;>>>Validity Control 1. Important!
inval_id1=where(mantleFeOabuarr lt 0., /null)
;if where(inval_id1 eq -1) eq -1 then begin ;;namely, exclude the case inval_id1= -1, namely, not fit to the conditions
mantleFeOabuarr(inval_id1)=nan ;; FeO cannot be negative 
mantleNiOabuarr(inval_id1)=nan ;; the amount of NiO making FeO to be negative is also invalid. 
coreFeabuarr(inval_id1)=nan
coreNiabuarr(inval_id1)=nan

;;>>>Validity Control 2. Important!
coreFe2Niarr=coreFeabuarr/coreNiabuarr
inval_id2=where(coreFe2Niarr lt coreFe2NiLL OR coreFe2Niarr gt coreFe2NiUL, /null)
;if where(inval_id2 eq -1) eq -1 then begin
coreFe2Niarr(inval_id2)=nan
coreNiabuarr(inval_id2)=nan
coreFeabuarr(inval_id2)=nan
mantleNiOabuarr(inval_id2)=nan
mantleFeOabuarr(inval_id2)=nan
;endif
;;

val_id=where(finite(mantleFeOabuarr) eq 1, Nval, /null) ;;FeO is the most stable reference as Fe must be present in planetabuN
print, 'Numbers of valid guess for Fe and Ni in the mantle and the core:', Nval, ' of', Nmc
if Nval eq 0 then begin
   print, 'NO valid UNIFORM-MC solutions for the given planet abundances. '
   print, 'output NAN results'
   goto, outputresults
endif

mantleFeOabuarrc=mantleFeOabuarr(val_id)
mantleNiOabuarrc=mantleNiOabuarr(val_id)
coreNiabuarrc=coreNiabuarr(val_id)
coreFeabuarrc=coreFeabuarr(val_id)

case Nval of
 ;  0: goto, modelcore
   1: begin
      compsabu(where(compsname eq 'FeO'))=mantleFeOabuarrc
      compsabu(where(compsname eq 'NiO'))=mantleNiOabuarrc         
      corecompsabu(where(corecompsname eq 'Fe'))=coreFeabuarrc
      corecompsabu(where(corecompsname eq 'Ni'))=coreNiabuarrc       
      goto, graphiteonly
      end
   else: begin
      compsabuarr=make_array(Ncomps, Nval, value=nan)
      corecompsabuarr=make_array(Ncorecomps, Nval, value=nan)
      compsabuarr(where(compsname eq 'FeO'), *)=mantleFeOabuarrc
      compsabuarr(where(compsname eq 'NiO'), *)=mantleNiOabuarrc
      compsabuarr(where(compsname eq 'SO3'), *)=replicate(nan, Nval)
      corecompsabuarr(where(corecompsname eq 'Fe'), *)=coreFeabuarrc
      corecompsabuarr(where(corecompsname eq 'Ni'), *)=coreNiabuarrc     
      corecompsabuarr(where(corecompsname eq 'S'), *)=replicate(nan, Nval)
      corecompsabuarr(where(corecompsname eq 'Si'), *)=replicate(nan, Nval) ;replicate(corecompsabu(where(corecompsname eq 'Si')), Nval)

      ;;>>processing for mantle compounds that are done already
      Ncompsdone=5 ; Na, Ca, Mg, Al, Si
      for i=0, 4 do compsabuarr[i,*]=replicate(compsabu[i], Nval)
      ;;>>processing for CO2 and C
      compsabuarr(where(compsname eq 'CO2'), *)=replicate(nan, Nval)
      compsabuarr(where(compsname eq 'C'), *)=replicate(planetabuN[nC], Nval)
      end
endcase

ENDELSE
END

1: BEGIN ;; IF THERE IS SULPHUR
print, 'Planet has S...'

IF planetOleft ge total([planetabuN[nFe], planetabuN[nNi], planetabuN[nS]*3.], /nan) THEN BEGIN ;;ASSUME METALS TO BE OXIDIZED IF THERE IS ENGOUGH OXYGEN. 
   print, 'A Coreless planet is formed!'
   
mantleFeOabu=planetabuN[nFe]
mantleNiOabu=planetabuN[nNi]
mantleSO3abu=planetabuN[nS]
compsabu[where(compsname eq 'FeO')]=mantleFeOabu
compsabu[where(compsname eq 'NiO')]=mantleNiOabu
compsabu[where(compsname eq 'SO3')]=mantleSO3abu

planetOleft=planetOleft-total([planetabuN[nFe], planetabuN[nNi], planetabuN[nS]*3.], /nan)


;;;>>>>>C as either CO2 or otherwise graphite/kerogen/diamond
if planetOleft/2. le planetabuN[nC] then begin
   compsabu[where(compsname eq 'CO2')]=planetOleft/2.
   compsabu[where(compsname eq 'C')]=planetabuN[nC] - planetOleft/2.
   print, 'O atoms have been exhusted after forming CO2!'
endif else begin
   compsabu[where(compsname eq 'CO2')]=planetabuN[nC]
   planetOleft_res=planetOleft - planetabuN[nC]*2.
   compsabu(where(compsname eq 'ExtraO'))=planetOleft_res
   print, 'There are extra O atoms of: ', planetOleft_res
   print, 'And they may be used for forming other minor oxides that are not considered yet - e.g. MnO, P2O5, K2O, etc.'
endelse
goto, chemsysdone
ENDIF ELSE BEGIN
mantleSO3abuarr=SO3min + (SO3max-SO3min)*randomu(seed, Nmc)
mantleNiOabuarr=NiOmin + (NiOmax-NiOmin)*randomu(seed, Nmc)
mantleFeOabuarr=planetOleft - mantleNiOabuarr - mantleSO3abuarr*3. ;; SO3

coreNiabuarr=planetabuN[nNi] - mantleNiOabuarr
coreFeabuarr=planetabuN[nFe] - mantleFeOabuarr
coreSabuarr=planetabuN[nS] - mantleSO3abuarr

;stop
;;>>>Validity Control 1. Important!
inval_id1=where(mantleFeOabuarr lt 0., /null)
;if where(inval_id1 eq -1) eq -1 then begin ;;namely, exclude the case inval_id1= -1, namely, not fit to the conditions
mantleFeOabuarr(inval_id1)=nan ;; FeO cannot be negative 
mantleNiOabuarr(inval_id1)=nan ;; the amount of NiO making FeO to be negative is also invalid. 
mantleSO3abuarr(inval_id1)=nan ;; the amount of SO4 making FeO to be negative is also invalid. 
coreFeabuarr(inval_id1)=nan ;; FeO cannot be negative 
coreNiabuarr(inval_id1)=nan ;; the amount of NiO making FeO to be negative is also invalid. 
coreSabuarr(inval_id1)=nan

;endif
;stop


;;>>>Validity Control 2. Important!
coreFe2Niarr=coreFeabuarr/coreNiabuarr
;stop
inval_id2=where(coreFe2Niarr lt coreFe2NiLL OR coreFe2Niarr gt coreFe2NiUL, /null)
;if where(inval_id2 eq -1) eq -1 then begin ;;namely, exclude the case inval_id1= -1, namely, not fit to the conditions
coreFe2Niarr(inval_id2)=nan
coreNiabuarr(inval_id2)=nan
coreFeabuarr(inval_id2)=nan
coreSabuarr(inval_id2)=nan ;; invalid for Fe and Ni, then also mean S is invalid
mantleNiOabuarr(inval_id2)=nan
mantleFeOabuarr(inval_id2)=nan
mantleSO3abuarr(inval_id2)=nan
;endif
;stop

;;>>>Validity Control 3. Important!
inval_id3=where(coreFeabuarr lt coreSabuarr, /null)
;if where(inval_id3 eq -1) eq -1 then begin
coreNiabuarr(inval_id3)=nan
coreFeabuarr(inval_id3)=nan
coreSabuarr(inval_id3)=nan
mantleNiOabuarr(inval_id3)=nan
mantleFeOabuarr(inval_id3)=nan
mantleSO3abuarr(inval_id3)=nan
;stop
;endif
;;
val_id=where(finite(mantleFeOabuarr) eq 1, Nval, /null) ;;Fe is the most stable reference as Fe must be present in planetabuN
print, 'Numbers of valid guess for Fe and Ni in the mantle and the core:', Nval, ' of', Nmc
if Nval eq 0 then begin
   print, 'NO valid UNIFORM-MC solutions for the given planet abundances. '
   print, 'output NAN results'
   goto, outputresults
endif

Fe_val_id=val_id
Ni_val_id=where(finite(mantleNiOabuarr) eq 1, Nval_Ni, /null)
S_val_id=where(finite(mantleSO3abuarr) eq 1, Nval_S, /null)
if Nval lt Nval_Ni or Nval lt Nval_S then stop


mantleFeOabuarrc=mantleFeOabuarr(val_id)
mantleNiOabuarrc=mantleNiOabuarr(val_id)
mantleSO3abuarrc=mantleSO3abuarr(val_id)
coreNiabuarrc=coreNiabuarr(val_id)
coreFeabuarrc=coreFeabuarr(val_id)
coreSabuarrc=coreSabuarr(val_id)


;stop

case Nval of
;   0: goto, modelcore
   1: begin
      compsabu(where(compsname eq 'FeO'))=mantleFeOabuarrc
      compsabu(where(compsname eq 'NiO'))=mantleNiOabuarrc      
      compsabu(where(compsname eq 'SO3'))=mantleSO3abuarrc   
      corecompsabu(where(corecompsname eq 'Fe'))=coreFeabuarrc
      corecompsabu(where(corecompsname eq 'Ni'))=coreNiabuarrc      
      corecompsabu(where(corecompsname eq 'S'))=coreSabuarrc 
      goto, graphiteonly
      end
   else: begin
      compsabuarr=make_array(Ncomps, Nval, value=nan)
      corecompsabuarr=make_array(Ncorecomps, Nval, value=nan)
      compsabuarr(where(compsname eq 'FeO'), *)=mantleFeOabuarrc
      compsabuarr(where(compsname eq 'NiO'), *)=mantleNiOabuarrc
      compsabuarr(where(compsname eq 'SO3'), *)=mantleSO3abuarrc
      corecompsabuarr(where(corecompsname eq 'Fe'), *)=coreFeabuarrc
      corecompsabuarr(where(corecompsname eq 'Ni'), *)=coreNiabuarrc     
      corecompsabuarr(where(corecompsname eq 'S'), *)=coreSabuarrc
      corecompsabuarr(where(corecompsname eq 'Si'), *)=replicate(nan, Nval) ;replicate(corecompsabu(where(corecompsname eq 'Si')), Nval)

      ;;>>processing for mantle compounds that are done already
      Ncompsdone=5 ; Na, Ca, Mg, Al, Si
      for i=0, 4 do compsabuarr[i,*]=replicate(compsabu[i], Nval)
      ;;>>processing for CO2 and C
      compsabuarr(where(compsname eq 'CO2'), *)=replicate(nan, Nval)
      compsabuarr(where(compsname eq 'C'), *)=replicate(planetabuN[nC], Nval)
      end
endcase
;STOP

ENDELSE
END
ENDCASE



;;>>>CAL MASS FRACTIONS IN CASE OF ARRAYS
compsmolarmass=make_array(Ncomps, Nval, value=nan)
mantletotalmolarmass=make_array(Nval, value=nan)
mantlecompsmassfra=make_array(Ncomps, Nval, value=nan)
;;molar
compsmolarfra=make_array(Ncomps, Nval, value=nan)
compsmolarmassfra=make_array(Ncomps, Nval, value=nan)
mantlemolarmeanmass=make_array(Nval, value=nan)
mantlemolarmassfrac=make_array(Ncomps, Nval, value=nan)

coremolarmass=make_array(Ncorecomps, Nval, value=nan)
coretotalmolarmass=make_array(Nval, value=nan)
corecompsmassfra=make_array(Ncorecomps, Nval, value=nan)
;;molar
corecompsmolarfra=make_array(Ncorecomps, Nval, value=nan) ;; relative molar fraction for each element in the core
corecompsmolarmassfra=make_array(Ncorecomps, Nval, value=nan)
coremolarmeanmass=make_array(Nval, value=nan)
corecompsmolarmassfrac=make_array(Ncorecomps, Nval, value=nan)

fcoremass=make_array(Nval, value=nan)
Rmolarmeanmass_m2c=make_array(Nval, value=nan)

;;;if meanmolarmass eq 'YES' then goto, alternativemean
;;;alternativemean:

for i=0, Nval-1 do begin
compsmolarmass[*,i] = compsabuarr[*,i]*compswt
coremolarmass[*,i] = corecompsabuarr[*,i]*corecompswt
endfor
compsmolarmass_mean = mean(compsmolarmass, dimension=2, /nan)
coremolarmass_mean = mean(coremolarmass, dimension=2, /nan)
mantletotalmolarmass = total(compsmolarmass_mean, /nan)
coretotalmolarmass = total(coremolarmass_mean, /nan)

mantlecompsmassfra = compsmolarmass_mean / mantletotalmolarmass
corecompsmassfra = coremolarmass_mean / coretotalmolarmass 
fcoremass = coretotalmolarmass / (coretotalmolarmass + mantletotalmolarmass)
Nval=1.


for i=0, Ncomps-1 do begin
   val_id_m=where(finite(mantlecompsmassfra[i]) eq 1, m_count, /null)
   if m_count ge 1 then  mantle_Nval[i]=1
endfor

for i=0, Ncorecomps-1 do begin
   val_id_c=where(finite(corecompsmassfra[i]) eq 1, c_count, /null)
   if c_count ge 1 then core_Nval[i]=1
endfor

val_id_cmf=where(finite(fcoremass) eq 1, cmf_count, /null)
if cmf_count ge 1 then cmf_Nval=1

goto, outputresults



modelcore:
corecompsabu[0:2]=[planetabuN[nFe], planetabuN[nNi], planetabuN[nS]] ;; in the case no constriants by Ni or S, then all Fe could be in the core, when no FeO in the mantle
;;Validity Control
IF finite(planetabuN[nNi]) eq 1 THEN BEGIN ;;ALL NI into the core, reaching the maximum of Fe and Ni alloy in the core, as well as the possible maximum S in the core (as no O left in mantle, preferably make the three elements into the core to the maximum)
if planetabuN[nFe] gt planetabuN[nNi]*coreFe2NiUL then begin 
corecompsabu[0]=corecompsabu[1]*coreFe2NiUL

case Siincorelable of
   'YES': begin
      compsabu[where(compsname eq 'Metals')]=[planetabuN[nFe]-corecompsabu[0]]
      compswt[where(compsname eq 'Metals')]=atomwt[nFe]
   end
   'NO': begin
      metalFeabu=planetabuN[nFe]-corecompsabu[0]
      metalabu=total([metalSiabu, metalFeabu], /nan)
      metalwt=total([metalSiabu/metalabu, metalFeabu/metalabu]*[atomwt[nSi], atomwt[nFe]], /nan)
      compsabu[where(compsname eq 'Metals')]=metalabu
      compswt[where(compsname eq 'Metals')]=metalwt
      ;stop
   end
endcase
  

endif else begin ;;making the maximum Fe under the assumed abundance of Ni
;when within the range of ratio; all Fe and Ni are just present in the core
if planetabuN[nFe] lt planetabuN[nNi]*coreFe2NiLL then begin
corecompsabu[1]=corecompsabu[0]/coreFe2NiLL

case Siincorelable of
   'YES': begin
compsabu[where(compsname eq 'Metals')]=[planetabuN[nNi]-corecompsabu[1]]
compswt[where(compsname eq 'Metals')]=atomwt[nNi]
end
   'NO': begin
      metalNiabu=planetabuN[nNi]-corecompsabu[1]
      metalabu=total([metalSiabu, metalNiabu], /nan)
      metalwt=total([metalSiabu/metalabu, metalNiabu/metalabu]*[atomwt[nSi], atomwt[nNi]], /nan)
      compsabu[where(compsname eq 'Metals')]=metalabu
      compswt[where(compsname eq 'Metals')]=metalwt
      ;stop
   end
endcase

endif
endelse
ENDIF

if corecompsabu[0] lt corecompsabu[2] then stop

graphiteonly:;;;;>>>incompleted oxiding case for mantle
;;;;>>>C should always present
compsabu[where(compsname eq 'C')]=planetabuN[nC]


chemsysdone:
Nval=1
;;Mantle
compsmolarmass=compsabu*compswt
mantletotalmolarmass=total(compsmolarmass, /nan)
mantlecompsmassfra=compsmolarmass/mantletotalmolarmass
mantlemolarmass=compsmolarmass

compsmolarfra=compsabu/total(compsabu, /nan) ;; relative molar fraction for each element in the core
mantlemolarmassfra=compsmolarfra*compswt
mantlemolarmeanmass=total(compsmolarfra*compswt, /nan)
mantlemolarmassfrac=mantlemolarmassfra/mantlemolarmeanmass ;;; confirmed to be identical to mantlecompsmassfra
;stop
;;Core
coremolarmass=corecompsabu*corecompswt ;; the average molar mass of core
coretotalmolarmass=total(coremolarmass, /nan)
corecompsmassfra=coremolarmass/coretotalmolarmass

corecompsmolarfra=corecompsabu/total(corecompsabu, /nan) ;; relative molar fraction for each element in the core
corecompsmolarmassfra=corecompsmolarfra*corecompswt
coremolarmeanmass=total(corecompsmolarfra*corecompswt, /nan)
corecompsmolarmassfrac=corecompsmolarmassfra/coremolarmeanmass ;; confirmed to be identical to corecompsmassfra

;;>>calCMF
fcoremass=coretotalmolarmass/(mantletotalmolarmass + coretotalmolarmass)
Rmolarmeanmass_m2c=mantlemolarmeanmass/coremolarmeanmass

if molarfralabel eq 'YES' then begin
mantlecompsmassfra=mantlemolarmassfrac
corecompsmassfra=corecompsmolarmassfrac
endif

;;count
for i=0, Ncomps-1 do $
   if finite(mantlecompsmassfra[i]) eq 1 then mantle_Nval[i]=1
for i=0, Ncorecomps-1 do $
   if finite(corecompsmassfra[i]) eq 1 then core_Nval[i]=1

if finite(fcoremass) eq 1 then cmf_Nval=1

;stop;---------------------------
outputresults:
;print, fcoremass
;stop

result={mantlecompsname:compsname, mantlecompsmassfra:mantlecompsmassfra, mantlestoicmass:mantlemolarmass, corecompsname:corecompsname, corecompsmassfra:corecompsmassfra, corestoicmass:coremolarmass, fcoremass:fcoremass, Nval:Nval, mantle_Nval:mantle_Nval, core_Nval:core_Nval, cmf_Nval:cmf_Nval}

TOC
return, result

;;______________________________
;TOC
;print, 'excecution done'
END
