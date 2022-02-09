PRO exoconstraints ;;THE ONLY PART OF THE CODES THAT YOU NEED TO RUN OR MODIFY
;;A Geo-&Cosmo-model [CODE-1] for estimating rocky exoplanetary
;;composition and structures: mantle and core composition, core mass fraction, core radius.
;;Created and maintained by Haiyang S. Wang
;;Citation: Wang, H. S., Liu, F., Ireland, T., Brasser, R., Yong, D.,
;;and Lineweaver, C. H. 2019. Enhanced constraints on the interior
;;composition and structure of terrestrial exoplanets. MNRAS 482:2222-2233. doi.org/10.1093/mnras/sty2749

TIC
;cgCleanup ;;This procedure will invoke X-WINDOW to clean-up and/or destroy any open graphics or widget windows on the display.
outpath='outputs/'

;;;>>>>>>>KEYS SET-UP

Nelems=83
nan=!values.f_nan

;;;>>>INPUT NECESSARY SUPPORTING DATA
readcol, 'data/atomwttc_new.txt', F='I,A,F,F', atomid, elemid, atomwt, elemtc
readcol, 'data/protosunwhy.txt', F='x,f,f,f', protosundex, protosundexerrup, protosundexerrdn
readcol, 'data/PEppmwhy.txt', F='x,d,d,d', bulkearthppm, bulkearthppmerrup, bulkearthppmerrdn
readcol, 'data/Asplund09dex.txt', F='x,x,f,f', A09solardex, A09solardexerr, skipline=3
readcol, 'data/Asplund2021dex.txt', F='x,x,f,f', A20solardex, A20solardexerr, skipline=3

;;;>>>INPUT TO-BE-DEVOLATILISED STELLAR ABUNDANCES
;;--Example star ([X/H], differntial abundances, by default)
readcol, 'data/examstar_XH.txt', F='A,f,f,f', elemid, examstardex, examstardexerrup, examstardexerrdn, skipline=1

;;;>>process invalid values
A09solardex(where(A09solardex eq -11))=nan
A09solardexerr(where(A09solardexerr eq -11))=nan
A20solardex(where(A20solardex eq -11))=nan
A20solardexerr(where(A20solardexerr eq -11))=nan

solardex=A09solardex
solardexerr=A09solardexerr
solardex20=A20solardex
solardex20err=A20solardexerr 


;;;>>>>EXAMPLE
starlabel='examstar'
stardiff=examstardex
stardifferr=examstardexerrup
stardifferrdn=examstardexerrdn ;;The parameter for the lower error bar is retained, in case that it is different from the upper error bar
Result=exomodel(stardiff, stardifferr, stardifferrdn, starlabel=starlabel, /withreferr, /residual_add, /quick) ;;/xplot, /F3sig, /F5sig, /upperlmt, /lowerlmt)
;;;/withreferr: a key (recommended to use) to determine if the errors of the reference solar abundances will be added in quadrature with your input (differential) stellar abundances.
;;;/residual_add: an optional key to determine if the residuals of the devolatilisation model will be added to the devolatilised stellar abundances to obtain the final planetary bulk composition (It is only recommended to use if it is to compare with the Earth value of Wang et al. 2018 for an acccurate comparison; it does not matter if the comparisions are among extrasolar planet cases only).
;;;/quick: a key for testing/debugging (the MC simuations are set to 100 times); omit it for full simulation (the default MC simulations are 2e4 times).
;;;extrakeys for potentially varying the devolatilisation scale ----
;;;/xplot: a plotting key (not recommended for user) that may be used together with /quick for debugging (but it's normally used for the developer only; also it may incur an X-window running error on a server if the server environment is not set properly for X11-forwarding)
;;;/F3sig or /F5sig (if set simultaneously, F5sig is activated) -- amplify the uncertainty range of the standard (1sigma) devolatilisation trend (Wang et al. 2019 Icarus) to 3 sigma or 5 sigma, respectively
;;;/upperlmt, /lowerlmt (if set simultaneously, lowerlmt is activated) -- when /F3sig or /F5sig is set, you can also choose to only apply either the upper or lower limit to devolatilise your stars only of the applied devolatilisation trend (this may be useful for doing some limit tests)

goto, outputresults ;; this line is useful for outputing results individually for multiple samples of stars. 

;;;>>>RUN FOR YOUR STARS
;youstar:
;starlabel='yourstar'
;stardiff=yourstardex
;stardifferr=yourstardexerrup
;stardifferrdn=yourstardexerrdn ;;The parameter for the lower error bar is retained, in case that it is different from the upper error bar
;Result=exomodel(stardiff, stardifferr, stardifferrdn, starlabel=starlabel, /withreferr, /residual_add, /quick)
;;;please refer to the example above for these keywords -- /withreferr, /residual_add, /quick -- and extra keys (/F3sig, /F5sig, /upperlmt, /lowerlmt)
;goto, outputresults 


;;>>output
outputresults:
print, 'Output UNCONSTRAINED (RAW) results to files for '+starlabel
writecol, outpath+starlabel+'_mantlecomp_raw.txt', Result.mantlecompsname, Result.mantlecompsmassfraBEST, Result.mantle_perc[*,0], Result.mantle_perc[*,1], Result.mantle_perc[*,2], Result.R_Nsim_m, fmt='(A, F, F, F, F, F)'
writecol, outpath+starlabel+'_corecomp_raw.txt',  Result.corecompsname, Result.corecompsmassfraBEST, Result.core_perc[*, 0], Result.core_perc[*, 1], Result.core_perc[*, 2], Result.R_Nsim_c, fmt='(A,  F, F, F, F, F)' 


print, 'Output CONSTRAINED (FINAL) results to files for '+starlabel
mantlec_per50 = Result.mantlec_perc[*,1]/total(Result.mantlec_perc[*,1], /nan) ;;; normalise to 1
mantlec_errup = (Result.mantlec_perc[*,2] - Result.mantlec_perc[*,1])/Result.mantlec_perc[*,1]*mantlec_per50
mantlec_errdn = (Result.mantlec_perc[*,1] - Result.mantlec_perc[*,0])/Result.mantlec_perc[*,1]*mantlec_per50
corec_per50 = Result.corec_perc[*,1]/total(Result.corec_perc[*,1], /nan) ;;; normalise to 1
corec_errup = (Result.corec_perc[*,2] - Result.corec_perc[*,1])/Result.corec_perc[*,1]*corec_per50 
corec_errdn = (Result.corec_perc[*,1] - Result.corec_perc[*,0])/Result.corec_perc[*,1]*corec_per50 

writecol, outpath+starlabel+'_mantlecomp_final.txt', Result.mantlecompsname, mantlec_per50, mantlec_errup, mantlec_errdn, fmt='(A, F, F, F)'
writecol, outpath+starlabel+'_corecomp_final.txt',  Result.corecompsname, corec_per50, corec_errup, corec_errdn, fmt='(A, F, F, F)' 

print, 'Output CMF to file for '+starlabel
cmf_per50 = Result.cmf_perc[1]
cmf_per50_errup = Result.cmf_perc[2] - cmf_per50
cmf_per50_errdn = cmf_per50 - Result.cmf_perc[0]
writecol, outpath+starlabel+'_fcoremass_final.txt',  cmf_per50, cmf_per50_errup, cmf_per50_errdn, fmt='(F,F, F)'

;if starlable ne 'examstar' then goto, yourstar ;;; This line is used for calling for your sample of stars and may be replicated accordingly.



;stop;---------------------------
TOC
print, 'excecution done!'
stop
END
