import numpy as np
from astropy.table import Table
from chemsysmodel import chemsysmodel


def load_data(filename, asymerr=False):
    """
    Loads the stellar abundances from 'filename' and returns the
    stellar abundances (and its upper and lower uncertainties)
    in log-units relative to the solar abundances (dex)

    PRE:
        filename        ... file containing the stellar abundances
        asymerr         ... (optional) True if the model shall use upper and lower
                            uncertainties; if False, the lower uncertainty will be
                            equal to the upper (symmetric error). Default is False.

    POST:
        stardex         ... stellar abundances, in dex
        stardexerr      ... stellar abundances uncertainty, in dex; if asymerr==True,
                            it will denote the upper uncertainties
        stardexerrdn    ... stellar abundances uncertainty, in dex; if asymerr==True,
                            it will denote the lower uncertainties; if asymerr==False,
                            it will be a copy of the upper uncertainties
    """

    tmp_table = Table.read(filename, format='ascii')

    stardex = np.asarray(tmp_table["abu_[X/H]"])
    stardexerr = np.asarray(tmp_table["err_abu_[X/H]"])
        
    if asymerr == True:
        stardexerrdn = np.asarray(tmp_table["errdn_abu_[X/H]"])
    else:
        stardexerrdn = stardexerr
          
    return stardex, stardexerr, stardexerrdn
    
    

def produce_planetabuN(stardex, stardexerr, stardexerrdn, refsolar='A21'):
   
    """
    PRE:
        stardex         ... array with log[X/H] for every element H to U, np.nan allowed
        stardexerr      ... array with err of stardex, np.nan values allowed
        stardexerrdn    ... in the case of asymetric error bars, this corresponds to the lower error bar of stardex, np.nan values allowed
        planetabuN      ... (optional) bool, True if 'stardex' is already a properly normalized planetary bulk chemistry; default is False
        refsolar        ... (optional) switch between differen solar abundances; options are:
                                A21     Asplund et al. 2021
                                A09     Asplund et al. 2009
                                L09     Lodders et al. 2009
                                W19     Wang et al 2019
                            default is 'A21'
        devol           ... (optional) bool, True if the abundances in 'stardex' shall be devolatilised according to the model of Wang et al. 2019
                            Set to 'False' if the abundances in 'stardex' are already devolatilised. Default is 'True'.

    POST:
        mantle_best ... dictionary with the best estimate for mantle compositon; each value is scalar
        core_best   ... same as mantle_best, but for the core
        fcoremass_best  ... best estimate of coremassfraction
        mantle      ... dictionary with mantle compositon; each entry is an array
        core        ... same as mantle, but for the core
        fcoremass   ... array with estimates for the coremassfraction
    """
    
    #=================== Input ==================#
    Nelems=83
      



    
    #================== Preparation ======================#

    tmp = Table.read('data/atomwttc_new.txt', format='ascii')
    atomid=np.asarray(tmp['atom'])
    elemid=np.asarray(tmp['elem'])
    atomwt=np.asarray(tmp['atwt'])
    elemtc=np.asarray(tmp['50%Tc'])

    tmp = Table.read('data/protosunppmwhy.txt', format='ascii')
    solarppm=np.asarray(tmp['abu/ppm'])
    solarppmerrup=np.asarray(tmp['abuerrup/ppm'])
    solarppmerrdn=np.asarray(tmp['abuerrdn/ppm'])

    tmp = Table.read('data/PEppmwhy.txt', format='ascii')
    bulkearthppm=np.asarray(tmp['abu/ppm'])
    bulkearthppmerrup=np.asarray(tmp['errUP/ppm'])
    bulkearthppmerrdn=np.asarray(tmp['errDN/ppm'])
    
    tmp = Table.read('data/Asplund09dex.txt', format='ascii')
    A09solardex=np.asarray(tmp['abu/dex'])
    A09solardexerr=np.asarray(tmp['err/dex'])
    A09solardex[np.where(A09solardex==-11)] = np.nan
    A09solardexerr[np.where(A09solardexerr==-11)] = np.nan        

    tmp = Table.read('data/Asplund2021dex.txt', format='ascii')
    A21solardex=np.asarray(tmp['abu/dex'])
    A21solardexerr=np.asarray(tmp['err/dex'])
    A21solardex[np.where(A21solardex==-11)] = np.nan
    A21solardexerr[np.where(A21solardexerr==-11)] = np.nan
    
    tmp = Table.read('data/Lodd09-origin.txt', format='ascii')
    L09solardex=np.asarray(tmp['recommabu'])
    L09solardexerr=np.asarray(tmp['recommerr'])
      
    tmp= Table.read('data/protosunwhy.txt', format='ascii')
    W19solardex=np.asarray(tmp['abu/dex'])
    W19solardexerr=np.asarray(tmp['abuerrup/dex'])
    W19solardexerrdn=np.asarray(tmp['abuerrdn/dex'])
    
    tmp = Table.read('data/devolmodel_lin_pub.txt', format='ascii')
    ymodelN=np.asarray(tmp['col2']); ymodelNsdup=np.asarray(tmp['col3'])
    ymodelNsddn=np.asarray(tmp['col4'])

    tmp = Table.read('data/devolmodel_log_pub.txt', format='ascii')
    ymodellog=np.asarray(tmp['col2'])
    ymodeldexsdup=np.asarray(tmp['col3'])
    ymodeldexsddn=np.asarray(tmp['col4'])
    del tmp
    
    
    ### replace the -11 data with np.nan
    A09solardex[np.where(A09solardex==-11)] = np.nan
    A09solardexerr[np.where(A09solardexerr==-11)] = np.nan
    A21solardex[np.where(A21solardex==-11)] = np.nan
    A21solardexerr[np.where(A21solardexerr==-11)] = np.nan
    stardex[np.where(stardex==-11)] = np.nan
    stardex[np.where(stardex==-11)] = np.nan
    ymodelN[np.where(ymodelN==-11)] = np.nan
    ymodelNsdup[np.where(ymodelNsdup==-11)] = np.nan
    ymodelNsddn[np.where(ymodelNsddn==-11)] = np.nan
    ymodellog[np.where(ymodellog==-11)] = np.nan
    ymodeldexsdup[np.where(ymodeldexsddn==-11)] = np.nan
    ymodeldexsddn[np.where(ymodeldexsddn==-11)] = np.nan
    
    #get array-indices for the important elements:
    for i in range (Nelems-1):
        if elemid[i] == 'C': nC=i
        if elemid[i] == 'O': nO=i
        if elemid[i] == 'S': nS=i
        if elemid[i] == 'Hg': nHg=i
        if elemid[i] == 'Na': nNa=i
        if elemid[i] == 'Si': nSi=i
        if elemid[i] == 'Fe': nFe=i
        if elemid[i] == 'Mg': nMg=i
        if elemid[i] == 'Ca': nCa=i
        if elemid[i] == 'Ti': nTi=i
        if elemid[i] == 'Al': nAl=i
        if elemid[i] == 'Cr': nCr=i
        if elemid[i] == 'Mn': nMn=i
        if elemid[i] == 'Ni': nNi=i
        if elemid[i] == 'K': nK=i
        if elemid[i] == 'P': nP=i


    #=========== process stardex to get planet bulk chemistry ==========#

    #---- choose the reference solar abundances; Asplund+2021 is the default reference
    if refsolar == 'A21':
        solardex = A21solardex
        solardexerr = A21solardexerr
        solardexerrdn = A21solardexerr    

    elif refsolar=='A09':
        solardex = A09solardex
        solardexerr = A09solardexerr
        solardexerrdn = A09solardexerr    
    elif refsolar == 'L09':
        solardex = L09solardex
        solardexerr = L09solardexerr
        solardexerrdn = L09solardexerr
    elif refsolar == 'W19':
        solardex = W19solardex
        solardexerr = W19solardexerr
        solardexerrdn = W19solardexerrdn
    elif refsolar == 'custom':
        tmp = Table.read('data/refsolar_dex.txt', format='ascii')
        selfsolardex=np.asarray(tmp['abu/dex'])
        selfsolardexerr=np.asarray(tmp['err/dex'])
        del tmp
        selfsolardex[np.where(A09solardex==-11)] = np.nan
        selfsolardexerr[np.where(A09solardexerr==-11)] = np.nan
        solardex = selfsolardex
        solardexerr = selfsolardexerr
        solardexerrdn = selfsolardexerr        
    else:
        raise NotImplementedError("Invalid keyword for the reference solar abundance 'refsolar'; please choose one from 'A09'- Asplund et al. 2009, 'l09' - Lodders et al. 2009, and 'W19' - Wang et al. 2019; Or otherwise leave out the keyword, and the default (recommended) Asplund et al. 2021 solar abundance will be adopted; 'custom' -- your customed reference solar abundances" )

    #------convert the differential abundances to the abunance in dex with solar abundances considered
    starabudex = stardex + solardex
    starabudexerr_up = np.sqrt(stardexerr**2 + solardexerr**2)
    starabudexerr_dn = np.sqrt(stardexerr**2 + solardexerrdn**2)
    
    #----- choose reference element ------#
    if not np.isnan(starabudex[nAl]):
        nref = nAl
    elif not np.isnan(starabudex[nCa]):
        nref = nCa
    else:
        nref = nFe        
    
    #------ convert stellar abu dex into number of atoms --------#
    starabuN = 10**(starabudex - starabudex[nref]) * 100 #use nref as reference metal, and normalize the abundancies to 100 atoms of the reference element
    starabuNerr_up = 10**(starabudexerr_up - 1)*starabuN
    starabuNerr_dn = (1-10**(-starabudexerr_dn))*starabuN
    

    #----- devolatilise -------#
    planetabuN = starabuN * ymodelN
    planetabuNerr_up = np.sqrt(starabuNerr_up**2.*ymodelN**2. + ymodelNsdup**2.*starabuN**2.)
    planetabuNerr_dn = np.sqrt(starabuNerr_dn**2.*ymodelN**2. + ymodelNsddn**2.*starabuN**2.)
    
    
    return dict(zip(elemid, planetabuN)), dict(zip(elemid, planetabuNerr_up)), dict(zip(elemid, planetabuNerr_dn))


def single_simulation_report(planetabuN):
    
    mantle, core, fcoremass, Nvalc, mantle_Nval, core_Nval = chemsysmodel(planetabuN)
    print("------------- MANTLE (wt%) ------------- ")
    for elem in mantle:
        if np.any(~np.isnan(mantle[elem])):
            mantle[elem] = np.nanmean(mantle[elem])
            print(elem, ':', np.round(mantle[elem]*100, 2))
        else:
            mantle[elem] = np.nan
    print("-------------- CORE (wt%) -------------- ")
    for elem in core:
        if np.any(~np.isnan(core[elem])):
            core[elem] = np.nanmean(core[elem])
            print(elem, ':', np.round(core[elem]*100, 2) )
        else:
            core[elem] = np.nan
    print("------- CORE MASS FRACTION (wt%) ------- ")
    fcoremass = np.nanmean(fcoremass)
    print("CMF: ", np.round(fcoremass*100, 2))
    #print("Nvalc: ", Nvalc, '\n')
    
    



def generate_asymmetric_random(mu, sigma_up, sigma_dn):
    p1 = np.random.rand()
    if p1<0.5:
        p2 = np.random.normal(loc=0, scale=sigma_dn)
        return mu-abs(p2)
    else:
        p2 = np.random.normal(loc=0, scale=sigma_up)
        return mu+abs(p2)
