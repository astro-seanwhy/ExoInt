import numpy as np
from astropy.table import Table

#======= Constants =========#
chemtable = Table.read('data/atomwttc_new.txt', format='ascii')
chemtable.add_index('elem')
atwt = dict(zip(np.asarray(chemtable['elem']), np.asarray(chemtable['atwt'])))


#======= Function =========#
def chemsysmodel(planetabuN, verbose=False, Si_in_core=True, Nmc=1e5):
    """
    PRE:
        planetabuN  ... dictionary with planetary elemental abundances, nan values allowed and used for 
                        elements which abundance is not known
        verbose     ... if True, then chemsysmodel comments on its progress and results.
    POST:
        mantle      ... dictionary with mantle compositions and their respective mass fractions, e.g.
                        mantle['SiO2'] gives you the mass-fraction of SiO2 in the mantle. Contains following components: SiO2, MgO, CaO, Na2O, Al2O3, FeO, NiO, SO3, CO2, C, O and metals [which may include Si, Mg, Ca, Fe, ... as well as less important elements such as Cr, V, Ti, Sc,...]. The values of the dictionary might be arrays.
        core        ... similar to mantle, but with the core composition. Includes Fe, Ni and S. Values can be
                        arrays.
        fcoremass   ... fraction of the core to the total mass
        Nvalc       ... number of valid draws chemsysmodel makes. If Nvalc==1, then the values of"mantle" and
                        "core" are scalar numbers. If Nvalc>1, "mantle" and "core" contain arrays, and if Nvalc == np.nan, the computation failed and can/should be dismissed.
    """

    #print('planetabuN from chemsysmodel: ', planetabuN)
    for key in planetabuN:
        if planetabuN[key] < 0:
            planetabuN[key] = np.nan   

    Nmc = int(Nmc)

    #============== chemcial network setup ==============#
    mantle = {'Na2O':np.nan, 'CaO':np.nan, 'MgO':np.nan, 'Al2O3':np.nan, 'SiO2':np.nan, 'FeO':np.nan, 'NiO':np.nan, 'SO3':np.nan, 'CO2':np.nan, 'C':np.nan, 'Metals':np.nan, 'ExtraO':np.nan}
    mantle_Nval = {'Na2O':np.nan, 'CaO':np.nan, 'MgO':np.nan, 'Al2O3':np.nan, 'SiO2':np.nan, 'FeO':np.nan, 'NiO':np.nan, 'SO3':np.nan, 'CO2':np.nan, 'C':np.nan, 'Metals':np.nan, 'ExtraO':np.nan}

    if Si_in_core == True:
        core = {'Fe':np.nan, 'Ni':np.nan, 'S':np.nan, 'Si':np.nan}
        core_Nval = {'Fe':np.nan, 'Ni':np.nan, 'S':np.nan, 'Si':np.nan}
    else:
        core = {'Fe':np.nan, 'Ni':np.nan, 'S':np.nan}
        core_Nval = {'Fe':np.nan, 'Ni':np.nan, 'S':np.nan}

    fcoremass=np.nan
    Nvalc = np.nan
    
    
    #========== assert presence of important elements =========#
    # important planet forming elements: O,Si,Mg,Fe
    if np.isnan(planetabuN['O']) or np.isnan(planetabuN['Si']) or np.isnan(planetabuN['Mg']) or np.isnan(planetabuN['Fe']):
        if verbose:
            print("Important elements missing !!")
            print("Calculation is DISMISSED; N/A Result for all are returned")
     
        mantle = np.nan
        core = np.nan
        fcoremass = np.nan
        Nvalc = np.nan
        return mantle, core, fcoremass, Nvalc, mantle_Nval, core_Nval
    
    
    #========== catch and dismiss carbide planets =========#
    if planetabuN['C']/planetabuN['O'] >= 0.8:
        if verbose:
            print("Carbide Planet !!")
            print("It may form C0, SiC, Mg2C, Fe3C, and CaC2 ... beyond the scope of a silicate planet!")
            print("Calculation is DISMISSED; N/A Result for all are returned")
     
        mantle = np.nan
        core = np.nan
        fcoremass = np.nan
        Nvalc = np.nan
        return mantle, core, fcoremass, Nvalc, mantle_Nval, core_Nval
    
    
    #================= MANTLE FORMATION ===================#

    planetOleft = planetabuN['O']

    # ---------- elements that will always fully oxidize --------#
    litho_order = ['Na', 'Ca', 'Mg', 'Al']
    elem_stoech = [2,1,1,2]
    oxy_stoech = [1,1,1,3]
    oxides = ['Na2O', 'CaO', 'MgO', 'Al2O3']
    OtoLithR = planetabuN['O']/np.nansum([planetabuN['Na']/2., planetabuN['Ca'],  planetabuN['Mg'], planetabuN['Al']*3/2])
    if OtoLithR <= 1:
        if verbose:
            print("O/(Na/2+Ca+Mg+Al*3/2) lt 1 -- An invalid MC draw of planet abundances")
            print("Re-draw abundances...")
        planetOleft = 0
        return mantle, core, fcoremass, Nvalc, mantle_Nval, core_Nval
    else:
        for n in range(len(litho_order)):
            elem = litho_order[n]
            compsabu = planetabuN[elem]/elem_stoech[n]
            mantle[oxides[n]] = compsabu*(elem_stoech[n]*atwt[elem] + oxy_stoech[n]*atwt['O'])
            mantle_Nval[oxides[n]] = 1
            if not np.isnan(planetabuN[elem]): planetOleft -= compsabu*oxy_stoech[n]
         
    # ---------- process SiO2; Si can be shared between core & mantle, or reside in a metal phase -----------#
    if planetOleft/2. <= planetabuN['Si']: 
        if verbose: print("O has been used up after forming SiO2") 
        compsabu=planetOleft/2.
        mantle['SiO2'] = compsabu*(atwt['Si'] + 2*atwt['O'])
        if Si_in_core == True:
            core['Si'] = (planetabuN['Si'] - compsabu) * atwt['Si']
            core_Nval['Si'] = 1
        else:
            mantle['Metals'] = (planetabuN['Si'] - compsabu) * atwt['Si']
            mantle['Metals'] = 1
        planetOleft = 0
    else: 
        compsabu = planetabuN['Si']
        mantle['SiO2'] = compsabu*(atwt['Si'] + 2*atwt['O'])
        planetOleft -= compsabu*2
    
    mantle_Nval['SiO2'] = 1




    #=============== CORE FORMATION ================#
    
    if planetOleft>0:
        """
        after oxidizing the mantle, if there is still oxygen left it will oxidize core elements
        """
        
        if planetOleft > np.nansum([planetabuN['Fe'], planetabuN['Ni'], 3*planetabuN['S']]):
            """
            there is sufficient oxygen to fully oxidize the core; we will get a coreless planet
            """

            
            if verbose: print("A coreless planet!")
            mantle['FeO'] = planetabuN['Fe']*(atwt['Fe'] + atwt['O'])
            mantle['NiO'] = planetabuN['Ni']*(atwt['Ni'] + atwt['O'])
            mantle['SO3'] = planetabuN['S']*(atwt['S'] + 3*atwt['O'])
            mantle_Nval['FeO'] = 1
            mantle_Nval['NiO'] = 1
            mantle_Nval['SO3'] = 1
            planetOleft -= np.nansum([planetabuN['Fe'], planetabuN['Ni'], 3*planetabuN['S']])
            
            #Here also CO2 and C
            if not np.isnan(planetabuN['C']):
                if planetOleft>=(planetabuN['C']*2):
                    mantle['CO2'] = planetabuN['C']*(atwt['C'] + 2*atwt['O'])
                    planetOleft -= 2*planetabuN['C']
                    if planetOleft>0:
                        if verbose: print("Remaining oxygen: ", planetOleft)
                        mantle['ExtraO'] = planetOleft * atwt['O']
                        mantle_Nval['ExtraO'] = 1
                else:
                    if verbose: print("No remaining oxygen")
                    planetCO2 = planetOleft/2
                    mantle['CO2'] = planetCO2*(atwt['C'] + 2*atwt['O'])
                    mantle['C'] = (planetabuN['C'] - planetCO2)*atwt['C']
                    mantle_Nval['C'] = 1
                
                mantle_Nval['CO2'] = 1  ## Should C abundance be available, CO2 will be present in this scenario 

            else:
                mantle['ExtraO'] = planetOleft * atwt['O']
                mantle_Nval['ExtraO'] = 1  
            
            core['Fe'] =  np.nan
            core['Ni'] = np.nan
            core['S'] = np.nan
            core['Si'] = np.nan
            fcoremass = 0
            
            #### turn molar mass into massfractions
            mantlemass = np.nansum(np.asarray(list(mantle.values())) )
            for comp in mantle:
                mantle[comp] /= mantlemass
            Nvalc = 1
            
            

        
        else:
            """
            This is the case when there is oxygen to oxidize some core material, but not all of it.
            """
            if verbose: print("Partially oxidise Fe, Ni, and S") 
            N_NiO_mantle = np.random.uniform(size=Nmc)*planetabuN['Ni']
            N_SO3_mantle = np.random.uniform(size=Nmc)*planetabuN['S']
            N_FeO_mantle = np.nansum([planetOleft*np.ones(Nmc), -N_NiO_mantle, -3*N_SO3_mantle], axis=0)
            
            
            #### first validity check: Filter all negative abundances for FeO in mantle
            first_validity_mask = np.where(N_FeO_mantle >= 0)
            N_NiO_mantle = N_NiO_mantle[first_validity_mask]
            N_SO3_mantle = N_SO3_mantle[first_validity_mask]
            N_FeO_mantle = N_FeO_mantle[first_validity_mask]
                        
            N_Ni_core = planetabuN["Ni"] - N_NiO_mantle
            N_S_core = planetabuN["S"] - N_SO3_mantle
            N_Fe_core = planetabuN["Fe"] - N_FeO_mantle
            
            Nvalc = first_validity_mask[0].size

            #### second validity check: Filter all core abundances where iron is not in balance with nickel
            if not np.isnan(planetabuN['Ni']):
                second_validity_mask = np.where(abs(N_Fe_core/N_Ni_core - 18) <= 4)
                
                N_NiO_mantle = N_NiO_mantle[second_validity_mask]
                N_SO3_mantle = N_SO3_mantle[second_validity_mask]
                N_FeO_mantle = N_FeO_mantle[second_validity_mask]
                                
                N_Ni_core = N_Ni_core[second_validity_mask]
                N_S_core = N_S_core[second_validity_mask]
                N_Fe_core = N_Fe_core[second_validity_mask]
            
                Nvalc = second_validity_mask[0].size
            
            #### third validity check: There should not be more sulfur than iron
            if not np.isnan(planetabuN['S']):
                third_validity_mask = np.where(N_Fe_core>=N_S_core)
                
                N_NiO_mantle = N_NiO_mantle[third_validity_mask]
                N_SO3_mantle = N_SO3_mantle[third_validity_mask]
                N_FeO_mantle = N_FeO_mantle[third_validity_mask]
                                
                N_Ni_core = N_Ni_core[third_validity_mask]
                N_S_core = N_S_core[third_validity_mask]
                N_Fe_core = N_Fe_core[third_validity_mask]
            
                Nvalc = third_validity_mask[0].size
            
            if verbose: print('Numbers of valid guess for Fe and Ni in the mantle and the core: ', Nvalc, ' of', Nmc)

            if Nvalc==0:
                #case where the Monte Carlo draws return zero valid core models:
                mantle = np.nan
                core = np.nan
                fcoremass = np.nan
                for key in mantle_Nval:
                    mantle_Nval[key] = np.nan
                for key in core_Nval:
                    core_Nval[key] = np.nan
                return mantle, core, fcoremass, Nvalc, mantle_Nval, core_Nval
            else:
                #case where the Monte Carlo draws return at least one valid core model:
                 
                #### turn the molar masses in the "mantle" dictionary into arrays
                for comp in mantle: 
                    mantle[comp] *= np.ones(Nvalc)
                
                
                mantle['FeO'] = N_FeO_mantle * (atwt['Fe'] + atwt['O'])
                mantle['NiO'] = N_NiO_mantle * (atwt['Ni'] + atwt['O'])
                mantle['SO3'] = N_SO3_mantle * (atwt['S'] + 3*atwt['O'])
                
                core['Fe'] = N_Fe_core * atwt['Fe']
                core['Ni'] = N_Ni_core * atwt['Ni']
                core['S'] = N_S_core * atwt['S']
                if Si_in_core == True: core['Si'] = np.ones_like(N_S_core)*core['Si']
                


                ## Here there will be only graphite, no CO2, since O has been used up by Fe, Ni, and S in this case
                mantle['C'] = planetabuN['C']*atwt['C'] * np.ones(Nvalc)
                mantle_Nval['C'] = 1


                #============= mean of every array ==================#
                # The output of ExoInt would be dominated by the MC draws if there is no summary
                # statistics computed on them, and the results would be wrong. Hence, we compute the
                # mean of every compound, but the median would work as well and would return similar
                # results.

                for mantlecompound in mantle.keys():
                    if not np.any(np.isnan(mantle[mantlecompound])):
                        #mantle[mantlecompound] = np.nanmedian( mantle[mantlecompound] )
                        mantle[mantlecompound] = np.nanmean( mantle[mantlecompound] )
                    else:
                        mantle[mantlecompound] = np.nan
                for corecompound in core.keys():
                    if not np.any(np.isnan(core[corecompound])):
                        #core[corecompound] = np.nanmedian( core[corecompound] )
                        core[corecompound] = np.nanmean( core[corecompound] )
                    else:
                        core[corecompound] = np.nan


               

                #======= turn molar mass into massfractions ==========#
                mantlemass = np.nansum(np.asarray(list(mantle.values())), axis=0 )
                coremass = np.nansum(np.asarray(list(core.values())), axis=0 )


               
                for comp in mantle:
                    mantle[comp] /= mantlemass
                for comp in core:
                    core[comp] /= coremass
                fcoremass = coremass/(mantlemass + coremass)


                #===== store occurance rates =======#
                mantle_Nval['FeO'] = 1
                mantle_Nval['NiO'] = 1
                mantle_Nval['SO3'] = 1
                core_Nval['Fe'] = 1
                core_Nval['Ni'] = 1
                core_Nval['S'] = 1

 
                
    else:
        # here is the case were all oxygen has been used up to form the mantle
                
        core['Fe'] =  planetabuN['Fe'] * atwt['Fe']
        core['Ni'] = planetabuN['Ni'] * atwt['Ni']
        core['S'] = planetabuN['S'] * atwt['S']

        core_Nval['Fe'] = 1
        core_Nval['Ni'] = 1
        core_Nval['S'] = 1

        if np.isnan(planetabuN['Ni']) != True:
            # constrain Fe/Ni to 18 +- 4
            if planetabuN['Fe']>22*planetabuN['Ni']:
                # case where we have too much iron; overabundant iron is transfered into the mantle
                core['Fe'] = planetabuN['Ni']*22 * atwt['Fe']
                mantle['Metals'] = np.nansum([mantle['Metals'], planetabuN['Fe']*atwt['Fe'] - core['Fe']])
                mantle_Nval['Metals'] = 1
            elif planetabuN['Fe']<14*planetabuN['Ni']:
                # case where we have too much nickel; overabundant nickel is transfered to the mantle
                core['Ni'] = planetabuN['Fe']/14*atwt['Ni']
                mantle['Metals'] = np.nansum([mantle['Metals'], planetabuN['Ni']*atwt['Ni'] - core['Ni'] ])
                mantle_Nval['Metals'] = 1
            else:
                # case where Fe/Ni is between 18 +- 4; don't do anything
                pass
        else:
            pass

        mantle['C'] = planetabuN['C']*atwt['C']
        mantle_Nval['C'] = 1




        #### turn molar mass into massfractions
        mantlemass = np.nansum(np.asarray(list(mantle.values())) )
        coremass = np.nansum(np.asarray(list(core.values())) )
        for comp in mantle:
            mantle[comp] /= mantlemass
        for comp in core:
            core[comp] /= coremass
        fcoremass = coremass/(mantlemass + coremass)

        Nvalc = 1



    for key in mantle:
        if np.isnan(mantle[key]):
            mantle_Nval[key] = np.nan
        else:
            mantle_Nval[key] - 1

    for key in core:
        if np.isnan(core[key]):
            core_Nval[key] = np.nan
        else:
            core_Nval[key] = 1
    
    return mantle, core, fcoremass, Nvalc, mantle_Nval, core_Nval
        
        
        
        
