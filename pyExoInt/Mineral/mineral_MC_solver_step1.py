"""
This program performs a Monte-Carlo simulation of exoplanet mantle
mineralogy.

Copyright (C) 2021-2023  Haiyang S. Wang, Fabian L. Seidler

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as np
from Utilities import Planet
from multiprocessing import Pool
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas
from scipy.interpolate import UnivariateSpline, interp1d
from tqdm import tqdm
import os
import sys
from astropy.table import Table
import signal
from contextlib import contextmanager
import math

#====path setup =====#
datapath = '../output/' # the stoichiometric compositional output of pyExoInt; adapt it, if necessary.

outpath = 'mineral_output/MC_files/'
if not os.path.isdir(outpath):
    os.system('mkdir -p '+ outpath)

#====define a function for MC draw on data with asymetric error bars ====#
def random_asymmetric(mu, sigma_up, sigma_dn):
    p1 = np.random.rand()
    if p1<0.5:
        p2 = np.random.normal(loc=0, scale=sigma_dn)
        return mu-abs(p2)
    else:
        p2 = np.random.normal(loc=0, scale=sigma_up)
        return mu+abs(p2)

def ordinal(n):
    if n%10==1: ordsign = 'st'
    elif n%10==2: ordsign = 'nd'
    elif n%10==3: ordsign = 'rd'
    else: ordsign = 'th'
    return ordsign


class TimeoutException(Exception): pass
@contextmanager
def time_limit(seconds):
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)

#===== Function ====#
def f(n): #R, CMF, mantle_chemsys, core_chemsys):
    """
    Evaluates one interior structure model from a random draw
    of mantle and core compostion

    Parameters
    ----------
        n : int
            index of thread (via pool.map)
    """
 
    global nrun   ##globalise the i ordinal number of the run
    np.random.seed(n)

    ID = str(np.random.randint(int(1e9)))

    with open(logpath+str(nrun)+"_"+ID+".log", 'w') as logfile:
        logfile.write("starting the process of " + str(ID) + "\n")


        ##--Mantle
        logfile.write("Mantle \n")
        mantle_chemsys_array = np.asarray([np.nan]*Nmantle)
        while True:
            for i in range(Nmantle): mantle_chemsys_array[i] = random_asymmetric(mantle_loc[i], mantle_scale_up[i], mantle_scale_dn[i])
            if np.all(mantle_chemsys_array>=0):
                #print('mantle_chemsys: ', mantle_chemsys_array)
                break
        
        mantle_chemsys = dict(zip(mantlechemsys_names, mantle_chemsys_array))
      
       
        logfile.write("CORE \n")
        ##--Core    
        core_chemsys_array = np.asarray([np.nan]*Ncore)
        while True:
            for i in range(Ncore): core_chemsys_array[i] = random_asymmetric(core_loc[i], core_scale_up[i], core_scale_dn[i])
            if np.all(core_chemsys_array>=0):
                #print('core_chemsys: ', core_chemsys_array)
                break
     
        core_chemsys = dict(zip(corechemsys_names, core_chemsys_array))
     
        
        logfile.write("CMF \n")
        ##--CMF
        while True:
            CMF = random_asymmetric(cmf_loc, cmf_scale_up, cmf_scale_dn) 
            if CMF>=0 and CMF<1:
                #print('CMF: ', CMF)
                break
    
        ##--R
        while True:
            R = np.random.normal(loc=R_loc, scale=R_scale)  #IF THE ERRORS ON 'R' ARE ASYMMETRIC, YOU CAN ALSO GENERATE THE ASYMMETRIC MC DRAW BY REFERRING TO THE EXAMPLE ABOVE FOR 'CMF'
            if R>=0 and R<1.6:  ## 2 is a conservative upper limit of the radius of a likely rocky planet
               # print('R: ', R)
                break
   
       
        # evaluate the interior structure model
        #--------------------------------------
        logfile.write("computing mineralogy  \n")
        print("computing mineralogy for " + str(n) + ordinal(n) + " MC draw ...")
        
        df2 = [] ##define an empty df2, in case of timeout
        #Ntimeout = 0
        try:
            with time_limit(180):
                P = Planet(R, CMF, mantle_chemsys, core_chemsys)
                P._Delta_T_CMB = "NL20"
    
                logfile.write("solve structure  \n")

                subpath_name = str(nrun)+'_'+str(n)               # where Perple_X stores temporary files
                P.solve_structure(subpath_name)

 
                logfile.write("store files  \n")
                # store in pandas DataFrame
                radius_array = np.append(P.R_core, P.R_mantle)
                pressure_array = np.append(P.P_core, P.P_mantle)
                temperature_array = np.append(P.T_core, P.T_core[-1])
                temperature_array = np.append(temperature_array, P.T_mantle)
                density_array = np.append(P.RHO_core, P.RHO_core[-1])
                density_array = np.append(density_array, P.RHO_mantle)
                density_array = np.append(density_array, P.RHO_mantle[-1])
                
                
                data = np.vstack((radius_array, pressure_array, temperature_array, density_array)).T
                
                
                    # interpolate minerals
                mineral_array = np.asarray([[np.nan]*len(P.mineral_phases)]*(P.N_mantle+1))  #np.zeros((P.N_mantle+1, len(P.mineral_phases)))
                pressure = np.array(P.mantle_dict["pressure"]) * 1e5
                logfile.write("final for loop  \n")
                for i in range(P.mineralogy.shape[1]):
                     mineralfit = interp1d(
                        pressure,
                        P.mineralogy[:, i],
                        kind="slinear",
                        bounds_error=False,
                        fill_value="extrapolate",
                        )
                     mineral = mineralfit(P.P_mantle)
                        #mineral = np.flip(mineral, axis=0)
                     mineral_array[:,i] = mineral
                
                mineral_array = np.vstack((np.zeros((P.N_core+1, len(P.mineral_phases))), mineral_array))
                #mineral_array = np.vstack((np.asarray([[np.nan]*len(P.mineral_phases)]*(P.N_core+1)), mineral_array))
                data = np.hstack((mineral_array, data))
                
                logfile.write("write to the file  \n")
                columns = np.hstack((P.mineral_phases, ['radius', 'pressure', 'temperature', 'density']))
                df2 = pandas.DataFrame(data, columns=columns)
                df2.to_csv(samplepath+str(nrun)+'_'+ID+'.csv')
                
                #print("phasename: ", P.mineral_phases)
                #print("mineral_array: ", mineral_array)

        except TimeoutException as skiperr:
            Ntimeout = Ntimeout + 1
            print(str(Ntimeout) + ordinal(Ntimeout) + " TimeOut when computing mineralogy/structure!")
            logfile.write("Timed out when computing mineralogy/structure! NO OUTPT for this process\n")
    
    # print('df2', df2)
    return df2


#======== Main =======#
if __name__ == '__main__':
    ##---Instructions---:
    ## please provide the input MRfile from the command line, as sys.argv[1]
    ## by default, the MRfile is stored in the "MRfiles" folder, in which an example file is given.
    
    MRdata = Table.read(sys.argv[1], format='ascii')
    planet_labels = MRdata['planet'].tolist()
    planet_R = MRdata['R'].tolist()
    planet_Rerr = MRdata['Rerr'].tolist()

    num_files = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[2].isdigit() else 100  ##the number of MC files to be generated (100, by default); the final number of files may be larger than the given one, since the round is towards the bigger side, if the rounding result is not an integer.
    num_cores = int(sys.argv[3]) if len(sys.argv) > 3 and sys.argv[3].isdigit() else 25  ##the number of cores to be used (25, by default) 

    if num_files < num_cores:
        print("the number of files, num_files, is < the number of cores, num_cores; num_files is increased to be equivalent to num_cores")
        num_files = num_cores
    
    Nrun = math.ceil(num_files / num_cores) #round te number of runs -- of each run generate n_cores number workers -- towards a bigger side, if not integer. 
    num_workers = num_cores 
   
    print("For each case, generate (approximate) " + str(num_files) + " files with " + str(num_cores) + " cores")

   

    tmpdir = './Utilities/tmp_files'
    for sname in planet_labels:
        print(r"Run for " + sname + "...")

        ##--create a folder, if not available, for an individual under the outpath ('mineral_output/MC_files/')
        samplepath = outpath + sname+'/'
        logpath = samplepath + 'logs/'
        if not os.path.isdir(samplepath):
            os.system("mkdir -p "+str(samplepath))
        if not os.path.isdir(logpath):
            os.system('mkdir -p '+ str(logpath))  ##put log files into the subdirectory 

        ##--always clean up the temporary files of PerPlex upon running for different cases
        if not os.path.isdir(tmpdir):
            os.system('mkdir ' + tmpdir)
        else:
            os.system('rm -rf ' + tmpdir)
            os.system('mkdir ' + tmpdir)

        ##--Mantle composition from the stoichemtric output of ExoInt
        tmp = Table.read(datapath+sname+'_mantlecomp_final.txt', format='ascii')
        mantlechemsys_names = ['SiO2', 'MgO', 'FeO', 'Al2O3', 'CaO', 'Na2O']
        Nmantle = len(mantlechemsys_names)

        nSiO2 = np.where(tmp['param'] == 'SiO2')
        nMgO = np.where(tmp['param'] == 'MgO')
        nFeO = np.where(tmp['param'] == 'FeO')
        nAl2O3 = np.where(tmp['param'] == 'Al2O3')
        nCaO = np.where(tmp['param'] == 'CaO')
        nNa2O = np.where(tmp['param'] == 'Na2O')
        mantlesum = np.nansum([float(tmp['median'][nSiO2]), float(tmp['median'][nMgO]), float(tmp['median'][nFeO]), float(tmp['median'][nAl2O3]), float(tmp['median'][nCaO]), float(tmp['median'][nNa2O])])

        mantle_loc = [float(tmp['median'][nSiO2]), float(tmp['median'][nMgO]), float(tmp['median'][nFeO]), float(tmp['median'][nAl2O3]), float(tmp['median'][nCaO]), float(tmp['median'][nNa2O])] / mantlesum

            #scale = [max(float(tmp['errup'][nSiO2]), float(tmp['errdn'][nSiO2])), max(float(tmp['errup'][nMgO]),float(tmp['errdn'][nMgO])), max(float(tmp['errup'][nFeO]), float(tmp['errdn'][nFeO])),  max(float(tmp['errup'][nAl2O3]), float(tmp['errdn'][nAl2O3])), max(float(tmp['errup'][nCaO]), float(tmp['errdn'][nCaO])), max(float(tmp['errup'][nNa2O]), float(tmp['errdn'][nNa2O]))] / mantlesum
        mantle_scale_up = [float(tmp['errup'][nSiO2]), float(tmp['errup'][nMgO]), float(tmp['errup'][nFeO]), float(tmp['errup'][nAl2O3]), float(tmp['errup'][nCaO]), float(tmp['errup'][nNa2O])] / mantlesum
        mantle_scale_dn = [float(tmp['errdn'][nSiO2]), float(tmp['errdn'][nMgO]), float(tmp['errdn'][nFeO]), float(tmp['errdn'][nAl2O3]), float(tmp['errdn'][nCaO]), float(tmp['errdn'][nNa2O])] / mantlesum

        mantle_loc[np.isnan(mantle_loc)] = 0
        mantle_scale_up[np.isnan(mantle_scale_up)] = 0
        mantle_scale_dn[np.isnan(mantle_scale_dn)] = 0

        ##--Core composition from the stoichemtric output of ExoInt
        tmp = Table.read(datapath+sname+'_corecomp_final.txt', format='ascii')
        corechemsys_names = ['Fe', 'Ni', 'Si', 'S']
        Ncore = len(corechemsys_names)

        nFe = np.where(tmp['param'] == 'Fe')
        nNi = np.where(tmp['param'] == 'Ni')
        nSi = np.where(tmp['param'] == 'Si')
        nS = np.where(tmp['param'] == 'S')
        coresum = np.nansum([float(tmp['median'][nFe]), float(tmp['median'][nNi]), float(tmp['median'][nSi]), float(tmp['median'][nS])])

        core_loc = [float(tmp['median'][nFe]), float(tmp['median'][nNi]), float(tmp['median'][nSi]), float(tmp['median'][nS])] / coresum
            #scale = [max(float(tmp['errup'][nFe]), float(tmp['errdn'][nFe])), max(float(tmp['errup'][nNi]), float(tmp['errdn'][nNi])), max(float(tmp['errup'][nSi]), float(tmp['errdn'][nSi])), max(float(tmp['errup'][nS]), float(tmp['errdn'][nS]))] / coresum
             
        core_scale_up = [float(tmp['errup'][nFe]), float(tmp['errup'][nNi]), float(tmp['errup'][nSi]), float(tmp['errup'][nS])] / coresum 
        core_scale_dn = [float(tmp['errdn'][nFe]), float(tmp['errdn'][nNi]), float(tmp['errdn'][nSi]), float(tmp['errdn'][nS])] / coresum        

        core_loc[np.isnan(core_loc)] = 0
        core_scale_up[np.isnan(core_scale_up)] = 0
        core_scale_dn[np.isnan(core_scale_dn)] = 0
            
        ##--Core mass fraction from the stoichemtric output of ExoInt
        tmp = Table.read(datapath+sname+'_fcoremass_final.txt', format='ascii')
        cmf_loc = float(tmp['median'])
        #scale = max(float(tmp['errup']), float(tmp['errdn']))
        cmf_scale_up = float(tmp['errup'])
        cmf_scale_dn = float(tmp['errdn'])
  
        ##--Planet radius 
        R_loc  = planet_R[planet_labels.index(sname)]
        R_scale  = 0.  #error on R has been set to be 0; otherwise, the MC plotting will average results (ad different various) incorrectly. A known issue to be fixed in the future

##=====multi process====##
        Ntimeout = 0
        for i in range(Nrun):
            nrun = i+1
            print("Run for " + sname + " for " +str(nrun) + ordinal(nrun) + " round >>>")
            print("by multiprocessing of "+ str(num_workers) + " worker processes on "+ str(num_cores)+" cores...")
            with Pool(processes = num_cores) as pool:
                result = pool.map(f, range(num_workers), chunksize=1)
            
        Nfiles = num_workers*Nrun - Ntimeout
        print("There were " + str(Ntimeout) + " TimeOut!")
        print(str(Nfiles) + " MC output files generated for uncertainty assessement for " + sname)
        
        #breakpoint()




