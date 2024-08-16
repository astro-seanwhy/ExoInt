"""
This program evaluates the mineralogy for a single planet, at 
a specified chemical composition.

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
from Utilities.beautiful_figures import fancy_mineralogy_plot
import matplotlib.pyplot as plt
import os
import sys
from astropy.table import Table

##cleanup any previous temporary running files from PerPlex, if exist
temppath = 'Utilities/tmp_files/perpleX_tmp'
if os.path.isdir(temppath):
    os.system('rm -r ' + temppath)

##make mineral_output folder if not exist yet.
outpath = 'mineral_output/'
if not os.path.isdir(outpath):
    os.system("mkdir -p "+str(outpath))
    
datapath = '../output/'  # the data path to the stoichiometric compositional output of pyExoInt; adapt it, if necessary.

MRdata = Table.read(sys.argv[1], format='ascii')
planet_labels = MRdata['planet'].tolist()
planet_R = MRdata['R'].tolist()  ##please note that the radius range should be restrined within [0.5, 1.5] and it is best tested for Earth-sized planets.
planet_M = MRdata['M'].tolist()  ##mass is not directly used for simulation, but you should ensure that it is less than 10 Earth_mass and the planet is, in principle, rocky. 


###calculate the 'mean' results with the 'mean' input
for sname in planet_labels:
    print(r'Run for '+sname+'...')

    tmp = Table.read(datapath+sname+'_mantlecomp_final.txt', format='ascii')
    nSiO2 = np.where(tmp['param'] == 'SiO2')
    nMgO = np.where(tmp['param'] == 'MgO')
    nFeO = np.where(tmp['param'] == 'FeO')
    nAl2O3 = np.where(tmp['param'] == 'Al2O3')
    nCaO = np.where(tmp['param'] == 'CaO')
    nNa2O = np.where(tmp['param'] == 'Na2O')
    mantlesum = np.nansum([float(tmp['median'][nSiO2]), float(tmp['median'][nMgO]), float(tmp['median'][nFeO]), float(tmp['median'][nAl2O3]), float(tmp['median'][nCaO]), float(tmp['median'][nNa2O])])
    mantle_chemsys = {'SiO2':float(tmp['median'][nSiO2])/mantlesum, 
                  'MgO':float(tmp['median'][nMgO])/mantlesum,
                  'FeO':float(tmp['median'][nFeO])/mantlesum,
                  'Al2O3':float(tmp['median'][nAl2O3])/mantlesum,
                  'CaO':float(tmp['median'][nCaO])/mantlesum,
                  'Na2O':float(tmp['median'][nNa2O])/mantlesum }
    
    tmp = Table.read(datapath+sname+'_corecomp_final.txt', format='ascii')
    nFe = np.where(tmp['param'] == 'Fe')
    nNi = np.where(tmp['param'] == 'Ni')
    nSi = np.where(tmp['param'] == 'Si')
    nS = np.where(tmp['param'] == 'S')
    coresum = np.nansum([float(tmp['median'][nFe]), float(tmp['median'][nNi]), float(tmp['median'][nSi]), float(tmp['median'][nS])])     
    core_chemsys = {'Fe':float(tmp['median'][nFe])/coresum,
                    'Ni':float(tmp['median'][nNi])/coresum, 
                    'Si':float(tmp['median'][nSi])/coresum, 
                    'S': float(tmp['median'][nS])/coresum }
    tmp = Table.read(datapath+sname+'_fcoremass_final.txt', format='ascii')
    CMF = float(tmp['median'])
    #breakpoint()

    
    # ======evaluate the interior structure model
    R = planet_R[planet_labels.index(sname)]
    M = planet_M[planet_labels.index(sname)] ##Please note that M (mass) is not used for the simulation, but for comparing with the modeled mass only.
    P = Planet(R, CMF, mantle_chemsys, core_chemsys)
    P._Delta_T_CMB = "NL20"

    P.solve_structure() 
    print('M/Mearth: ', P.M) ##in Mearth = 5/974e24 kg
    f = open(outpath+'runpara_R_M_'+sname+'.txt', "w")
    f.writelines(['R: ', str(R), '\nM_model: ', str(P.M), '\nM_data: ', str(M), '\nCRF: ', str(P.CRF)]) #, '\nDeltaT: ', str(P.Delta_T_CMB)])
    f.close()

    #Generate the mineralogy plot
    fig, ax, ax2 = fancy_mineralogy_plot(P) #, refRhofile=True, PREM=True)
    plt.savefig(outpath+'mineralogy_'+sname+'.pdf')
  
    os.system('cp Utilities/mineralogy.tab '+outpath+'mineralogy_'+sname+'.txt')

    #breakpoint()

