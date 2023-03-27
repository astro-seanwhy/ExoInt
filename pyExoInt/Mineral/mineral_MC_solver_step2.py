"""
This program plots the output of a mineralgoy Monte-Carlo simulation.

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
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
from scipy.interpolate import UnivariateSpline, interp1d
from scipy.integrate import quad
from scipy.stats import gaussian_kde
#from sklearn.neighbors import KernelDensity
from tqdm import tqdm
import os
import sys
from astropy.table import Table 

#==== Sample (please change the following two lines accordingly) =====#
#planet_labels = ['examplestar1', 'examplestar2']  ## PLEASE MODIFY THESE LABELS TO FIT TO YOUR INPUT SAMPLE
#planet_R = [1.0, 1.0] #in Rearth

MRdata = Table.read(sys.argv[1], format='ascii')
planet_labels = MRdata['planet'].tolist()
planet_R = MRdata['R'].tolist()

##===Evaluate the uncertainties upon the output files of Step 1===##
print("Evaluating the uncertainties...")

##===define paths====
datapath = 'mineral_output/MC_files/'
outpath =  'mineral_output/'

#==== Constants ====#
Rearth = 6.371e3    # km

##====define a plotting function used subsequently===#
def plot_line(ax, array, colour, label=None, loc='best'):
    #---- percentiles ---#
    sixteen, median, eightyfour = np.nanpercentile(array, q=[15.9, 50, 84.1], axis=0)
    ax.fill_between(radii/Rearth, sixteen, eightyfour, color=colour, alpha=0.5, zorder=-1)
    ax.plot(radii/Rearth, median, colour, linewidth=3, zorder=1, label=label, alpha=0.7)
    ax.legend(loc=loc)
    return np.round(sixteen,8), np.round(median,8), np.round(eightyfour,8)


for sname in planet_labels:
    print(r'Run for '+sname+'...')

    samplepath = datapath + sname + '/'
    #list_of_csv = os.listdir(samplepath)
    list_of_csv = next(os.walk(samplepath))[2]  #only list files in the top directory and 1-dirpath, 2-dirnames, 3-filenames

    radii = np.linspace(0, Rearth*planet_R[planet_labels.index(sname)], 1000) #np.linspace(0, Rearth, 1000)
    Nradii = len(radii)
    Nfile = len(list_of_csv)
    plg = np.array([[np.nan]*Nradii]*Nfile) 
    O = np.array([[np.nan]*Nradii]*Nfile) 
    Wad = np.array([[np.nan]*Nradii]*Nfile)
    Ring = np.array([[np.nan]*Nradii]*Nfile)
    Opx = np.array([[np.nan]*Nradii]*Nfile)
    Cpx = np.array([[np.nan]*Nradii]*Nfile) 
    Hpcpx = np.array([[np.nan]*Nradii]*Nfile) 
    Aki = np.array([[np.nan]*Nradii]*Nfile) 
    Gt =np.array([[np.nan]*Nradii]*Nfile) 
    Pv = np.array([[np.nan]*Nradii]*Nfile) 
    capv = np.array([[np.nan]*Nradii]*Nfile) 
    Ppv = np.array([[np.nan]*Nradii]*Nfile) 
    Wus = np.array([[np.nan]*Nradii]*Nfile)
    CF = np.array([[np.nan]*Nradii]*Nfile)
    q = np.array([[np.nan]*Nradii]*Nfile) 
    coe = np.array([[np.nan]*Nradii]*Nfile) 
    Sp = np.array([[np.nan]*Nradii]*Nfile) 
    stv = np.array([[np.nan]*Nradii]*Nfile) 
    seif = np.array([[np.nan]*Nradii]*Nfile) 
    density = np.array([[np.nan]*Nradii]*Nfile) 
    pressure = np.array([[np.nan]*Nradii]*Nfile) 
    temperature = np.array([[np.nan]*Nradii]*Nfile) 
    
    Rcore = np.zeros((len(list_of_csv)))


#======== Do stuff =============#     

    for i,file in enumerate(list_of_csv):
        
        if file=='.directory':
            continue
       
        #print('file:', file)
        df = pd.read_csv(samplepath+file)
            
        R = (df['radius']/1e3)     # km
        Rcore[i] = R[501]/Rearth
            
        def eval_spline(name, array):
            if name in df.keys():
                f = interp1d(R, df[name], fill_value=float("nan"), bounds_error=False)
                cut = np.argmin(abs(radii-R[0])) + 1
                array[i, :cut] = np.nan 
                array[i, cut:] = f(radii[cut:])
                    
    
            
        eval_spline('O', O)
        eval_spline('Wad', Wad)
        eval_spline('Ring', Ring)
        eval_spline('Opx', Opx)
        eval_spline('Cpx', Cpx)
        eval_spline('C2/c', Hpcpx)
        eval_spline('Gt', Gt)
        eval_spline('Aki', Aki)
        eval_spline('st', stv)
        eval_spline('Pl', plg)
        eval_spline('Sp', Sp)
        eval_spline('P', Pv)
        eval_spline('Pp', Ppv)
        eval_spline('ca-p', capv)
        eval_spline('Cf', CF)
        eval_spline('Wus', Wus)
        eval_spline('density', density)
        eval_spline('pressure', pressure)
        eval_spline('temperature', temperature)
 

    #br    eakpoint()
    #========= Plot =======#
    print("Visualising the results...")
    
    fig, ax = plt.subplots(4,2, figsize=(16, 16))
    for i in range(2):
        for j in range(4):
            _ax = ax[j,i]

            #### core patch
            core_percentiles = np.percentile(Rcore, [15.9, 50, 84.1])
            _ax.fill_betweenx([0,12000], 0, core_percentiles[2], color='grey', zorder=-1, alpha=0.2)
            _ax.fill_betweenx([0,12000], 0, core_percentiles[0], color='darkgrey', zorder=-1, alpha=0.2)
            #_ax.fill_between(R_array, np.zeros_like(R_array), cummulative_core_func(R_array), color='lightgrey', zorder=1)
            _ax.axvline(core_percentiles[1], ls='--', color='darkgrey', linewidth=1)


            if j==3:
                _ax.set_xlabel('Radius [R$_\oplus$]')
            else:
                _ax.set_xlabel(' ')

            if j==2:
                if i==1:
                    # Density
                    _ax.set_ylabel('Density [g/cm$^3$]')
                    vertical_position = 3.6
                    vertical_position_text = 3.9
                if i==0:
                    # Perovskites
                    _ax.set_ylabel('Volume fraction [%]')
                    vertical_position = 90
                    vertical_position_text = 93
            elif j==3:
                if i==0:
                    _ax.set_ylabel('Pressure [GPa]')
                    vertical_position = 50
                    vertical_position_text = 65
                if i==1:
                    _ax.set_ylabel('Temperature [K]')
                    vertical_position = 2000
                    vertical_position_text = 2200
            else:
                _ax.set_ylabel('Volume fraction [%]')
                vertical_position = 85
                vertical_position_text = 90

            _ax.errorbar(core_percentiles[1], vertical_position, xerr=[[core_percentiles[1]-core_percentiles[0]],[core_percentiles[2]-core_percentiles[1]]], color='k',capsize=5, elinewidth=0.5, marker='o')
            _ax.text(core_percentiles[1], vertical_position_text, '$R_{core}$')

            _ax.set_ylim([0, 100])
            _ax.set_xlim([core_percentiles[0]-0.05, planet_R[planet_labels.index(sname)]])

    ax[0,0].set_title('')
    O_p16, O_p50, O_p84 = plot_line(ax[0,0], O, 'tab:blue', label='Ol', loc='upper center')
    Wad_p16, Wad_p50, Wad_p84 = plot_line(ax[0,0], Wad, 'tab:orange', label='Wad', loc='upper center')
    Ring_p16, Ring_p50, Ring_p84 = plot_line(ax[0,0], Ring, 'tab:green', label='Ring', loc='upper center')
    #ax[0,0].plot(df_planet.radius/1e3, df_planet.O, 'darkkhaki', ls='dashed')
    
    ax[0,1].set_title('')
    
    Opx_p16, Opx_p50, Opx_p84 = plot_line(ax[0,1], Opx, 'tab:blue', label='Opx', loc='upper center')
    Cpx_p16, Cpx_p50, Cpx_p84 = plot_line(ax[0,1], Cpx, 'tab:orange', label='Cpx', loc='upper center')
    Hpcpx_p16, Hpcpx_p50, Hpcpx_p84 = plot_line(ax[0,1], Hpcpx, 'tab:green', label='hp-cpx', loc='upper center')
    
    ax[1,0].set_title('') 
    Gt_p16, Gt_p50, Gt_p84 = plot_line(ax[1,0], Gt, 'tab:blue', label='Gt')
    stv_p16, stv_p50, stv_p84 = plot_line(ax[1,0], stv, 'tab:orange', 'Stv')
    if np.isnan(plg).all() == False: plg_p16, plg_p50, plg_p84 = plot_line(ax[1,0], plg, 'tab:green', label='Plg')
    if np.isnan(Sp).all() == False: Sp_p16, Sp_p50, Sp_p84 = plot_line(ax[1,0], Sp, 'tab:purple', label='Sp')

    ax[1,1].set_title('')
    Wus_p16, Wus_p50, Wus_p84 = plot_line(ax[1,1], Wus, 'tab:blue', label='Wus')
    CF_p16, CF_p50, CF_p84 = plot_line(ax[1,1], CF, 'tab:orange', label='CF')

    ax[2,0].set_title('')
    Pv_p16, Pv_p50, Pv_p84 = plot_line(ax[2,0], Pv, 'tab:blue', label='Mg-pv (Bridgm)') #, loc='center')
    Ppv_p16, Ppv_p50, Ppv_p84 = plot_line(ax[2,0], Ppv, 'tab:orange', label='Mg-postpv') #, loc='center')
    capv_p16, capv_p50, capv_p84 = plot_line(ax[2,0], capv, 'tab:green', label='Ca-pv') #, loc='center')


    ax[2,1].set_title('')
    den_p16, den_p50, den_p84 = plot_line(ax[2,1], density/1e3, 'k', label='Density') 
    ax[2,1].set_ylim([3, 13]) 

    ax[3,0].set_ylabel('')
    P_p16, P_p50, P_p84 = plot_line(ax[3,0], pressure/1e9, 'b', label='Pressure')
    ax[3,0].set_ylim([0, 300])
   
    
    ax[3,1].set_ylabel('')
    T_p16, T_p50, T_p84 = plot_line(ax[3,1], temperature, 'r', label='Temperature')
    ax[3,1].set_ylim([1500, 5500])

    
    ##------output the MC evaluation results----------
    MC_eval = pd.DataFrame({'R/Rearth': np.flipud(radii/Rearth),
                            'Ol_p50': np.flipud(O_p50), 'Ol_p16': np.flipud(O_p16), 'Ol_p84': np.flipud(O_p84),
                            'Opx_p50': np.flipud(Opx_p50), 'Opx_p16': np.flipud(Opx_p16), 'Opx_p84': np.flipud(Opx_p84),
                            'Cpx_p50': np.flipud(Cpx_p50), 'Cpx_p16': np.flipud(Cpx_p16), 'Cpx_p84': np.flipud(Cpx_p84),
                            'Hpcpx_p50': np.flipud(Hpcpx_p50), 'Hpcpx_p16': np.flipud(Hpcpx_p16), 'Hpcpx_p84': np.flipud(Hpcpx_p84),
                            'Wad_p50': np.flipud(Wad_p50), 'Wad_p16': np.flipud(Wad_p16), 'Wad_p84': np.flipud(Wad_p84),
                            'Ring_p50': np.flipud(Ring_p50), 'Ring_p16': np.flipud(Ring_p16), 'Ring_p84': np.flipud(Ring_p84),
                            'Gt_p50': np.flipud(Gt_p50), 'Gt_p16': np.flipud(Gt_p16), 'Gt_p84': np.flipud(Gt_p84),
                            'Pv_p50': np.flipud(Pv_p50), 'Pv_p16': np.flipud(Pv_p16), 'Pv_p84': np.flipud(Pv_p84),
                            'Ppv_p50': np.flipud(Ppv_p50), 'Ppv_p16': np.flipud(Ppv_p16), 'Ppv_p84': np.flipud(Ppv_p84),
                            'Wus_p50': np.flipud(Wus_p50), 'Wus_p16': np.flipud(Wus_p16), 'Wus_p84': np.flipud(Wus_p84),
                            'den_p50': np.flipud(den_p50), 'den_p16': np.flipud(den_p16), 'den_p84': np.flipud(den_p84),
                            'P_p50': np.flipud(P_p50), 'P_p16': np.flipud(P_p16), 'P_p84': np.flipud(P_p84),
                            'T_p50': np.flipud(T_p50), 'T_p16': np.flipud(T_p16), 'T_p84': np.flipud(T_p84),
                            })

    MC_eval.to_csv(outpath+'mineralogy_'+sname+'_MC_eval.csv')

    ##-----generate the plot--------------------
    plt.tight_layout()
    plt.savefig(outpath+'mineralogy_'+sname+'_MC_eval.pdf') 

print("Evaluation done.")

