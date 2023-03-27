"""
This program presents a low-level "wrapper" to Perple_X (Connolly 2005,2009), which
computes the mineral phase assemblage given composition, pressure and temperature.
Originally adapted from the earlier - and now deprecated - version of ExoPlex 
(https://github.com/amloren1/ExoPlex), which is superseded by https://github.com/CaymanUnterborn/ExoPlex.
ExoPlex is released under the GNU General Public License version 3. 
This version is further modified by Fabian L. Seidler and Haiyang S. Wang

Copyright (C) 2018-2023 Alejandro Lorenzo, Cayman Unterborn, Fabian L. Seidler, Haiyang S. Wang

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as
published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import os
import sys
import pexpect as pe
import numpy as np
import time
from subprocess import Popen, PIPE, STDOUT

PerPlex_path = os.path.dirname(os.path.realpath(__file__))+"/PerPlex"
 
    
    
def calculate_mineralogy(mantle_chemsys, path_name, verbose=False):
    """
    Computes the mineralogy using Perple_X (Connolly 2005, 2009).
    This means calling the Perple_X executables in the following order:

    build : Creates the parameter files for the simualtion; it contains the
            bulk chemistry, the thermodynamic data for the mineral phases of
            interest and the path to the expected pressure-temperature gradient
            of the mantle (unfortunately, we cannot compute this self-consistently;
            however, Perple_X could allow for a parametrized gradient, but we choose
            a pre-defined gradient instead).
    vertex : Constructs the pseudocompounds and performs the free energy minimization
    werami : Computes thermodynamic properties of interest (e.g., density, mineralogy)

    Parameters
    ----------
        mantle_chemsys : dict
            Dictionary containing mantle chemistry
        path_name : str
            where to store the temporary files and results of Perple_X
    """
    
    abundances = str(mantle_chemsys.get('SiO2')) + ' ' + str(mantle_chemsys.get('MgO')) + ' ' \
             + str(mantle_chemsys.get('FeO')) + ' ' + str(mantle_chemsys.get('Al2O3')) \
             + ' ' + str(mantle_chemsys.get('CaO'))+ ' ' + str(mantle_chemsys.get('Na2O'))

    #print(abundances)

    projectname = 'mantlephys'
    filename = path_name + '/' + projectname
    
    assert (len(filename)<91), "Path to Perple_X must not have more than 90 characters!"  
    
    
    #### Build ####
    p = Popen(PerPlex_path+'/./build', stdout=PIPE, stdin=PIPE, stderr=STDOUT, encoding='utf8')
    
    stdin= str(filename)+'\n'\
    +PerPlex_path+'/stx11ver.dat\n' \
    +PerPlex_path+'/perplex_option.dat\n' \
    'N\n3\nN\nN\n' \
    'SIO2\nMGO\nFEO\nAL2O3\nCAO\nNA2O\n\n'\
    'Y\n'\
    +path_name+'/mantlegrid.dat\n'\
    '1\nY\n'\
    +abundances+'\n' \
    'N\nN\nY\n'\
    +PerPlex_path+'/stx11_solution_model.dat\n'\
    'C2/c\nWus\nPv\nPl\nSp\nO\nWad\nRing\nOpx\nCpx\nAki\nGt\nPpv\nCF\n\n'\
    +filename+'calc\n'
   
    #print(stdin)
    
    stdout = p.communicate(input=stdin, timeout=60)[0]
    
    #logfile = open(path_name+'/build.log', 'w')
    #for line in stdout:
    #    logfile.write(line)
    p.wait()
    
    if verbose: print("Done with Build, moving on to Vertex")

    
    #### Vertex ####
    p = Popen(PerPlex_path+'/./vertex', stdout=PIPE, stdin=PIPE, stderr=STDOUT, encoding='utf8')
    ####
    
    stdin= str(filename)+'\n0\n'
    
    stdout = p.communicate(input=stdin, timeout=90)[0]
    
    logfile = open(path_name+'/vertex.log', 'w')
    for line in stdout:
        logfile.write(line)
    p.wait()

    
    if verbose: print('Finished with Vertex, beginning Werami')

    
    #### Werami ####
    p = pe.spawn(PerPlex_path+"/./werami",timeout=2400)


    p.sendline(filename)
    # select 1D path
    p.sendline('3')
    # Below, select parameters density, Grüneisen parameter, adiabatic bulk modulus.
    # Ns for no calculating individual phase properties
    p.sendline('2') # density
    p.sendline('N')
    p.sendline('9') # Grüneisen parameter
    p.sendline('N')
    p.sendline('10')    # adabatic bulk modulus
    p.sendline('N')
    p.sendline('13')
    p.sendline('N')
    p.sendline('14')
    p.sendline('N')
    
    ##### the next lines will pass requests to perplex to print phases and their proportions into the .tab file

    phases = ['Wus', 'Pv', 'Ppv', 'CF', 'O', 'Wad', 'Ring', 'Pl', 'C2/c', 'Opx', 'Cpx', 'Aki', 'Gt', 'st', 'q', 'ca-pv', 'coe', 'ky', 'Sp', 'neph', 'seif' ]
    
    p.sendline('7')
    for phase in phases:
        p.sendline(phase)
        ii = p.expect(['try again:','proportions keyword'])
        if ii == 0:
            continue
            #sys.exit()
        elif phase != 'seif' and ii == 1:
            p.sendline('7')
            continue
        else:
            continue
        
    #exit parameter choosing

    p.sendline('0')
    p.sendline('0')
    p.expect('0 - EXIT')
    #p.sendline('0')
    #p.logfile = open(path_name+'/werami.log','wb')
    p.terminate()
    
    if verbose: print("Done with PerPlex")
    
    


def make_mantle_grid(path_name):
    

    projectname = 'mantlephys'
    filename = path_name + '/' + projectname

    file = open(filename+'_1.tab')

    tmp_file = file.readlines()
    num_rows = len(tmp_file[9:])
    num_columns = len(tmp_file[8].split())

    header = tmp_file[8].strip('\n').split()

    data = tmp_file[9:]
    grid = np.zeros((num_rows,num_columns))

    for i in range(num_rows):
        #for j in range(num_columns):
        columns = data[i].strip('\n').split()
        grid[i] = [float(j) for j in columns]
        

    temperature_grid = [row[2] for row in grid]
    pressure_grid = [row[1] for row in grid]
    density_grid = [row[3] for row in grid]
    gr_grid = [row[4] for row in grid]
    KS_grid = [row[5] for row in grid]

    keys = ['temperature','pressure','density','Grueneisen','KS']
    mantle_dict = dict(list(zip(keys,[temperature_grid,pressure_grid,density_grid,gr_grid,KS_grid])))

    return mantle_dict


def make_mineralogy_grid(path_name):
    """
    Reads Perple_X output files and retrieves important structural 
    properties, i.e. density, Grüneisen-parameter (gr) and isentropic
    bulk modulus (KS), as well as the mineralogy (vol%)

    Parameters
    ----------
        path_name : str
            path to Perple_X output files
    Returns
    -------
        mantle_dict : dict
            dictionary with mantle properties as numpy arrays; 
            keys = ['temperature','pressure','density','Grueneisen','KS', 'speeds', 'phases']
            units:
                temperature : K
                pressure : bar
                density : kg/m^3
                Grueneisen: dimensionless
                KS : Pa
                speeds : m/s
                phases : vol%
        Phases : list
            names of the mineral phases

    """

    projectname = 'mantlephys'
    filename = path_name + '/' + projectname

    file = open(filename+'_1.tab')

    tmp_file = file.readlines()
    num_rows = len(tmp_file[9:])
    num_columns = len(tmp_file[8].split())

    header = tmp_file[8].strip('\n').split()

    data = tmp_file[9:]
    grid = np.zeros((num_rows,num_columns))

    for i in range(num_rows):
        #for j in range(num_columns):
        columns = data[i].strip('\n').split()
        grid[i] = [float(j) for j in columns]
        

    temperature_grid = [row[2] for row in grid]
    pressure_grid = [row[1] for row in grid]
    density_grid = [row[3] for row in grid]
    gr_grid = [row[4] for row in grid]
    KS_grid = [row[5] for row in grid]

    keys = ['temperature','pressure','density','Grueneisen','KS']
    mantle_dict = dict(list(zip(keys,[temperature_grid,pressure_grid,density_grid,gr_grid,KS_grid])))
    
    speed_grid = [[row[6],row[7]] for row in grid]
    
    mantle_dict['speeds'] = speed_grid

    Phases = header[8:]

    for i in range(len(Phases)):
        Phases[i] = Phases[i].strip(",vo%")

    num_phases = len(grid[0][8:])
    phases_grid = np.zeros((num_rows,num_phases))
    for i in range(num_rows):
        phases_grid[i] = grid[i][8:]
    phase_grid = [row[8:] for row in grid]
    
    mantle_dict['phases'] = phase_grid

    return mantle_dict,Phases
    
