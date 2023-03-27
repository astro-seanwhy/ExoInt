"""
This program calls the three related mineralogy codes as an one go or otherwise run individually by giving a flag

Copyright (C) 2021-2023  Haiyang S. Wang

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

import os
import sys
import warnings
"""
---Instructions---:
## 1. please provide the input MRfile from the command line, as sys.argv[1]; by default, the MRfile is stored in the "MRfiles" folder, in which an example file is given.
## 2. you can choose to run the three individual codes individually, by giving a flag "-quick", "-MCstep1", or "-MCstep2". These flags should be used individually; namely, it's not designed for a combined use (which is also unnecessary).
      [note: without any flag, the ...OneGo.py code will run the three codes of the 'Mineral' module as a streamline automatically, i.e. '...quick_solver.py' -> '...MC_solver_step1.py' -> '...MC_solver_step2.py'.]
## 3. you can also add two optional args -- (approximate) number of files to be generated (num_files = 100 by default) and the number of cores to be used (num_cores = 25 by default).
      E.g. "python3 mineral_MC_solver_step2.py MRfiles/MR_example.csv 200 30" ; here you indicate (approximate) 200 files to be generated, which will be running on 30 cores. 
      [pleaste note that:
      i. these two arguments are recommended to be used jointly and in the order of num_files (first) followed by num_cores. In case of only one number given, it is interpreted as the number of files to be generated, with the default number (25) of cores.
      ii. these two arguments can be used together with one of the flags above, and the flag can be placed in any order beyond the input MRfile.]
"""

if len(sys.argv) >2:
    for arg in sys.argv[2:]:
        if arg.isdigit() == False and (arg != '-quick' and arg != '-MCstep1' and arg != '-MCstep2'):
            #warnings.warn(arg + " is not a permitted flag; the flag is ignored and the code is continuing")
            print(arg)
            raise Exception( arg + "is not a permitted flag")

if '-quick' in sys.argv:
    ##---Run for the mean mineralogy, as a quick solution ---
    print('>>> Run mineral_quick_solver.py...')
    os.system("python3 mineral_quick_solver.py " + sys.argv[1])

elif '-MCstep1' in sys.argv:
    ##---Run for the Monte Carlo simulations to asess the mineralogy uncertainty ---
    print('>>> Run mineral_MC_solver_step1.py ...') # for ' + str(Nrun) + " round(s), with " + str(Nrun*39) + " MC output for evaluation...")

    if len(sys.argv) > 2:
        optargs = sys.argv[2:]
        #for arg in optargs:
        #    if type(arg) == str and arg!= '-MCstep1':
        #        raise Exception(arg+ " is not a permitted flag!")
        if "-MCstep1" in optargs:
            optargs.remove("-MCstep1")
        if len(optargs) == 0:
            os.system("python3 mineral_MC_solver_step1.py " + sys.argv[1])
        if len(optargs) == 1 and optargs[0].isdigit():
            print("[ACTION REQUIRED]: you have only provide one number, which will be interpreted as the (approximate) number of files to be generated. The number of cores to use is 25, by default.")
            print ("if that's what you wanted, type 'coninue' to continue the sim")
            print ("otherwise, please provide two numbers (for 'files', followed by 'cores') to be clearer")
            breakpoint()
            os.system("python3 mineral_MC_solver_step1.py " + sys.argv[1] + " " + optargs[0])
        if len(optargs) ==2 and optargs[0].isdigit() and optargs[1].isdigit():
            os.system("python3 mineral_MC_solver_step1.py " + sys.argv[1] + " " + optargs[0] + ' ' + optargs[1]) # + " " + arg_option)
        if len(optargs) >= 3:
            print("here!")
            raise Exception("Sorry, you provided more numbers than the required (2 max.)!")
    else:
        os.system("python3 mineral_MC_solver_step1.py " + sys.argv[1])


elif '-MCstep2' in sys.argv:
    ##---Run for generating the MC uncertainty plot ---
    print('>>> Run mineral_MC_solver_step2.py...')
    os.system("python3 mineral_MC_solver_step2.py " + sys.argv[1])

else:
    ##--Run for all the related mineralogy codes automatically as a sequence ---
    print('>>> Run mineral_quick_solver.py...')
    os.system("python3 mineral_quick_solver.py " + sys.argv[1])

   
    print('>>> Run mineral_MC_solver_step1.py ...')
    if len(sys.argv) > 2 and (sys.argv[2].isdigit() and sys.argv[3].isdigit()):
        os.system("python3 mineral_MC_solver_step1.py " + sys.argv[1] + " " + sys.argv[2] + ' ' + sys.argv[3]) # + " " + arg_option)
    else: 
        os.system("python3 mineral_MC_solver_step1.py " + sys.argv[1])



    print('>>> Run mineral_MC_solver_step2.py...')
    os.system("python3 mineral_MC_solver_step2.py " + sys.argv[1])


