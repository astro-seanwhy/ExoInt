pyExoInt is a python version of ExoInt_v1.2 (IDL) -- A devolatilization and interior modeling package for rocky planets, which is created and maintained by Haiyang S. Wang (HSW) and is copyrighted with the GNU GPLv3 open-acess license (unless specified otherwise for some modules). 

Disclaimer: 
"""
The python version is not be guaranteed to be identical to the IDL version of the original ExoInt_v1.2, but it contains all the main features of the original code and may be useful to users who are more familiar with the Python language.

It was co-developed by Haiyang S. Wang and Fabian L. Seidler @ 2021-2022.
It has been further modified and maintained by Haiyang S. Wang @ 2022-2023.
"""

------------------- Running ExoInt --------------------
To run pyExoInt, start with 'exoconstraints.py', which calls the function 'exomodel' from 'exomodel.py'. "chemsysmodel" is further called by "exomodel".

----render stellar abundances------------
--Please format your stellar abundances to be differential, i.e., [X/H], as the example files (including the header)
--examplestar1 has symmetric error bars while examplestar2 has asymmetric error bars 
--and place them in the folder 'data/stellar_abu/' (by default)

----output ------
The final results are stored in the "output" folder, and are named by the filename of the input stellar abundances, followed with suffixes '_results.csv' and '_results_formated.txt' -- i.e. two different files with identical results.
In addition, three seperated files -- ...mantlecomp_final.txt, corecomp_final.txt, and fcoremass_final.txt -- are also generated, which may be called for the subsequent mineral modellling through the code "Mineral" (read further instructions there) 

------------------- Citations ----------------------
The primary citation should be given to:
Wang, H. S., Liu, F., Ireland, T., Brasser, R., Yong, D., and Lineweaver, C. H. 2019. Enhanced constraints on the interior composition and structure of terrestrial exoplanets. MNRAS 482:2222-2233. doi.org/10.1093/mnras/sty2749

Relevant references for this Python version: 
Wang, H. S., Quanz, S. P., Yong, D., Liu, F., Seidler, F., Acuna, L., Mojzsis, S. J. 2022. Detailed chemical compositions of planet hosting stars: II. Exploration of the interiors of terrestrial-type exoplanets. MNRAS 513:5829-5846. https://doi.org/10.1093/mnras/stac1119
Wang, H. S., Linweaver, C. H., Quanz, S. P., Mojzsis, S. J., Ireland, T. R., Sossi, P. A., Seidler, F., and Morel, T. 2022. A model Earth-sized planet in the habitable zone of α Centauri A/B. ApJ 927:134. https://doi.org/10.3847/1538-4357/ac4e8c


Contact:
Questions, comments, and suggestions may be raised in the 'Issues' tab, or otherwise directly sent to haiwang@phys.ethz.ch (for requests, in particular).
