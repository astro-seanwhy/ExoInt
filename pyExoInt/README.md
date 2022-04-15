# pyExoInt
pyExoInt is a python version of ExoInt_v1.2 (IDL) -- A devolatilization and interior modeling package for rocky planets, which is created and maintained by Haiyang S. Wang (HSW) and is copyrighted with the GNU GPLv3 open-acess license.

Disclaimer: The python version is contributed by Fabian Seidler, under the guidance of HSW. It is not guaranteed to be identical to the IDL version, but it contains all the main features of the original code and may be useful to users who are more familiar with the Python language.

------------------- Running ExoInt --------------------

To run pyExoInt, start with 'exoconstraints.py', which calls the function 'exomodel' from 'exomodel.py'. "chemsysmodel" is further called by "exomodel".

The minimum working example on how to do this is illustrated in 'exoconstraints.py', where you can adjust these Keyword arguments:
    -Nsim      (optional) number of random draws from the stellar abundance distribution; default is 1e4; best results with high number of runs, but takes longer
    -refsolar  (optional) Solar reference abundances. Options are:
                    'A21'     Asplund et al. 2021, A&A, default
                    'A09'     Asplund et al. 2009, ARAA
                    'L09'     Lodders et al. 2009, in the Solar System (Trumper ed.) 10.1007/978-3-540-88055-4_34
                    'W19'     Wang et al. 2019, Icarus

The final results are stored in the "output" folder, and are named by the filename of the input stellar abundances, followed with suffixes '_results.csv' and '_results_formated.txt' -- i.e. two different files with identical results.

------------------- IMPORTANT ----------------------

Please note that the input stellar abundances have to be 'differential' abundances ([X/H]) in dex, i.e. relative to a reference solar abundance. It is crucial to format them in the same way as in the given example file of stellar abundances in the "stellar_abu_files" folder, i.e. have the 'abu' (the abundance) and 'err' (the uncertainty) column present.


------------------- Citations ----------------------

The primary citation should be given to:
Wang, H. S., Liu, F., Ireland, T., Brasser, R., Yong, D., and Lineweaver, C. H. 2019. Enhanced constraints on the interior composition and structure of terrestrial exoplanets. MNRAS 482:2222-2233. doi.org/10.1093/mnras/sty2749

Relevant reference for this Python version:
Wang, H. S., Quanz, S. P., Yong, D., Liu, F., Seidler, F., Acuna, L., Mojzsis, S. J. Detailed chemical compositions of planet hosting stars: II. Exploration of the interiors of terrestrial-type exoplanets. MNRAS, in press.

Contact:
Questions, comments, and suggestions may be raised in the 'Issues' tab, or otherwise directly sent to haiwang@phys.ethz.ch (for requests, in particular).
