from exomodel import exomodel

####----render stellar abundances------------
## Please format your stellar abundances to be differential, i.e., [X/H], as the example files (including the header)
## examplestar1 has symmetric error bars while examplestar2 has asymmetric error bars 
## and place them in the folder 'data/stellar_abu/' (by default)

datapath = 'data/stellar_abu/'
starnames=['examplestar1', 'examplestar2']
filetype='.txt'  ##specify your data file extension

for star in starnames:
    print(r'Run for '+star+'...')

    result_df = exomodel(datapath+'abu_'+star+filetype,
                     asymerr = False , ## Asymetric error bars on your stellar abundances (in [X/H])? No - False (default); Yes - True
                     Nsim = 1e3,  ## the number of Monte-Carlo simulations for drawing stellar abundances from your input file
                     refsolar = 'A21',   ## Please choose one of the following reference solar abunances (A(X), in dex)
                                         ## 'A09'- Asplund et al. 2009 ARAA 
                                         ## 'l09' - Lodders et al. 2009 Book: Landolt- Börnstein, New Series, Vol. VI/4B
                                         ## 'W19' - Wang et al. 2019 ICARUS
                                         ## 'A21' - Asplund et al. 2021 A&A solar abundance will be adopted."
                                         ## 'custom' -- if opt for 'custom', please make sure you have formatted your reference solar abundances as one of the given solar abundances above and have placed the file under the 'data' folder as well; the error bars on your solar abundances (A(X), in dex) should be symetric
                     verbose=True  ## automatically print out the simulation results to the screen (True, by default) or not (False) 
    )
