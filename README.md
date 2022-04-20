# ExoInt
Rocky planets as devolatilized stars (Source codes for Wang et al. 2019, MNRAS, doi.org/10.1093/mnras/sty2749). Purpose: devolatilize stellar abundances to produce rocky exoplanetary bulk composition, with which constrain the modelling of the exoplanet interiors. 

The codes, named “ExoInt”, are created and maintained by Haiyang S. Wang. 

Primary citation: Wang, H. S., Liu, F., Ireland, T., Brasser, R., Yong, D., and Lineweaver, C. H. 2019. Enhanced constraints on the interior composition and structure of terrestrial exoplanets. MNRAS 482:2222-2233. https://doi.org/10.1093/mnras/sty2749

A moderate updated version (v1.2) in IDL is available under the folder "v1.2". 

A Python version corresponding to the v1.2 (IDL) version is available under the folder "pyExoInt".

To devolatilize your sample of stellar abundances, please prepare them as differential abundances (i.e. [X/H], dex) according to the format of the example given and then place your datafiles under the folder “data”. To run the codes, simply type “idl exoconstraints.pro” (or run it with an IDL interface), and the modelling results for both the planetary bulk composition (devolatilized stellar composition) and interior compositions (mantle, core, and core mass fraction) will be automatically stored under the folder “outputs”. 

Copyright@Haiyang S. Wang, with the following permissions: 

The code "exoconstraints.pro" is free to download, use and modify, should it be for academic and non-profit purpose. 

The codes "exomodel.pro" and "chemsysmodel.pro" are free to download and use, but not allowed to modify without the permission of Haiyang S. Wang. Please be compliant with the GNU GPLv3 open-acess license if you would make any modification to these codes. And this license applies to all versions of ExoInt. 

The use of the codes or part of the codes should have the citation given to the paper mentioned above. 

Questions, comments, and suggestions may be raised in the 'Issues' tab, or otherwise directly sent to haiwang@phys.ethz.ch (for requests, in particular). 

