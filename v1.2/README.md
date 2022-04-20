# ExoInt_v1.2
This "v1.2" version is a moderate update to “ExoInt” that is created and maintained by Haiyang S. Wang.
The updates include:
1. Allow silicon to be partially in the core, depending on the planet's oxidation state
2. Final output is constrained by a P-value that reflects the modeling degeneracy of the fractionation of Fe and Ni (and light elements, e.g. Si and S) between mantle and core.
3. The possibility to adjust the devolatilisation levels by 3sigma or 5sigma (or their upper/lower limits)

The codes are developed in IDL, and a Python version (corresponds to this v1.2 version) has also been available here https://github.com/astro-seanwhy/ExoInt/tree/master/pyExoInt 

Please run the codes starting with ‘exoconstraints.pro’.
An example has been given on how to implement the code.

Please hands on by testing with any other stellar abundances (in dex, relative to the Sun) of your interest (But the predictions work best for potential silicate planets around Sun-like stars).

In the folder "outputs", you can find your modelling results of the mantle and core composition as well as core mass fraction, which approximate the most plausible properties of a habitable-zone rocky planet around your sample of star.


The primary citation should be given to:
Wang, H. S., Liu, F., Ireland, T., Brasser, R., Yong, D., and Lineweaver, C. H. 2019. Enhanced constraints on the interior composition and structure of terrestrial exoplanets. MNRAS 482:2222-2233. doi.org/10.1093/mnras/sty2749

Relevant references for this version include: 
Wang, H. S., Quanz, S. P., Yong, D., Liu, F., Seidler, F., Acuna, L., Mojzsis, S. J. Detailed chemical compositions of planet hosting stars: II. Exploration of the interiors of terrestrial-type exoplanets. MNRAS, in press.

Wang, H. S., Lineweaver, C. H., Quanz, S. P., Mojzsis, S. J., Ireland, T. R., Sossi, P. A., Seidler, F., and Morel T. 2022. A model Earth-sized planet in the habitable zone of α Centauri A/B. ApJ 927:134. https://doi.org/10.3847/1538-4357/ac4e8c

Copyright@Haiyang S. Wang, with the following permissions:

The code "exoconstraints.pro" is free to download, use and modify, should it be for academic and non-profit purpose.

The codes "exomodel.pro" and "chemsysmodel.pro" are free to download and use, but not allowed to modify (except for choosing different keywords if necessary) without the permission of Haiyang S. Wang.

The use of any part of the codes should have the citation given to the papers mentioned above.

Questions, comments, and suggestions may be raised in the 'Issues' tab, or otherwise directly sent to haiwang@phys.ethz.ch (for requests, in particular).
