This "v1.2" version is a moderate update to “ExoInt” that is created and maintained by Haiyang S. Wang.
The updates include:
1. Allow silicon to be partially in the core, depending on the planet's oxidation state
2. Final output is constrained by a P-value that reflects the modeling degeneracy of the fractionation of Fe and Ni (and light elements, e.g.Si and S) between mantle and core.
3. The possibility to adjust the devolatilisation levels by 3sigma or 5sigma (or their upper/lower limits)

The citation should still be given to:
Wang, H. S., Liu, F., Ireland, T., Brasser, R., Yong, D., and Lineweaver, C. H. 2018. Enhanced constraints on the interior composition and structure of terrestrial exoplanets. MNRAS, in press. doi.org/10.1093/mnras/sty2749

The codes are written in IDL, and a Python version has been developed and will be online as well, soon.

Please run the codes starting with ‘exoconstraints.pro’.
An example has been given how to implement the code.

Please hands on by testing with any other stellar abundances (in dex, relative to the Sun) of your interest (But the predictions work best for potential silicate planets around Sun-like stars).

In the folder "outputs", you can find your modelling results of the mantle and core composition as well as core mass fraction, which approximate the most plausible properties of a habitable-zone rocky planet around your sample of star.

Copyright@Haiyang S. Wang, with the following permissions:

The code "exoconstraints.pro" is free to download, use and modify, should it be for academic and non-profit purpose.

The codes "exomodel.pro" and "chemsysmodel.pro" are free to download and use, but not allowed to modify (except for choosing different keywords if necessary) without the permission of Haiyang S. Wang.

The use of any part of the codes should have the citation given to the paper mentioned above.

Questions, comments, and requests for any particular permission should be addressed to haiwang@ethz.ch.
