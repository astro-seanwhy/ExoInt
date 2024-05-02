# ExoInt
Rocky planets as devolatilized stars (Source codes for Wang et al. 2019, MNRAS, doi.org/10.1093/mnras/sty2749). Purpose: devolatilize stellar abundances to produce rocky exoplanetary bulk composition, which is used togegther with mass/radius information to model the exoplanetary interior structure and mineralogy. 

The codes, named “ExoInt”, are created and maintained by Haiyang S. Wang, Copyright @ 2019-2024. 

---Updates---

A moderate updated version (v1.2) in IDL is available under the folder "v1.2". 

A Python version corresponding to the v1.2 (IDL) version is available under the folder "pyExoInt".

The 'Mineral' moduel has been added to output the mineralogy based on the stoichemitric output of mantle and core compositions, core mass fraction, along with the given mass and radius information.

---How to run ExoInt ---

To devolatilize your sample of stellar abundances, please prepare them as differential abundances (i.e. [X/H], dex) according to the format of the example given and then place your datafiles under the folder “data/stellar_abu”. To run the codes, simply type “idl exoconstraints.pro” (or run it with an IDL interface), and the modelling results for both the planetary bulk composition (devolatilized stellar composition) and interior compositions (mantle, core, and core mass fraction) will be automatically stored under the folder “output”. 

---Licence ---

The free programme is distributed under the GNU GPLv3 open-acess license (unless otherwise specified for some modules), in the hope that it will be useful.

You are particularly welcome to make pull requests and contribute to the further development of the programme. 

---Citations---

The use of the programme should have the primary citation given to:

Wang, H. S., Liu, F., Ireland, T., Brasser, R., Yong, D., and Lineweaver, C. H. 2019. Enhanced constraints on the interior composition and structure of terrestrial exoplanets. MNRAS 482:2222-2233. https://doi.org/10.1093/mnras/sty2749 

For using the 'Mineral' module, please refer to the IDL v1.2 version or the pyExoInt for more appropriate citations.

---Contact---

Questions, comments, and suggestions may be raised in the 'Issues' tab, or otherwise directly sent to haiyang.wang@sund.ku.dk or haiwang@ethz.ch (for requests, in particular). 

