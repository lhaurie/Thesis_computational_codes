# Thesis_computational_codes

For start, please run Minsquares_2x2_ with modules_2x2_ as it describes the simpler case of a paramagnetic metal.
xx contains the three self-consistent parameters: mu (chemical potential), and the parameter e and p (cf phd manuscript for theoretical details : https://theses.fr/2024UPASP073).
modules contains the computation of the correlation function in momentum space and their average. The file is called by Minsquares for a broydn minimization.
The parameters to modify are the following:
electron density n, temperature temp, Coulomb interaction U, tight-binding parameter t, initial guesses x01, x02, x03
These codes must be compiled using llapack library.
After convergence, the code can generate a text file to store the datas, please uncomment in minsquares the saving paragraph if required.


Minsquares_superconductivity_longer_ranged and modules_superconductivity_longer_ranged are an extension of the composite operators method to superconductivity and
longer ranged hoppings. Note that putting theta (the superconducting order parameter) to zero and longer ranged hoppings to zero gives back the same results as the 
paramagnetic metallic case.

Min_SO and mod_SO contains the study of orbital selective mott phase (OSMP) and orbital uniform phase (OUP). 
The python code display_bands can plot for a given set of converged parameters the bands and Fermi surface of the two orbital Hubbard model.
