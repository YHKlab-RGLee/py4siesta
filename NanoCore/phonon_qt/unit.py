#
# UNIT conversion factors
#

from math import pi

fu2THz = 2*pi*15.6330214     # THz
fu2eV  = fu2THz * 4.136*10**-3  # eV
fu2meV = fu2THz * 4.136         # meV
fu2cm1 = fu2THz * 33.356
vasp2cm1 = 521.4708336735473    # ev/ang -> cm^-1
factor = vasp2cm1

hbar = 4.1356675*10**-15 / (2*pi) # eV.s
k_boltz = 1.38062*10**-23 #J/K
k_boltz = k_boltz/(1.602*10**-19) #eV/K 

conv_factor_to_thz = {
'VASP'      : 15.633302,
'WIEN2k'    : 3.44595837,
'QE'        : 108.97077,
'ABINIT'    : 21.49068,
'SIESTA'    : 21.49068,
'Elk'       : 154.10794,
'CRYSTAL'   : 15.633302,
'DFTB+'     : 154.10794,
'TURBOMOLE' : 154.10794,
}
