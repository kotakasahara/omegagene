import numpy as np

class Const(object):
    DEBUG = True
    elem_charge = 1.60217657e-19 # [As]
    avogadro = 6.0221413e23 # [mol^-1]
    permittivity = 8.85418782e-12 # [m^-3 kg^-1 s^4 A^2]
    #charge_coeff = elem_charge**2 * avogadro / ( 4 * np.pi * permittivity * 1e-10 * 4.184 * 1e3) 
    ## [Angstrome cal mol-1]
    charge_coeff =  332.06378
    eps = 1e-8
    force_vel = 4.184e-4
    gas_const = 8.31451

    joule_cal = 4.184
    kinetic_coeff = 1e7 /(joule_cal*1e3)*0.5

