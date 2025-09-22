#
# Thermoelectic parts
#

import numpy as np 
import scipy.interpolate as intp
import scipy.integrate as intg
from math import pi

factor = (1.60219 * 10**-19) * (10**12) # eV -> J -> nJ


def read_etrns(file_name):
    # read transiesta AVTRANS file
    # unit eV, 1
    E = []; T = []
    f = open(file_name)
    lines = f.readlines()
    for line in lines:
        if not line: break
        if not line.split()[0] == '#':
            E.append(float("%8.6f"%float(line.split()[0])))
            T.append(float("%8.6f"%float(line.split()[1])))
    return np.array(E), np.array(T)


def read_ptrns(trns_file):
    # read own phonon transmission calculation result
    # unit cm-1, 1
    f = open(trns_file)
    lines = f.readlines()
    T = []; omega = []
    for line in lines:
        t, omg = line.split()
        t = float(t); omg = float(omg)
        T.append(t)
        omega.append(omg)
    return np.array(omega), np.array(T)


def fd(E, mu, T):
    # Fermi-Dirac distribution function
    # IN: eV, K
    # OUT: (no dim)
    k = 1.38062*10**-23 #J/K
    k = k/(1.602*10**-19) #eV/K
    if T == 0.0:
        temp = []
        i = 0
        for e in E:
            if e < mu: temp.append(1.0)
            elif abs(e-mu) < 10**-20: temp.append(0.5)
            else: temp.append(0.0)
            i += 1
        return np.array(temp)    
    else: return (1+np.exp((E-mu)/(k*T)))**-1


def fd1(E, mu, T):
    # 1st derivative of Fermi-Dirac distribution function
    # IN: eV, K
    # OUT: eV^-1 
    k = 1.38062*10**-23 #J/K
    k = k/(1.602*10**-19) #eV/K
    FD1 = -1./(k*T) * (1+np.exp((E-mu)/(k*T)))**-2 * np.exp((E-mu)/(k*T))
    # numerical error treatment
    index = 0    
    for prob in FD1:
        if str(prob) == 'nan': FD1[index] = 0.0
        elif prob > 0: FD1[index] = 0.0
        index += 1
    return FD1


def be(omega, Temp):
    # Bose-Einstein distribution function
    # IN: cm-1, K
    # OUT: (no dim)
    k = 1.38062*10**-23 #J/K
    k = k/(1.602*10**-19) #eV/K 
    hb = 4.1356675*10**-15 / (2*pi) # eV.s
    beta = 1./(k*Temp) # eV-1
    omega = (0.02997961 * 10 ** 12) * omega
    #omega = (1./fu2cm1) * fu2THz * 10**12 * np.array(omega) # cm-1 --> s-1
    return 1. / ( np.exp(beta*hb*omega)-1 ) # used on cm-1 grid


def be1(omega, Temp):
    # 1st derivative of Bose-Einstein distribution function by Temperature
    # IN: cm-1, K
    # OUT: K^-1
    k = 1.38062*10**-23 # J/K
    k = k/(1.602*10**-19) # eV/K 
    hb = 4.1356675*10**-15 / (2*pi) # eV.s
    beta = 1./(k*Temp) # eV-1
    omega = (0.02997961 * 10 ** 12) * omega
    #omega = (1./fu2cm1) * fu2THz * 10**12 * np.array(omega) # cm-1 --> s-1
    A = hb * omega / (k * Temp**2)
    B = np.exp(beta*hb*omega) / ( np.exp(beta*hb*omega) - 1 )**2
    BE1 = A * B
    # numerical error treatment
    index = 0
    for prob in BE1:
        if str(prob) == 'nan': BE1[index] = 0.0
        index += 1
    return BE1


def five_point_diff(arr_x, arr_y, intp_mag=1):
    # interpolation
    arr_y = np.interp(np.linspace(arr_x[0], arr_x[-1],
                                  intp_mag*(len(arr_x)-1)+1),
                      arr_x, arr_y)
    arr_x = np.linspace(arr_x[0], arr_x[-1], intp_mag*(len(arr_x)-1)+1)
    h = arr_x[1] - arr_x[0] # same interval h
    # extrapolation using spline
    extrapolator = intp.UnivariateSpline(arr_x, arr_y, k=3)
    x_l1 = arr_x[0]  - (arr_x[1]  - arr_x[0])
    x_l2 = arr_x[0]  - (arr_x[2]  - arr_x[0])
    x_r1 = arr_x[-1] + (arr_x[-1] - arr_x[-2])
    x_r2 = arr_x[-1] + (arr_x[-1] - arr_x[-3])
    arr_xl = np.array([x_l2, x_l1]); arr_yl = extrapolator(arr_xl)
    arr_xr = np.array([x_r1, x_r2]); arr_yr = extrapolator(arr_xr)
    # expaned data
    arr_x = np.array( list(arr_xl) + list(arr_x) + list(arr_xr) )
    arr_y = np.array( list(arr_yl) + list(arr_y) + list(arr_yr) )
    # first derivative (5pt. formula)
    y_1 = []
    i = 2
    for x in arr_x:
        if i > len(arr_y)-3: break
        y_ = (arr_y[i-2]-8*arr_y[i-1]+8*arr_y[i+1]-arr_y[i+2])/(12.*h)
        y_1.append(y_)
        i += 1
    return np.array(y_1)


def simpson_38(arr_x, arr_y, intp_mag=5):
    # interpolation
    arr_y = np.interp(np.linspace(arr_x[0], arr_x[-1], intp_mag*len(arr_x)), arr_x, arr_y)
    arr_x = np.linspace(arr_x[0], arr_x[-1], intp_mag*len(arr_x))
    # simpson, trapezoidal sum for the remainders
    S = 0; S_r = 0
    if len(arr_x) % 3 == 2:
        S_r = (arr_y[-2]+arr_y[-1])*(arr_x[-1]-arr_x[-2])/2.
        arr_x = arr_x[:-1]; arr_y = arr_y[:-1]
    elif len(arr_y) % 3 == 0:
        S_r = ((arr_y[-1]+arr_y[-2])*(arr_x[-1]-arr_x[-2])/2.) + \
              ((arr_y[-2]+arr_y[-3])*(arr_x[-2]-arr_x[-3])/2.)
        arr_x = arr_x[:-2]; arr_y = arr_y[:-2]
    # integration
    i = 0
    for x in arr_x:
        if i > len(arr_x)-2: break
        a = arr_x[i]; b = arr_x[i+3]; h = b - a
        S += (h/8)*(arr_y[i]+3*arr_y[i+1]+3*arr_y[i+2]+arr_y[i+3])
        i += 3
    return S + S_r


def Ln(E, T, n, Temp, Ef=0.0):
    # unit : eV**n
    e = 1.6022 * 10**-19 # C
    h = 4.1356675*10**-15 # eV.s
    # FD function
    #FD = fd(E, 0.0, Temp)
    # 1st derivative of FD function
    FD1 = fd1(E, Ef, Temp)
    #FD1 = five_point_diff(E, T, intp_mag=5)
    # Ln integral: CHECK
    #curnt = simpson_38(E, (E**n) * T * -FD1) # own simpson3/8
    curnt = intg.simps( ((E-Ef)**n) * T * -FD1, E-Ef, even='avg') # simpson in scipy
    return curnt


def S(E, T, Temp, Ef=0.0):
    # Calculate Seeback Coefficient / unit: Volt/K = J/C.K
    e = 1.6022 * 10**-19 # C
    L1 = Ln(E, T, 1, Temp, Ef) # eV**2
    L0 = Ln(E, T, 0, Temp, Ef) # eV
    return -1./(e*Temp) * (L1/L0) # C**-1 K**-1 eV


def k_el(E, T, Temp):
    # unit : eV / s.K
    e = 1.6022 * 10**-19 # C
    h = 4.1356675*10**-15 # eV.S
    # T conversion
    L2 = Ln(E, T, 2, Temp) # eV**3
    L1 = Ln(E, T, 1, Temp) # eV**2
    L0 = Ln(E, T, 0, Temp) # eV
    return 2./(h*Temp) * (L2 - (L1**2)/L0) # ev s**-1 K**-1 


def ZT_el(E, T, Temp):
    # unit : 1
    L2 = Ln(E, T, 2, Temp)
    L1 = Ln(E, T, 1, Temp)
    L0 = Ln(E, T, 0, Temp)
    ZT_el_1 = ( L0 * L2 / (L1**2) ) - 1
    return ZT_el_1**-1


def k_thermal(omega, T, Temp):
    # unit : eV / s.K
    hb = 4.1356675*10**-15 / (2*pi) # eV.s
    # BE` distribution
    BE1 = be1(omega, Temp) # be` on cm-1 grid
    # cm-1 to s-1
    omega = (0.02997961 * 10 ** 12) * omega # cm-1 --> s-1
    curnt = intg.simps( (1./(2*pi)) * hb * omega * BE1 * T, omega, even='avg') # eV s-1 K-1
    return curnt


def thermal_broadening_function(E, Temp, mu):
    k = 1.38062*10**-23 #J/K
    k = k/(1.602*10**-19) #eV/K
    return 1./(k*Temp) * (np.cosh(E/(2*k*Temp)))**-1


def conductance(E, Temp, mu, T_e):
    q = 1.6022 * 10**-19 # C
    h = 4.1356675*10**-15 # eV.s
    integral = intg.simps(T_e * thermal_broadening_function(E, Temp, mu), E, even='avg')
    #return (q**2) / h * integral
    return integral


def thermal_current(omega, T, Temp1, Temp2):
    # unit : eV/s
    hb = 4.1356675*10**-15 / (2*pi) # eV.s
    BE12 = be(omega, Temp1) - be(omega, Temp2) # (no dim)
    omega = (0.02997961 * 10 ** 12) * omega # cm-1 --> s-1
    J = intg.simps( (1./(2*pi)) * hb * omega * BE12 * T, omega, even='avg') # eV s-1
    return J
