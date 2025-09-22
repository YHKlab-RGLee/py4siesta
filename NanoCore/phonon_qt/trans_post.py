#
# Post processing
#

import numpy as np
import matplotlib.pyplot as plt
from .unit import *


def mycmp(a, b): return cmp(a[0], b[0])


def write_oneD_one_line(file_name, array_x, array_y, label_x, label_y):
    oneD = open(file_name, 'w')
    if len(array_x) != len(array_y):
        raise ValueError("array_x and array_y should have same length.")
    head1 = "#NumField:1\n"
    head2 = "#LabelX:%s,LabelY:%s\n" % (label_x, label_y)
    oneD.write(head1); oneD.write(head2)
    head3 = "#Field1:%s,NumPoint:%i\n" % (file_name, len(array_x))
    oneD.write(head3)
    i_y = 0
    for x in array_x:
        line = "%12.8f    %12.8f\n" % (x, array_y[i_y]); i_y +=1
        oneD.write(line)
    oneD.close()


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


def be(omega, Temp):
    # Bose-Einstein distribution function
    # IN: cm^-1, K
    # OUT: (no dim)
    beta = 1./(k_boltz*Temp) # eV^-1
    omega = (0.02997961 * 10 ** 12) * omega
    #omega = (1./fu2cm1) * fu2THz * 10**12 * np.array(omega) # cm-1 --> s-1
    return 1. / ( np.exp(beta * hbar * omega)-1 ) # used on cm-1 grid


def be1(omega, Temp):
    # 1st derivative of Bose-Einstein distribution function by Temperature
    # IN: cm^-1, K
    # OUT: K^-1
    beta = 1./(k_boltz * Temp) # eV^-1
    omega = (0.02997961 * 10 ** 12) * omega
    #omega = (1./fu2cm1) * fu2THz * 10**12 * np.array(omega) # cm-1 --> s-1
    A = hbar * omega / (k_boltz * Temp**2)
    B = np.exp(beta * hbar * omega) / ( np.exp(beta * hbar * omega) - 1 )**2
    BE1 = A * B
    # numerical error treatment
    index = 0
    for prob in BE1:
        if str(prob) == 'nan': BE1[index] = 0.0
        index += 1
    return BE1


class VisTrans(object):

    __slots__ = ['omega_grid', 'q_mesh', 'T_omega_q', 'W_omega_q']

    def __init__(self, T_all_, q_mesh, omega_grid):
        self.omega_grid = omega_grid
        self.q_mesh = q_mesh
        # empty container
        self.T_omega_q = np.zeros( (len(omega_grid), len(q_mesh)) )
        self.W_omega_q = np.zeros( (len(omega_grid), len(q_mesh)) )
        # for all q-trans,
        for T_all in T_all_:
            T_all.sort(key=lambda x: x[0])
            # for each (q, omega),
            for omega, t, qpt, wei in T_all:
                index_omega = list(omega_grid).index(omega)
                index_qpt   = list(q_mesh).index(qpt)
                self.T_omega_q[index_omega, index_qpt] = t
                self.W_omega_q[index_omega, index_qpt] = wei

    def get_qavtrns(self, factor=vasp2cm1, is_savefig=0):
        # q-averaged transmission
        i = 0
        q_av_T = [] 
        while i < len(self.omega_grid):
            j = 0
            t = 0.; wsum = 0
            while j < len(self.q_mesh):
                tq = self.T_omega_q[i,j]
                w  = self.W_omega_q[i,j]
                t += w * tq; wsum += w
                j += 1
            q_av_T.append(t/wsum)
            i += 1
        if is_savefig:
            fig = plt.figure(1, figsize=(9,4))
            fig1 = fig.add_subplot(111)
            fig1.plot(self.omega_grid*factor, q_av_T)
            fig.savefig('q_averaged_transmission.png')
        write_oneD_one_line('q_averaged_transmission.txt', self.omega_grid*factor, q_av_T, 'omega[cm-1]', 'T[omega]')
        return self.omega_grid*factor, q_av_T

    def get_qtrns(self, q=[0.,0.,0.], factor=vasp2cm1, is_savefig=0):
        index = list(self.q_mesh).index(q)
        T = self.T_omega_q[:,index].T
        if is_savefig:
            fig = plt.figure(1, figsize=(9,4))
            fig1 = fig.add_subplot(111)
            fig1.plot(self.omega_grid*factor, T)
            fig.savefig('q_transmission.png')
        return self.omega_grid*factor, T

    def get_qdivtrns(self, direction=None, factor=vasp2cm1, is_savefig=0):
        if is_savefig:
            fig = plt.figure(1, figsize=(9,4))
            fig1 = fig.add_subplot(111)
            i_q = 0
            for T in self.T_omega_q.T:
                fig1.plot(self.omega_grid*factor, T, label='%s' % self.q_mesh[i_q])
                write_oneD_one_line('q_transmission_%f_%f_%f.txt' % tuple(self.q_mesh[i_q]),
                                    self.omega_grid*factor, T, 'omega[cm-1]', 'T[omega]')
                i_q += 1
            fig.savefig('q_transmission.png')
        return self.omega_grid*factor, self.T_omega_q.T

    def get_thermal_conductivity(self, Temp=300.0):
        # unit : eV / s.K
        BE1 = be1(self.omega_grid, Temp) # be` on cm-1 grid
        omega = (0.02997961 * 10 ** 12) * self.omega_grid # cm-1 --> s-1
        #return intg.simps( (1./(2*pi)) * hbar * omega * BE1 * T, omega, even='avg') # eV s-1 K-1
        return simpson_38( omega, (1./(2*pi)) * hbar * omega * BE1 * self.get_qavtrns()[-1]) # eV s-1 K-1

    def get_thermal_current(self, Temp1, Temp2):
        # unit : eV / K
        BE12 = be(self.omega_grid, Temp1) - be(self.omega_grid, Temp2) # (no dim)
        omega = (0.02997961 * 10 ** 12) * self.omega_grid # cm-1 --> s-1
        #return intg.simps( (1./(2*pi)) * hbar * omega * BE12 * T, omega, even='avg') # eV s-1
        return simpson_38( omega, (1./(2*pi)) * hbar * omega * BE12 * self.get_qavtrns()[-1]) # eV s-1
