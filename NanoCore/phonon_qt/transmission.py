#
# SOLVE surface Green`s function
#

import numpy as np

def surf2green(k00, k11, k01, omega, tol = 10**-5):
    k00 = np.matrix(k00)
    k11 = np.matrix(k11)
    k01 = np.matrix(k01)
    Omega2 = np.matrix(omega**2 * np.diag(np.ones(len(k00))))
    s = np.matrix(k00)
    e = np.matrix(k11)
    alpha = np.matrix(k01)
    max_iter = 1000
    i_iter = 1
    while abs(alpha).max() > tol:
        #print i_iter, abs(alpha).max()
        g = (Omega2 - e)**-1
        beta = alpha.T
        s = s + alpha * g * beta
        e = e + alpha * g * beta + beta * g * alpha
        alpha = alpha * g * alpha
        i_iter += 1
        if i_iter > max_iter:
            print ("surface Green`s function is not converged.")
            return 'None'
    #print i_iter
    return (Omega2 - s)**-1

#
# TRANSMISSION 
#

def caroli(omega, kL00, kL11, kL01, kR00, kR11, kR01, VLC, VRC, kC, eta=(10**-6)*1j, return_aux=0):
    gL = surf2green(kL00, kL11, kL01, omega+eta)
    gR = surf2green(kR00, kR11, kR01, omega+eta)
    if gL == 'None' or gR == 'None':
        print ("surf. Green`s function is not converged.", omega)
        return 0.0
    VLC = np.matrix(VLC)
    VRC = np.matrix(VRC)
    kC = np.matrix(kC)
    SigmaL = VLC.T * gL * VLC
    SigmaR = VRC.T * gR * VRC
    #SigmaL = VLC.conj() * gL * VLC
    #SigmaR = VRC.conj() * gR * VRC
    m = len(kC)
    Omega2 = ( (omega)**2 * np.matrix(np.diag(np.ones(m))) )
    Gr = (Omega2 - kC - SigmaR - SigmaL)**-1
    GammaL = -2.0 * SigmaL.imag
    GammaR = -2.0 * SigmaR.imag
    if return_aux:
        return Gr, SigmaL, SigmaR
    GrGLGrGR = np.array(Gr * GammaL * Gr.conj().T * GammaR)
    t = GrGLGrGR.trace().real
    return t

