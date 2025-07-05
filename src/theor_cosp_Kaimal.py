# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 12:30:08 2022

@author: Ariane Faures

Expression of the theoretical cospectra according to Kaimal 1972
This can be used as the cospectra of referece as it is computed for a "perfect"
site : fully flat and homogeneous.

Inputs:
    - Stability class (zL)
    - Natural frequency (nf)
    - Normalised frequency (kf)

Outputs:
    - Theoretical Kaimal cospectral as fct of stability
"""


def Kaimal_cosp(nf, kf, zL):

    if zL > 0:  # Stable conditions
        kf0 = 0.23*(1 + 6.4*zL)**0.75
        kaimal_cosp=0.81*kf/kf0/(1 + 1.5*(kf/kf0)**(2.1))/nf
        # reynolds stresses
        Au = 0.124 * ((1 + 7.9 * zL)**0.75)
        Bu = 2.34 * (Au**(-1.1))
        Cospwu = kf / (nf * (Au + Bu * kf**2.1))
    else:  # Unstable conditions
        if kf <= 1:
            kaimal_cosp=11 *kf/(1 + 13.3*kf)**(7/4)/nf
        else:
            kaimal_cosp=4 *kf/(1 +3.8*kf)**(7/3)/nf

    return kaimal_cosp


def Kaimal_cosp_EP(fnorm, zL):  # According to the EddyPro source code (G. Fratini)
    if zL > 0:  # Stable conditions
        Ak = 0.284 * ((1 + 6.4 * zL)**0.75)
        Bk = 2.34  * (Ak**(-1.1))
        kaimal = fnorm / (Ak + Bk * fnorm**2.1)
    else : # Unstable conditions
        if fnorm <= 0.54:
            kaimal = 12.92 * fnorm / ((1 + 26.7 * fnorm)**1.375)
        else:
            kaimal = 4.378 * fnorm / ((1 +  3.8 * fnorm)**2.4)

    return kaimal
