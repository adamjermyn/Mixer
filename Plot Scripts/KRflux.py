from numpy import arctan
from numpy import sin
from numpy import cos
import numpy as np

# KR theory
aKR=2. # Anisotropy parameter

# KR theory
def psi0(o): # Psi_ortho (KR 1993)
    return 15./128./o**4*(o**2-21.-8.*o**2/(1.+o**2)+(21.+14.*o**2+o**4)/o*arctan(o))
def psi1(o): # Psi_parallel (KR 1993)
    return 15./32./o**4*(21.-3.*o**2+4.*o**4/(1.+o**2)-(21.+4.*o**2-o**4)/o*arctan(o))
# Lambda effects (KR 2005)
def J0(o): 
    return 1./2./o**4*(9.-2.*o**2/(1.+o**2)-(o**2+9.)/o*arctan(o))
def J1(o):
    return 1./2./o**4*(45.+o**2-4.*o**2/(1.+o**2)+(o**4-12.*o**2-45.)/o*arctan(o))
def I0(o):
    return 1./4./o**4*(-19.-5./(1.+o**2)+(3.*o**2+24.)/o*arctan(o))
def I1(o):
    return 3./4./o**4*(-15.+2.*o**2/(1.+o**2)+(3.*o**2+15.)/o*arctan(o))
def V0(o):
    return J0(o)+aKR*I0(o)
def H1(o):
    return J1(o)+aKR*I1(o)

# Heat flux anisotropy (KPR 1994)
def phi0(o):
    return 3./4./o**2*(1.+(o**2-1.)/o*arctan(o))
def phi1(o):
    return 3./4./o**2*(-3.+(o**2+3.)/o*arctan(o))

# Fluxes
def flux(tt,om,tag,mode):
    ost=2.*om
    ff = 0
    if tag=='13':
        if mode=='Omr':
            ff=psi0(ost)+cos(tt)**2*psi1(ost)
        if mode=='Omt':
            ff=-sin(tt)*cos(tt)*psi1(ost)
        if mode=='zero':
            ff=-V0(ost)+cos(tt)**2*H1(ost)
    elif tag=='23':
        if mode=='Omr':
            ff=-sin(tt)*cos(tt)*psi1(ost)
        if mode=='Omt':
            ff=psi0(ost)+sin(tt)**2*psi1(ost)
        if mode=='zero':
            ff=-cos(tt)*sin(tt)*H1(ost)
    elif tag=='14':
        if mode=='Sr':
            ff=phi0(ost)+phi1(ost)*cos(tt)**2
        if mode=='St':
            ff=-cos(tt)*sin(tt)*phi1(ost)
    elif tag=='24':
        if mode=='St':
            ff=phi0(ost)+phi1(ost)*sin(tt)**2
        if mode=='Sr':
            ff=-cos(tt)*sin(tt)*phi1(ost)            
    if tag[1]=='4':
        ff=-1./12.*ff
    else:
        ff=-1./15.*ff
    return 9*ff

# theta is co-latitude. Omega is rotation period.
# Two usage examples:
#  <v_r v_phi>=<u'^2>*  (flux(theta,Omega*tau,'13','Omr') *dOmega_dlnr*tau
#                      + flux(theta,Omega*tau,'13','Omt') *dOmega_dtheta*tau
#                      + flux(theta,Omega*tau,'13','zero')*Omega*tau)
#
#  <v_r r_r>= <u'^2>tau/3*(flux(theta,omega_over_N,'14','Sr') *entHat_r
#                      +   flux(theta,omega_over_N,'14','St') *entHat_t)
# 
# I think more or less <u'^2>~ (N lambda)**2 and N~1/tau
# tag combines correlations between vr-1, vtheta-2, vphi-3, or r_r-4.
# entHat_r and _t are the r theta components of your entHat entropy
# gradient unit vector.
