import numpy as np

def fmb10_mah(aexp, z0mass=1e12, xi_in=np.logspace(-4, 0), reds=False):

    """
    Integration of the mean mass growth rate given in Fakhouri, Ma, \&
    Boylan-Kolchin 2010 (Eqn. 2)
    """

    # if reds=False (default), use expansion factor. If reds=True, use redshift
    # rather than expansion factor.
    if reds:
        myaexp=1./(1+aexp)
    else:
        myaexp=aexp.copy()
        
    A=46.1
    M0=1e12
    alpha1=1.1
    zpar1=1.11
    H0=70.2

    alpha2=alpha1-1.0
    # fz=zpar1*(1+reds)-(zpar1-1.0)*log(1+reds)-zpar1
    fz=zpar1/myaexp+(zpar1-1.0)*np.log(myaexp)-zpar1
    prefac=A/(10*H0*(M0/1e12))
    return z0mass*(1+prefac*(z0mass/M0)**alpha2*fz)**(-1./(alpha1-1))
