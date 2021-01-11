import astropy.constants as ac
import numpy as np
import math
# planck     =  6.62606957e-27_rk         ! Planck constant in (non-SI) erg*second
# avogno     =  6.02214129e+23_rk         ! Avogadro constant
# vellgt     =  2.99792458e+10_rk         ! Speed of light constant in (non-SI) cm/second
# boltz      =  1.3806488e-16_rk          ! Boltzmann constant in (non-SI) erg/Kelvin
# bohr       =  0.52917720859_rk          ! bohr constant in Angstroms
# hartree    =  219474.6313705_rk         ! hartree in cm-1
# ev         =  8065.54465_rk             ! ev in cm-1
# uma        =  1.660538921e-24_rk        ! unified atomic mass unit [=mass(C-12)/12 ] in grams
# aston      =  planck/(8._rk*PI**2*vellgt*uma*1e-16_rk)  !rotational factor in cm-1 amu Ang^2
# todebye    =  2.541765_rk               ! a.u. in debye
# c2         =  1.43877736d0              ! second radiative constant NIST http://physics.nist.gov/cgi-bin/cuu/Value?c22ndrc
# R_         =  8.3144598d0               ! Molar gas constant R, J/mol/K
RT2LN2 = math.sqrt(2*math.log(2))
SQRT2PI = math.sqrt(2*math.pi)
PLANCK = ac.h.to('erg s').value
AVGNO = ac.N_A.value
SPDLIGT = ac.c.to('cm/s').value
KBOLTZ = ac.k_B.to('erg/K').value
c2 = PLANCK*SPDLIGT/KBOLTZ
PI = np.pi