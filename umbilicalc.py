# -*- coding: utf-8 -*-
"""
Written in April, 2020 by D. J. Wilke.
ver0.0
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
 _    _  __  __  ____  _____  _      _____      _____            _       _____ 
| |  | ||  \/  ||  _ \|_   _|| |    |_   _|    / ____|    /\    | |     / ____|
| |  | || \  / || |_) | | |  | |      | | ___ | |        /  \   | |    | |     
| |  | || |\/| ||  _ <  | |  | |      | ||___|| |       / /\ \  | |    | |     
| |__| || |  | || |_) |_| |_ | |___  _| |_    | |____  / ____ \ | |___ | |____ 
 \____/ |_|  |_||____/|_____||_____||_____|    \_____|/_/    \_\|_____| \_____|

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
 
This script takes as input the cord length (L) and coiling width (w), both in 
cm, in addition to the number of coils (N) measured from the human umbilical 
cord. It computes a range of parameters relevant to the estimation of the total 
vessel pressure drop in both the vein and artery of the umbilical cord. The 
procedure follows that outlined in Wilke et. al. (2020) (in submission), and 
calculates the umbilical pressure index (PX) presented in that work. For cords 
with distinctly varying sections the script can be run for each part, with the 
pressure drop and PX simply the sum of constituents.


Inputs:
- As a minimum, N, L and w are required as input, however, further details on 
the vessel size, as well as the blood flow-rate, density and viscosity may be 
added if known*. 

- The boolean flag 'output_csv', when set to 'True', outputs data to 
'uc_output.csv' within the working directory.


Calculation:
- See Wilke et. al. (2020) for full details. The pressure calculations are 
completed using the formulae of Liu and Masliyah (1993) for vessels satisfying 
gamma < 0.1. For other vessels (with Re = 100), the pressure gradient is 
determined using a radial basis interpolation of additional numerics, contained 
within the data file, 'numerics.dat'. Note that this script uses the 'scipy' 
package 'interpolate' which differs slightly to the work presented in 
Wilke et. al. (2020), which was completed using MATLAB. 

*Note that additional blood-flow information may require the re-calculation of 
the Reynolds number, and if the vessel has gamma >= 0.1 there is currently no 
solution (see note under assumptions). The forumulae from Liu & Masliyah (1993)
also has constraints on the Dean number (Dn), however, these are unlikely to be 
breached by flow within the umbililcal vessels and so this constraint is not 
enforced here.


Dependencies:
- The script requires the 'numpy', 'pandas', 'scipy' and 'math' packages. 
Computation for cord vessels with gamma >= 0.1, Re = 100, requires the 
'numerics.dat' data file.


# -----------------------------------------------------------------------------
# References ------------------------------------------------------------------
# -----------------------------------------------------------------------------

# By the authors:

D. J. Wilke, J. P. Denier, T. Y. Khong, T. W. Mattner, Pressure and
flow in the umbilical cord, Journal of Biomechanics (79) (2018) 78–87.
doi:10.1016/j.jbiomech.2018.07.044.

D. J. Wilke, Pressure and flow within the umbilical vessels, Ph.D thesis,
The University of Adelaide (2016).



# From the literature:

G. Acharya, T. Wilsgaard, G.K.R. Berntsen, J.M. Maltau, T. Kiserud, Reference
ranges for serial measurements of blood velocity and pulsatility index at the
intra-abdominal portion, and fetal and placental ends of the umbilical artery, 
Ultrasound in Obstetrics and Gynecology 26 (2005a) 162-169. 

G. Acharya, T. Wilsgaard, G.K.R. Berntsen, J.M. Maltau, T. Kiserud, Reference 
ranges for umbilical vein blood flow in the second half of pregnancy based on 
logitudinal data, Prenatal Diagnosis 25 (2005b) 99-111.

C. Barbieri, J. G. Cecatti, F. G. Surita, E. F. Marussi, J. V. Costa, Sono-
graphic measurement of the umbilical cord area and the diameters of its
vessels during pregnancy, Journal of Obstetrics and Gynaecology 32 (2012)
230–236. doi:10.3109/01443615.2011.647129.

C. G. Caro, T. J. Pedley, C. W. Schroter, W. A. Seed, The mechanics of
the circulation, Oxford University Press, 1978.

A. Guettouche, J. C. Challier, Y. Ito, C. Papapanayotou, Y. Cherruault, 
A. Azancot-Benisty, Mathematical modeling of the human fetal arterial blood 
circulation, International Journal of Biomedical Computing 31 (1992) 127–139.

A. D. Kaplan, A. J. Jaffa, I. E. Timor, D. Elad, Hemodynamic analysis
of arterial blood flow in the coiled umbilical cord, Reproductive Sciences
17 (3) (2010) 258–268. doi:10.1177/1933719109351596.

S. Liu, J. H. Masliyah, Axially invariant laminar flow in helical pipes with
a finite pitch, Journal of Fluid Mechanics 251 (1993) 315–353.

L. Raio, F. Ghezzi, E. Di Naro, R. Gomez, M. Franchi, M. Mazor,
H. Brühwiler, Sonographic measurement of the umbilical cord and fetal
anthropometric parameters, European Journal of Obstetrics, Gynecology,
and Reproductive Biology 83 (2) (1999) 131–135.

M. Tahmasebi, R. Alighanbari, Evaluation of umbilical cord thickness,
cross-sectional area, and coiling index as predictors of pregnancy outcome,
Indian Journal of Radiology and Imaging 21 (3) (2011) 195–198.

S. L. Waters, C. Guiot, Flow in an elastic tube subject to prescribed forcing:
a model of umbilical venous flow, Journal of Theoretical Medicine 3 (4)
(2001) 287–298. doi:10.1080/10273660108833081.

"""

# -----------------------------------------------------------------------------
# Load packages ---------------------------------------------------------------
# -----------------------------------------------------------------------------

import math      
import numpy as np 
import pandas as pd 

from scipy.interpolate import Rbf

# -----------------------------------------------------------------------------
# User input ------------------------------------------------------------------
# -----------------------------------------------------------------------------

N = 10 # Number of coils (default 10)
L = 50 # Cord length in cm (default 50)
w = 1.6 # Coiling width in cm (default 1.6)

output_csv = True # Flag to output data.

# -----------------------------------------------------------------------------
# Assumptions -----------------------------------------------------------------
# -----------------------------------------------------------------------------

R_V = 0.35 # Vessel radius (vein) in cm, Waters et. al. (2001)
R_A = 0.2 # Vessel radius (artery) in cm, Kaplan et. al. (2010)
Q_V = 265 # Blood flow-rate (vein) in mL/min, Acharya et. al. (2005b)
Q_A = 145 # Blood flow-rate (artery) in mL/min, Acharya et. al. (2005a)
rho = 1060 # Blood density (kg/m^3), Guettouche et. al. (1992) 
mu = 4e-3 # Blood viscosity, Caro et. al. (1978)

# The Reynolds number can be determined from the above via:
# Re = (Q*rho)/(np.pi*R*mu). For the reference case, Re_V = 106, Re_A = 102. To
# simplify we take Re = 100 in this work. Note that if Re is not equal to 100
# the additional numerics for gamma >= 0.1 cannot be used.
Re = 100 # Reynolds number, Wilke (2016). 

# For the reference cord (to determine indices):
dP_dim_R_V = 1.9011605673665974 # Pressure drop (vein), Wilke et. al. (2020) 
dP_dim_R_A = 9.8291206851390350 # Pressure drop (artery), Wilke et. al. (2020) 

# -----------------------------------------------------------------------------
# Calculate -------------------------------------------------------------------
# -----------------------------------------------------------------------------
print('\nRunning with...')
print('N = ',  "%.1f" %N)
print('L = ',  "%.3f" %L, '(cm)')
print('w = ',  "%.3f" %w, '(cm)')
print('\n')

UCI = N/L # Umbilical coiling index (UCI)

R = np.array([R_V, R_A])
Q = np.array([Q_V, Q_A])/6e+7 # Convert to m^3/s
dP_dim_R = np.array([dP_dim_R_V, dP_dim_R_A])

# Non-dimensionalise:
L = L/R
w = w/R

# Compute pressure:
if N == 0.0:
    # Pressure drop for straight tube:
    Rh = np.array(['N/A', 'N/A']) # Helical radius of vessel centreline
    Ph = np.array(['N/A', 'N/A']) # Helical pitch of vessel centreline
    tau = np.array(['N/A', 'N/A']) # Helical torsion
    kappa = np.array(['N/A', 'N/A']) # Helical curvature
    Dn = np.array(['N/A', 'N/A']) # Deans parameter
    gamma = np.array(['N/A', 'N/A']) # Gamma parameter
    Lh = L # Helical arclength
    print('Calculating for a straight cord.')
    dP = 8.0*L # Vessel pressure drop (non-dimensional)
else:
    # Helical Section:
    Rh = w/2.0 - 1.0 
    Ph = L/N 
    tau = (Ph/(2*math.pi))/(Rh*Rh+(Ph/(2*math.pi))**(2))
    kappa = Rh/(Rh*Rh+(Ph/(2*math.pi))**(2)); 
    Dn = 2*Re*kappa**(0.5); 
    gamma = tau/(kappa**(3/4)*(2*Re)**0.5)
    
    if ( gamma[0] >= 0.1 ) | ( gamma[1] >= 0.1 ):
        # load in interpolant data and setup interpolation scheme:
        dat = pd.read_csv('numerics.dat', sep='\s+',header=None, names=["kappa", "tau", "dPdZ"])
        interp = pd.DataFrame(dat)
        rbfi = Rbf(interp['kappa'], interp['tau'], interp['dPdZ'])  # radial basis function interpolator instance

    fRe = np.array([0.0, 0.0])
    Lh = np.array([0.0, 0.0])
    dP = np.array([0.0, 0.0])
    
    if gamma[0] < 0.1:
        print('gamma (V) = ', "%.3f" %gamma[0])
        print('==> Using the forumluae of Liu & Masliyah (1993) for the vein.\n')
        fRe[0] = (16.0 + (((2.0*Re)**0.5)*0.378 + 12.1*(kappa[0]**(-0.5))*Dn[0]**(-0.5))*tau[0]**2.0)*(1.0 + ((0.0908+0.0233*kappa[0]**(0.5))*Dn[0]**(0.5) - 0.132*kappa[0]**(0.5) + 0.37*kappa[0] - 0.2)/(1+(49/Dn[0]))) # Flow friction factor from Liu and Masliyah (1993)
        Lh[0] = 2*math.pi*(N)*(Rh[0]*Rh[0]+(Ph[0]/(2*math.pi))*(Ph[0]/(2*math.pi)))**(0.5)
        dP[0] = (fRe[0]/2)*Lh[0]
    elif Re == 100:
        # Need to use interpolant
        print('gamma (V) = ', "%.3f" %gamma[0])
        print('==> Using interpolant from Wilke et. al. (2020) for the vein.\n')
        Lh[0] = 2*math.pi*(N)*(Rh[0]*Rh[0]+(Ph[0]/(2*math.pi))*(Ph[0]/(2*math.pi)))**(0.5)
        dPdZ = rbfi(kappa[0], tau[0]) # interpolated
        dP[0] = -1.0*dPdZ*Lh[0]
    else:
        print('gamma (V) = ', "%.3f" %gamma[0])
        print('==> Re not equal to 100, cannot compute pressure for the vein.\n')
        Lh[0] = np.nan
        dP[0] = np.nan
        
    if gamma[1] < 0.1:
        print('gamma (A) = ', "%.3f" %gamma[1])
        print('==> Using the forumluae of Liu & Masliyah (1993) for the artery.\n')
        fRe[1] = (16.0 + (((2.0*Re)**0.5)*0.378 + 12.1*(kappa[1]**(-0.5))*Dn[1]**(-0.5))*tau[1]**2.0)*(1.0 + ((0.0908+0.0233*kappa[1]**(0.5))*Dn[1]**(0.5) - 0.132*kappa[1]**(0.5) + 0.37*kappa[1] - 0.2)/(1+(49/Dn[1]))) # Flow friction factor from Liu and Masliyah (1993)
        Lh[1] = 2*math.pi*(N)*(Rh[1]*Rh[1]+(Ph[1]/(2*math.pi))*(Ph[1]/(2*math.pi)))**(0.5) # Helical arclength
        dP[1] = (fRe[1]/2)*Lh[1]
    elif Re == 100:
        # Need to use interpolant
        print('gamma (A) = ', "%.3f" %gamma[1])
        print('==> Using interpolant from Wilke et. al. (2020) for the artery.\n')
        Lh[1] = 2*math.pi*(N)*(Rh[1]*Rh[1]+(Ph[1]/(2*math.pi))*(Ph[1]/(2*math.pi)))**(0.5)
        dPdZ = rbfi(kappa[1], tau[1]) # interpolated
        dP[1] = -1.0*dPdZ*Lh[1]
    else:
        print('gamma (A) = ', "%.3f" %gamma[1])
        print('==> Re not equal to 100, cannot compute pressure for the artery.\n') 
        Lh[1] = np.nan
        dP[1] = np.nan

dP_dim = (mu*mu*Re*dP*0.0075)/(rho*(R/100)*(R/100)) # Vessel pressure drop in mmHg
PX = dP_dim/dP_dim_R # Pressure index (PX)
dP_s = (mu*mu*Re*8.0*Lh*0.0075)/(rho*(R/100)*(R/100)) # Vessel pressure drop in mmHg through a straightened vessel of the same arclength
CPA = dP_dim/dP_s # Coil pressure anomaly

# -----------------------------------------------------------------------------
# Output ----------------------------------------------------------------------
# -----------------------------------------------------------------------------

print('Pressure drop (V) = ',  "%.3f" %dP_dim[0], '(mmHg)')
print('Pressure index, PX (V) = ',  "%.3f\n" %PX[0])

print('Pressure drop (A) = ',  "%.3f" %dP_dim[1], '(mmHg)')
print('Pressure index, PX (A) = ',  "%.3f\n" %PX[1])


if output_csv == True:
    # Save data to .csv
    filename = 'uc_output.csv'
    data = {'param' : ['N','L','w','Re','R','Q','rho','mu','Lh','tau','kappa','gamma','dP','dP_R','PX','CPA','UCI'],
            'vein'  :  [N, L[0]*R[0], w[0]*R[0], Re, R_V, Q_V, rho, mu, Lh[0]*R[0], tau[0], kappa[0], gamma[0], dP_dim[0], dP_dim_R_V, PX[0],CPA[0],UCI],
            'artery':  [N, L[1]*R[1], w[1]*R[1], Re, R_A, Q_A, rho, mu, Lh[1]*R[1], tau[1], kappa[1], gamma[1], dP_dim[1], dP_dim_R_A, PX[1],CPA[1],UCI],
            'info':  ['no. of coils','cord length (cm)','coiling width (cm)','flow Reynolds number','vessel cross-sectional radius (cm)','cross-sectional blood flow-rate (mL/min)','blood density (kg/m^3)','blood viscosity (kg/ms)','vessel arclength (cm)','vessel torsion','vessel curvature','gamma parameter','vessel pressure drop (mmHg)','reference vessel pressure drop (mmHg)', 'pressure index', 'coil pressure anomaly: a ratio of pressure drop in the coiled and straightened vessels','umbilical coiling index']
        }
    df = pd.DataFrame(data, columns = ['param','vein','artery','info'])
    df.to_csv(filename, encoding='utf-8', index=False, header=True, float_format='%.3f')
    print('Output written to "',"%s" %filename,'"')

