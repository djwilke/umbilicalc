This script (umbilicalc.py) computes a range of parameters relevant to the 
estimation of the total vessel pressure drop in both the vein and artery of 
the human umbilical cord. As a minimum it requires as input the cord length 
(L) and coiling width (w), both in cm, in addition to the number of coils (N).
The procedure follows that outlined in Wilke et. al. (2020) (in submission), 
and calculates the umbilical pressure index (PX) presented in that work. For 
cords with distinctly varying sections the script can be run for each part, 
with the pressure drop and PX simply the sum of constituents.

# Inputs:
- As a minimum, N, L and w are required as input, however, further details on 
the vessel size, as well as the blood flow-rate, density and viscosity may be 
added if known*. 
- The boolean flag 'output_csv', when set to 'True', outputs data to 
'uc_output.csv' within the working directory.

# Dependencies:
- The script requires the 'numpy', 'pandas', 'scipy' and 'math' packages. 
Computation for cord vessels with gamma >= 0.1, Re = 100, requires the 
'numerics.dat' data file.

# Calculation:
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

# References:
By the authors:

D. J. Wilke, J. P. Denier, T. Y. Khong, T. W. Mattner, Pressure and
flow in the umbilical cord, Journal of Biomechanics (79) (2018) 78–87.
doi:10.1016/j.jbiomech.2018.07.044.

D. J. Wilke, Pressure and flow within the umbilical vessels, Ph.D thesis,
The University of Adelaide (2016).

From the literature:

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
