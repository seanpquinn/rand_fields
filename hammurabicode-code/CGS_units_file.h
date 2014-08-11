/* CGS-UNITS */
#ifndef CGS_UNITS_H
#define CGS_UNITS_H

#define CGS_U_erg 1.e-2
#define CGS_U_cm  1.e-3
#define CGS_U_sec 1.
#define CGS_U_Kelvin 1. /* set always to 1!! */

#define CGS_U_gram (CGS_U_erg*CGS_U_sec*CGS_U_sec/(CGS_U_cm*CGS_U_cm))
#define CGS_U_Joule (1.e7*CGS_U_erg)
#define CGS_U_Watt  (CGS_U_Joule/CGS_U_sec)

#define CGS_U_pi 3.14159265358979
#define CGS_U_arcmin (2.*CGS_U_pi/(360.*60.))
#define CGS_U_arcsec (2.*CGS_U_pi/(360.*60.*60.))
#define CGS_U_sterad (4.* CGS_U_pi)
#define CGS_U_esu sqrt(CGS_U_erg*CGS_U_cm)
#define CGS_U_Coulomb 2997924579.99957*CGS_U_esu
#define CGS_U_Year (365.25*24.*3600.*CGS_U_sec)
#define CGS_U_m (1.e2 * CGS_U_cm)
#define CGS_U_km  (1.e5 * CGS_U_cm)
#define CGS_U_pc (3.0856775806e+18 * CGS_U_cm)
#define CGS_U_kpc (3.0856775806e+21 * CGS_U_cm)
#define CGS_U_Mpc (1.e3 * CGS_U_kpc)
#define CGS_U_Gpc (1.e6 * CGS_U_kpc)
#define CGS_U_eV  (1.60217733e-12 * CGS_U_erg)
#define CGS_U_keV (1.e3 * CGS_U_eV)
#define CGS_U_MeV (1.e6 * CGS_U_eV)
#define CGS_U_GeV (1.e9 * CGS_U_eV)
#define CGS_U_Hz  (1./CGS_U_sec)
#define CGS_U_kHz (1.e3*CGS_U_Hz)
#define CGS_U_MHz (1.e6*CGS_U_Hz)
#define CGS_U_GHz (1.e9*CGS_U_Hz)
#define CGS_U_barn  (1.e-24*CGS_U_cm*CGS_U_cm)
#define CGS_U_mbarn (1.e-27*CGS_U_cm*CGS_U_cm)
#define CGS_U_Jansky (1.e-23*CGS_U_erg/CGS_U_sec /(CGS_U_cm*CGS_U_cm) /CGS_U_Hz)
#define CGS_U_Gauss   sqrt(CGS_U_erg/(CGS_U_cm*CGS_U_cm*CGS_U_cm))
#define CGS_U_muGauss (1.e-6 * CGS_U_Gauss)
#define CGS_U_ccm  (CGS_U_cm*CGS_U_cm*CGS_U_cm)

#define CGS_U_C_light   (2.99792458e+10* CGS_U_cm / CGS_U_sec)
#define CGS_U_h_planck   (CGS_U_hq*2.*CGS_U_pi) /* Planck constant */
#define CGS_U_hq  (1.05457266e-27*CGS_U_erg*CGS_U_sec)  /* Plancks constant/(2pi) */
#define CGS_U_qe  (4.8032068e-10*CGS_U_esu)
#define CGS_U_MEC2  (0.51099907e-3 * CGS_U_GeV)     /* Electron Mass times c^2 */
#define CGS_U_MPC2  (938.272310e-3 * CGS_U_GeV)     /* Proton   Mass times c^2 */
#define CGS_U_ME    (CGS_U_MEC2 /(CGS_U_C_light*CGS_U_C_light))
#define CGS_U_MP    (CGS_U_MPC2 /(CGS_U_C_light*CGS_U_C_light))
#define CGS_U_sigmaT (0.66524616*CGS_U_barn)
#define CGS_U_kB     (1.380622e-16*CGS_U_erg/CGS_U_Kelvin)
#define CGS_U_re     (2.81794092e-13*CGS_U_cm)

#define CGS_U_h_50  1.
#define CGS_U_Ho    (50.* CGS_U_km /CGS_U_sec /CGS_U_Mpc * CGS_U_h_50)

#define CGS_U_MSol (1.98892e+33* CGS_U_gram)
#define CGS_U_LSol (3.826e+33*CGS_U_erg/CGS_U_sec)
#define CGS_U_RSol (6.9599e+10*CGS_U_cm)

#endif
