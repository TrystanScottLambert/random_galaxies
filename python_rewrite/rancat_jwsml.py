"""
rancat_jswml module mimicing the fortran one
"""

from data_classes import SurveySpec, ModelParameters
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM

from histograms import wcic_histogram
from rancat_subroutines import tabulate_bin_pvolumes, initilize_bin_dataframe


redshifts = np.random.random(1000)
mags = np.random.random(1000) * 5 + 14

survey = SurveySpec(min(redshifts), max(redshifts), 19.00, 0.006, 0.006)
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
number_z_bins = 40
number_lf = 25


z_bins_limits = np.linspace(np.min(redshifts), np.max(redshifts), number_z_bins)
z_bins = initilize_bin_dataframe(z_bins_limits, cosmo, survey)
z_bins['pvol'] = tabulate_bin_pvolumes(
    z_bins_limits, survey.survey_fractional_area, 0.1, 0.3, cosmo)
z_bins['n'] = wcic_histogram(redshifts, z_bins['delta'], z_bins['mid'])


par = ModelParameters(a=0, u=0, mu=0, p_post=0)
par_previous = ModelParameters(a = 1e30, u=1e30, mu=1e30, p_post=0)
par_last = par

delta = np.ones(number_z_bins)
lf = np.zeros(number_lf)
lf_last = np.zeros(number_lf)


converged_pe = False
converged = False

converged_du = False
i_branch = 0
i_search = 0

iter = 1
N_SEARCH = 1000
NITERMAX = 400

tolE = 1.0e-3
tolR = 1.0e-4
tolL = 1.0e-6
tolD = 1.0e-4
tolM = 1.0e-4

while True:
    if i_search > N_SEARCH:
        raise TimeoutError("No Convergence, Try increasing N_SEARCH")

    nitermin = min(2+iter, NITERMAX)
    while not converged and iter < NITERMAX:
        if (iter > nitermin) and (abs(par.a - par_last.a) < tolE) and (abs(par.u - par_last.u) < tolE):
            if not converged_pe:
                print(f"Density and Luminosity converged to accuracy of err(a) = {abs(par.a - par_last.a)} and err(u) = {abs(par.u - par_last.u)}")
            converged_pe = True
        
        if (converged_pe and max(abs(delta - z_bins['delta'])) < tolD) and (max(abs(lf - lf_last)) < tolL):
            converged = True
            print('Density bins converged with an accuray: ', max(abs(delta - z_bins['delta'])))
            print("LF bins converged to an accuracy: ", max(abs(lf - lf_last)))

        
        z_bins['pvol'] = tabulate_bin_pvolumes(
            z_bins_limits, survey.survey_fractional_area, 0.1, par.a, cosmo)
        