"""
rancat_jswml module mimicing the fortran one
"""

import numpy as np

from histograms import wcic_histogram
from rancat_subroutines import tabulate_bin_pvolumes
from data_classes import SurveySpec, ModelParameters, BinnedZ


redshifts = np.random.random(1000)
mags = np.random.random(1000) * 5 + 14
survey = SurveySpec(min(redshifts), max(redshifts), 19.00, 0.006, 0.006)

number_z_bins = 40
number_lf = 25


z_bins = np.linspace(np.min(redshifts), np.max(redshifts), number_z_bins)
z_bin_midpoints = 0.5 * (z_bins[1:] + z_bins[:-1])
z_bin_deltas = np.ones(number_z_bins)
z_bin_n = np.zeros(number_z_bins)
z_bin_ran = np.ones(number_z_bins)
z_bin_vol = np.ones(number_z_bins)
z_bin_pvol = np.ones(number_z_bins)

volbin, pvolbin = tabulate_bin_pvolumes(z_bins, survey.solid_angle, 0.1, 0.3)
n_of_z = wcic_histogram(redshifts, np.ones(len(z_bins)), z_bins)


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
        
        if (converged_pe and max(abs(delta - z_bin_deltas)) < tolD) and (max(abs(lf - lf_last)) < tolL):
            converged = True
            print('Density bins converged with an accuray: ', max(abs(delta - z_bin_deltas)))
            print("LF bins converged to an accuracy: ", max(abs(lf - lf_last)))

        
        volbin, pvolbin = tabulate_bin_pvolumes(z_bins, 0.0001, 0.1, par.a, omega0=0.3)
        


"""
def rancat_jswml(ncat, cat, nmult, iseed, nran, rancat_out, nlf, magbin, lf_out, nzbin, zbins_out,
                 delta_out, zref, survey, omega0, lambda0, prior, par_it_out,
                 nitermax, search_for_max, iv):

    # Constants
    tolM = 1e-4
    tolL = 1e-6
    tolD = 1e-4
    tolE = 1e-3
    tolR = 1e-4
    n_srch = 1000

    # Redshift bins
    #dzbin = (survey.zmax - survey.zmin) / nzbin
    #zbin = [BinnedZ(z=survey.zmin + (i + 0.5) * dzbin) for i in range(nzbin)]
    #zbins_out[:] = [z.z for z in zbin]
    zbins = np.linspace(survey.zmin, survey.zmax, nzbin)

    # Check magnitude bins uniformity, 
    dmag = (magbin[nlf - 1] - magbin[0]) / (nlf - 1)
    for i in range(nlf):
        expected = magbin[0] + i * dmag
        if abs(magbin[i] - expected) > tolM * dmag:
            print(f"Warning: Resetting magnitude bin {i} to uniform spacing")
            magbin[i] = expected

    par = ModelParameters()

    if iv >= 3:
        print("Tabulating P-evolution weighted volume of each redshift bin")
    
    vol_bins, p_vol_bins = tabulate_bin_pvolumes(zbins, survey.solid_angle, zref, par.p_post)

    # Compute N(z)
    if iv >= 3:
        print("Computing N(z)")
    for g in cat:
        g.weight = 1.0
    wcic_histogram()

    # Allocate output arrays
    nrancat = int(1 + ncat * nmult * 1.1)
    rancat = [Galaxy(z=0.0) for _ in range(nrancat)]
    par_it = [ModelParameters() for _ in range(nitermax)]

    # Initialization
    par_prev = ModelParameters(u=1e30, a=1e30, mu=1e30)
    par_last = ModelParameters()
    lf = np.zeros(nlf)
    lf_last = np.zeros(nlf)
    delta = np.ones(nzbin)

    convergedPE = False
    converged = False
    convergedu = False
    ibranch = 0
    u_srch = np.zeros(n_srch)
    logpmax_srch = np.zeros(n_srch)
    u_srch[0] = par.u
    i_srch = 0
    iter = 1

    # Main iteration loop (skeleton)
    while True:
        i_srch += 1
        if i_srch > n_srch:
            raise RuntimeError("ERROR: rancat_jswml() n_srch too small")

        nitermin = min(2 + iter, nitermax)
        while not converged and iter <= nitermax:
            # Convergence checks and updates go here
            if iv >= 3:
                print(f"Iteration {iter}: Retabulating distance relation and P-weighted volbins")

            bins, survey_area_sr, zref, Pparam, omega0=0.3

            tabulate_bin_pvolumes(survey, nzbin, [z.z for z in zbin], zbin, zref)

            derived_gal_props(par.u, zref, survey, ncat, cat)

            lf_last[:] = lf[:]
            estimate_lf(ncat, cat, nzbin, zbin, dzbin, delta, nlf, magbin, lf, par.u)

            normalize_weights(cat, ncat, par)

            par.p_post = logp_post(par, prior, zbin, delta, nzbin, lf, nlf)
            par_it[iter - 1] = par
            if iv >= 2:
                print(f"iter={iter} Log(P_post)={par.p_post:.6f}")

            estimate_nofz(ncat, cat, par.mu, zbin, dzbin)

            ntot = sum(z.n for z in zbin)
            ntot_ran = sum(z.n_ran for z in zbin)
            if iv >= 3:
                print(f"Renormalization factor for rancat: {ntot / ntot_ran:.3f}")

            for i in range(nzbin):
                zbin[i].delta = delta[i]
            delta, delta_raw, sigma = estimate_delta(nzbin, zbin, prior.fourPiJ3)

            ntot_model = sum(zbin[i].n * delta[i] for i in range(nzbin))
            if iv >= 3:
                print(f"ntot={ntot:.3f}   ntot_model={ntot_model:.3f}")

            # Check convergence (to be implemented)
            # ...

            iter += 1

        if converged or iter > nitermax:
            break

    # Save final results
    rancat_out[:] = rancat
    nran[0] = nrancat
    par_it_out[:] = par_it
    lf_out[:] = lf
    delta_out[:] = delta
"""