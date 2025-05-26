"""
Recreation of the histograms fortran module
"""

import numpy as np

def whist(values, weights, bin_edges):
    """
    Simple weighted historgam
    """
    hist, _ = np.histogram(values, bins=bin_edges, weights=weights)
    return hist

def wcic_histogram(values, weights, bins):
    """
    Construct a cloud-in-cell (CIC) weighted histogram assuming uniformly spaced bins.

    Parameters:
        values (ndarray): Array of values to bin (length N)
        weights (ndarray): Weights associated with each value (length N)
        bins (ndarray): Bin centers (length B); assumed uniform spacing

    Returns:
        hist (ndarray): Weighted histogram (length B)
    """
    values = np.asarray(values)
    weights = np.asarray(weights)
    bins = np.asarray(bins)

    nbin = len(bins)
    hist = np.zeros(nbin)
    inv_dbin = (nbin - 1.0) / (bins[-1] - bins[0])  # Inverse bin spacing

    for v, w in zip(values, weights):
        ri = 1.0 + (v - bins[0]) * inv_dbin  # fractional bin index
        ibin = int(ri + 0.5)

        if ri <= ibin:
            ileft = ibin - 1
        else:
            ileft = ibin

        iright = ileft + 1
        wleft = iright - ri

        # Boundaries: underflow and overflow go to first or last bin
        if ileft < 0:
            ileft = 0
            iright = 0
        if iright >= nbin:
            iright = nbin - 1
            ileft = nbin - 1

        hist[ileft] += w * wleft
        hist[iright] += w * (1.0 - wleft)

    return hist


def wcic_histogram_deriv(values, dvalues_dx, weights, bins):
    """
    Construct a Cloud-In-Cell (CIC) weighted histogram and its derivative 
    with respect to a perturbation in the values.

    Parameters
    ----------
    values : ndarray
        Data values to be histogrammed (length N).
    dvalues_dx : ndarray
        Derivative of each value with respect to x (length N).
    weights : ndarray
        Weights associated with each value (length N).
    bins : ndarray
        Uniformly spaced bin centers (length B).

    Returns
    -------
    hist : ndarray
        Weighted histogram (length B).
    dhist_dx : ndarray
        Derivative of histogram with respect to x (length B).
    """
    values = np.asarray(values)
    dvalues_dx = np.asarray(dvalues_dx)
    weights = np.asarray(weights)
    bins = np.asarray(bins)

    nbin = len(bins)
    hist = np.zeros(nbin)
    dhist_dx = np.zeros(nbin)

    inv_dbin = (nbin - 1.0) / (bins[-1] - bins[0])

    for v, dv_dx, w in zip(values, dvalues_dx, weights):
        ri = 1.0 + (v - bins[0]) * inv_dbin  # fractional bin index
        dri_dx = dv_dx * inv_dbin            # derivative of ri w.r.t x
        dwleft_dx = -dri_dx                  # shift right reduces weight on left bin

        ibin = int(ri + 0.5)
        ileft = ibin - 1 if ri <= ibin else ibin
        iright = ileft + 1
        wleft = iright - ri

        # Clamp to bounds
        if ileft < 0:
            ileft = iright = 0
        elif iright >= nbin:
            ileft = iright = nbin - 1

        # Assign weights to left/right bins
        hist[ileft] += w * wleft
        hist[iright] += w * (1.0 - wleft)

        # Handle 100% right-bin weight edge case
        if wleft == 0.0 and dwleft_dx < 0:
            ileft += 1
            iright += 1
            if iright >= nbin:
                iright = nbin - 1
            if ileft >= nbin:
                ileft = nbin - 1

        # Accumulate derivatives
        dhist_dx[ileft] += w * dwleft_dx
        dhist_dx[iright] += w * (-dwleft_dx)

    return hist, dhist_dx
