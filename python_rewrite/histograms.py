"""
Recreation of the histograms fortran module
"""

import numpy as np

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
