"""
Module for  creating random galaxies that are used in FoFR. This method doesn't use the cole method
instead simply smoothing the data then fitting a skewed-normal pdf which is then used to generate
the underlying distribution. 

But Trystan, if we have the pdf, we don't need to make a distribution and pass that to another
program. Just pass the PDF. 

Shut up. 
"""

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import skewnorm
from scipy.interpolate import interp1d


def smooth_and_sample(redshifts: np.ndarray, smoothing_sigma: float = 4) -> np.ndarray:
    """
    Smooths a redshift distribution using a Gaussian filter and generates a random sample
    from the smoothed PDF as a continuous distribution.

    Parameters:
        redshifts (array-like): Input array of redshift values.
        smoothing_sigma (float): Standard deviation for the Gaussian filter.

    Returns:
        random_sample (ndarray): Random sample generated from the smoothed PDF.
    """
    hist, bin_edges = np.histogram(redshifts, bins="auto", density=True)
    smoothed_pdf = gaussian_filter1d(hist, sigma=smoothing_sigma)

    # Ensure the PDF is normalized
    smoothed_pdf /= np.sum(smoothed_pdf * np.diff(bin_edges))

    # Create a cumulative distribution function (CDF)
    cdf = np.cumsum(smoothed_pdf * np.diff(bin_edges))
    cdf = np.insert(cdf, 0, 0)  # Ensure the CDF starts at 0

    # Interpolate the inverse CDF for sampling
    inverse_cdf = interp1d(cdf, bin_edges, bounds_error=False, fill_value="extrapolate")

    # Generate random samples
    uniform_randoms = np.random.rand(
        len(redshifts)
    )  # Uniformly distributed random numbers
    random_sample = inverse_cdf(uniform_randoms)
    return random_sample


def fit_skewed_normal(
    redshifts: np.ndarray,
    num_samples: int = 1000,
    smoothing_sigma: float = 2,
    smooth: bool = True,
) -> np.ndarray:
    """
    Fits a skewed normal distribution to the data or its smoothed version and generates a random 
    sample. Ensures all generated samples are non-negative.

    Parameters:
        redshifts (array-like): Input array of redshift values.
        num_samples (int): Number of random samples to generate from the fitted distribution.
        smoothing_sigma (float): Standard deviation for the Gaussian filter (used if smooth=True).
        smooth (bool): Whether to smooth the data before fitting the skewed normal distribution.

    Returns:
        random_sample (ndarray): Random sample (non-negative) from the fitted skewed normal 
        distribution.
    """
    if smooth:
        smoothed_data = smooth_and_sample(redshifts, smoothing_sigma=smoothing_sigma)
        data_to_fit = smoothed_data
    else:
        data_to_fit = redshifts

    # Fit a skewed normal distribution to the data
    shape, loc, scale = skewnorm.fit(data_to_fit)

    # Generate random samples and filter out negative values
    random_sample = []
    batch_size = num_samples * 2  # Generate more samples to account for filtering
    while len(random_sample) < num_samples:
        samples = skewnorm.rvs(shape, loc=loc, scale=scale, size=batch_size)
        non_negative_samples = samples[samples >= 0]
        random_sample.extend(non_negative_samples)

    # Truncate to the desired number of samples
    random_sample = np.array(random_sample[:num_samples])
    return random_sample


def generate_randoms(
    redshifts: np.ndarray,
    min_ra: float,
    max_ra: float,
    number_of_samples: int,
    smoothing_sigma: float = 3,
    smooth: bool = False,
) -> np.ndarray:
    """
    Generates a random sample based on the data. Smoothed if the user wishes.
    """
    redshift_sample = redshifts[(redshifts < max_ra) & (redshifts > min_ra)]
    randoms = fit_skewed_normal(
        redshift_sample, number_of_samples, smoothing_sigma, smooth=smooth
    )
    return randoms

def main():
    """
    main function 
    """
    infile =  "../../GAMA_paper_plotter/GAMA_galaxies.dat"
    redshifts = np.loadtxt(infile, usecols=(-1), unpack=True, skiprows=1)
    redshifts = redshifts[(redshifts > 0.01) & (redshifts < 0.6)]

    # generating the gama randoms catalog:
    print("creating randoms")
    randoms = generate_randoms(redshifts, 0.01, 1, int(400 * len(redshifts)), 4, smooth=True)
    print("saving ...")
    outfile = "../../GAMA_paper_plotter/GAMA_randoms_skewed_normal.csv"
    np.savetxt(outfile, randoms, header='z', comments='')

if __name__ == "__main__":
    main()
