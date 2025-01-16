"""
Simple Gaussian smoothing of redshifts and fitting a sqewed normal pdf
"""

import numpy as np
from scipy.stats import skewnorm, gamma
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt

# Example redshift array
redshifts = np.loadtxt('cut_9.dat', usecols=(-2), skiprows=1)


# 1. Create a histogram
bins = np.arange(np.min(redshifts), np.max(redshifts), 0.01)  # Define bins for the histogram
hist, bin_edges = np.histogram(redshifts, bins=bins, density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2  # Get the center of each bin


# 2. Smooth the histogram using a Gaussian filter
smoothed_hist = gaussian_filter1d(hist, sigma=4)

# 3. Fit a skewed normal distribution
# Fit using scipy.stats.skewnorm
params = skewnorm.fit(redshifts)  # Returns shape (alpha), location (xi), scale (omega)
params = gamma.fit(redshifts)
alpha, loc, scale = params

# Generate fitted values
x_fit = np.linspace(0, 1, 500)  # Smooth range for the fit
fitted_values = gamma.pdf(x_fit, alpha, loc, scale)

# 4. Plot the results
gama_rand = np.loadtxt('randoms_unpublished.csv', skiprows=1)
plt.figure(figsize=(8, 5))
plt.hist(redshifts, bins=bins, density=True, alpha=0.5, label="Histogram")
plt.hist(gama_rand, bins=bins, density=True, histtype='step', lw=3, color='r')
plt.plot(bin_centers, smoothed_hist, label="Smoothed Histogram", lw=2)
plt.plot(x_fit, fitted_values, label=f"Skewed Normal Fit\nalpha={alpha:.2f}, loc={loc:.2f}, scale={scale:.2f}", lw=2)
plt.xlabel("Redshift")
plt.ylabel("Density")
plt.legend()
plt.show()
