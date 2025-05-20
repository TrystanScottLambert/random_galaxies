"""
ASCII version of rancat_example.
Inputs and outputs ASCII data files rather than hdf5.

This program reads an input catalogue with corresponding selection parameters,
finds the JSWML Luminosity Function (LF) and generates a corresponding random catalogue.
The random catalogue, the LF and various auxiliary parameters, together
with a copy of the original catalogue are written to a set of ascii files.
"""

import numpy as np
from astropy.cosmology import FlatLambdaCDM
import pandas as pd

# Import modules that we'll implement later
from datatypes import SurveySpec, PriorParameters, ModelParameters, Galaxy
from galaxy_datatype import GalaxyType
from rancat_subs import rancat_subs
from jswml import rancat_jswml


def main():
    # Cosmological parameters
    omega0 = 0.3
    lambda0 = 0.7
    cosmo = FlatLambdaCDM(H0=70, Om0=omega0, Ode0=lambda0)
    
    # Set the priors
    prior = PriorParameters()
    prior.fourPiJ3 = 5000.0  # Galaxy clustering parameter used in density fluctuation prior
    prior.u = 0.01  # Sigma of Gaussian prior on luminosity evolution parameter
    prior.a = 0.01  # Sigma of Gaussian prior on density evolution parameter
    prior.spline = True  # Use B-spline approximation instead of pure Gaussian
    
    # Survey parameters
    survey = SurveySpec()
    
    # Program parameters
    nmult = 400  # Factor by which random catalogue is larger than original
    nzbin = 40   # Number of redshift bins between min and max redshift
    nlf = 25     # Number of absolute magnitude bins in the LF
    dmag = 0.5   # Magnitude bin spacing
    nitermax = 400  # Maximum number of iterations
    iseed = -5875  # Random seed for generating random catalogue
    search_for_max = False  # Set to True for the method described in the paper
    verbosity = 2  # Controls amount of info written during execution
    zref = 0.0  # Reference redshift at which k+e=0
    
    # Initialize magnitude bins for luminosity function
    magbin = np.array([-25.0 + i * dmag for i in range(1, nlf+1)])
    
    # Initialize arrays
    lf = np.zeros(nlf)  # Luminosity function array
    zbin = np.zeros(nzbin)  # Redshift bins
    delta = np.zeros(nzbin)  # Overdensity in redshift bins
    
    # Read the input catalogue
    catfile = "TestData/catalogue.txt"
    print(f"Reading input catalogue from {catfile}")
    
    with open(catfile, 'r') as f:
        # Skip header line
        f.readline()
        
        # Read survey parameters
        line = f.readline().strip().split()
        survey.zmin = float(line[2])
        survey.zmax = float(line[4])
        print(f"zmin={survey.zmin:.4f} zmax={survey.zmax:.4f}")
        
        line = f.readline().strip().split()
        survey.magfaint = float(line[3])
        print(f"mag_faintlimit={survey.magfaint:.4f}")
        
        line = f.readline().strip().split()
        survey.solid_angle = float(line[3])
        print(f"solid angle= {survey.solid_angle:.4f} square degrees")
        
        line = f.readline().strip().split()
        ncat = int(line[3])
        
        # Skip another header line
        f.readline()
        
        # Initialize the catalogue array
        cat = [Galaxy() for _ in range(ncat)]
        
        # Read galaxy data
        for i in range(ncat):
            line = f.readline().strip().split()
            cat[i].mag = float(line[0])
            cat[i].z = float(line[1])
            cat[i].id = int(line[2])
    
    print(f"Catalogue read. ncat={ncat} maglimit={survey.magfaint}\n")
    
    # Call the main subroutine to find LF and generate random catalogue
    print("Calling rancat_jswml() to iterate to the ML solution")
    par_it = [ModelParameters() for _ in range(nitermax)]
    
    rancat, nran = rancat_jswml(
        ncat=ncat,
        cat=cat,
        nmult=nmult,
        iseed=iseed,
        nlf=nlf,
        magbin=magbin,
        lf=lf,
        nzbin=nzbin,
        zbin=zbin,
        delta=delta,
        zref=zref,
        survey=survey,
        omega0=omega0,
        lambda0=lambda0,
        prior=prior,
        par_it=par_it,
        nitermax=nitermax,
        search_for_max=search_for_max,
        verbosity=verbosity
    )
    
    # Save results to output files
    rancatfile = "TestData/rancat.txt"
    print(f"Saving results in {rancatfile}")
    
    # Save random catalogue
    with open(rancatfile, 'w') as f:
        f.write("# Parameters:\n")
        f.write(f"# Omega0= {omega0:.4f} Lambda0= {lambda0:.4f}\n")
        f.write(f"# zmin= {survey.zmin:.4f} zmax= {survey.zmax:.4f}\n")
        f.write(f"# faint maglimit= {survey.magfaint:.4f}\n")
        f.write(f"# solid angle= {survey.solid_angle:.4f} square degrees\n")
        f.write(f"# Catalogue: nran= {nran}\n")
        f.write("# mag    z         absmag     dm       v          vmax\n")
        
        for i in range(nran):
            f.write(f"{rancat[i].mag:.4f} {rancat[i].z:.4f} {rancat[i].absmag:.4f} "
                   f"{rancat[i].dm:.4f} {rancat[i].pv:.3f} {rancat[i].pvmax:.3f} {rancat[i].id}\n")
    
    # Save iteration history
    print("Saving iteration history")
    with open("TestData/iterations.txt", 'w') as f:
        f.write(f"# nitermax= {nitermax}\n")
        f.write("# u              a              mu             P_post\n")
        
        for i in range(nitermax):
            f.write(f"{par_it[i].u} {par_it[i].a} {par_it[i].mu} {par_it[i].p_post}\n")
    
    # Save luminosity function
    print("Saving Luminosity Function")
    with open("TestData/LF.txt", 'w') as f:
        f.write(f"# nmag= {nlf}\n")
        f.write("# mag            LF\n")
        
        for i in range(nlf):
            f.write(f"{magbin[i]} {lf[i]}\n")
    
    # Save delta(z)
    print("Saving Delta(z)")
    with open("TestData/delta.txt", 'w') as f:
        f.write(f"# nzbin= {nzbin}\n")
        f.write("# z            delta\n")
        
        for i in range(nlf):  # Note: Using nlf here, same as in original
            f.write(f"{zbin[i]} {delta[i]}\n")
    
    # Save copy of input catalogue with extra parameters
    print("Saving copy of input catalogue but with extra parameters, eg absmag")
    with open("TestData/catalogue_out.txt", 'w') as f:
        f.write("# Selection criteria:\n")
        f.write(f"# zmin= {survey.zmin:.4f} zmax= {survey.zmax:.4f}\n")
        f.write(f"# faint maglimit= {survey.magfaint:.4f}\n")
        f.write(f"# solid angle= {survey.solid_angle:.4f} square degrees\n")
        f.write(f"# Catalogue: ncat={ncat}\n")
        f.write("# mag    z         absmag     dm       v          vmax     id\n")
        
        for i in range(ncat):
            f.write(f"{cat[i].mag:.4f} {cat[i].z:.4f} {cat[i].absmag:.4f} "
                   f"{cat[i].dm:.4f} {cat[i].pv:.3f} {cat[i].pvmax:.3f} {cat[i].id}\n")


if __name__ == "__main__":
    main()