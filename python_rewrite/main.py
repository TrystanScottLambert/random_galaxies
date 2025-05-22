import numpy as np
from astropy.cosmology import FlatLambdaCDM
import pandas as pd

# Updated import from your new module
from data_classes import SurveySpec, PriorParameters, ModelParameters, Galaxy

# These will be implemented next
from rancat_subs import rancat_subs
from jswml import rancat_jswml


def main():
    # Cosmological parameters
    omega0 = 0.3
    lambda0 = 0.7
    cosmo = FlatLambdaCDM(H0=70, Om0=omega0, Ode0=lambda0)
    
    # Set the priors
    prior = PriorParameters()
    prior.fourPiJ3 = 5000.0
    prior.u = 0.01
    prior.a = 0.01
    prior.spline = True
    
    # Survey parameters
    survey = SurveySpec()
    
    # Program parameters
    nmult = 400
    nzbin = 40
    nlf = 25
    dmag = 0.5
    nitermax = 400
    iseed = -5875
    search_for_max = False
    verbosity = 2
    zref = 0.0
    
    # Initialize magnitude bins for luminosity function
    magbin = np.array([-25.0 + i * dmag for i in range(1, nlf+1)])
    
    # Initialize arrays
    lf = np.zeros(nlf)
    zbin = np.zeros(nzbin)
    delta = np.zeros(nzbin)
    
    # Read the input catalogue
    catfile = "TestData/catalogue.txt"
    print(f"Reading input catalogue from {catfile}")
    
    with open(catfile, 'r', 'utf8') as f:
        f.readline()  # Skip header
        
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
        
        f.readline()  # Skip another header line
        
        cat = [Galaxy() for _ in range(ncat)]
        
        for i in range(ncat):
            line = f.readline().strip().split()
            cat[i].mag = float(line[0])
            cat[i].z = float(line[1])
            cat[i].id = int(line[2])
    
    print(f"Catalogue read. ncat={ncat} maglimit={survey.magfaint}\n")
    
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
    
    rancatfile = "TestData/rancat.txt"
    print(f"Saving results in {rancatfile}")
    
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
    
    print("Saving iteration history")
    with open("TestData/iterations.txt", 'w') as f:
        f.write(f"# nitermax= {nitermax}\n")
        f.write("# u              a              mu             P_post\n")
        for i in range(nitermax):
            f.write(f"{par_it[i].u} {par_it[i].a} {par_it[i].mu} {par_it[i].p_post}\n")
    
    print("Saving Luminosity Function")
    with open("TestData/LF.txt", 'w') as f:
        f.write(f"# nmag= {nlf}\n")
        f.write("# mag            LF\n")
        for i in range(nlf):
            f.write(f"{magbin[i]} {lf[i]}\n")
    
    print("Saving Delta(z)")
    with open("TestData/delta.txt", 'w') as f:
        f.write(f"# nzbin= {nzbin}\n")
        f.write("# z            delta\n")
        for i in range(nlf):  # nlf used here to match original logic
            f.write(f"{zbin[i]} {delta[i]}\n")
    
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
