"""
Putting GAMA catalog into form that fortran can run on.
"""

import pandas as pd
import astropy.units as u
from astropy.units import Quantity
import numpy as np

from geom_calcs import calculate_area_of_rectangular_region


def main(
    z_min: float, z_max: float, mag_lim: float, names: list[str], areas=list[Quantity]
) -> None:
    """
    Main logic of the script. Loops through the gama catalogues and puts them in the correct
    format for the fortran code.
    """
    for area, name in zip(areas, names):
        df = pd.read_csv(
            f"~/Desktop/GAMA_paper_plotter/gama_galaxy_catalogs/{name}_galaxies.dat",
            sep="\s+",
        )
        df = df[df["Rpetro"] < mag_lim]
        with open(f"cole_fortran/TestData/{name}_catalogue.txt", "w", encoding="utf-8") as file:
            file.write("# Selection criteria:\n")
            file.write(f"# zmin=  {z_min:.4f} zmax=  {z_max:.4f}\n")
            file.write(f"# faint maglimit= {mag_lim:.4f}\n")
            file.write(f"# solid angle=  {area.value:.4f}  square degrees\n")
            file.write(f"# Catalogue: ncat=       {int(len(df))}\n")
            file.write("mag     z      id\n")
            mags = np.array(df["Rpetro"])
            redshifts = np.array(df["Z"])
            ids = np.arange(len(mags)) + 1
            for mag, redshift, _id in zip(mags, redshifts, ids):
                file.write(f"{mag} {redshift} {_id}\n")


if __name__ == "__main__":
    ZMIN = 0
    ZMAX = 0.6
    MAG_LIM = 19.65
    g09_area = calculate_area_of_rectangular_region(
        129 * u.deg, 141 * u.deg, -2 * u.deg, 3 * u.deg
    )
    g12_area = calculate_area_of_rectangular_region(
        174 * u.deg, 186 * u.deg, -3 * u.deg, 2 * u.deg
    )
    g15_area = calculate_area_of_rectangular_region(
        211.5 * u.deg, 223.5 * u.deg, -2 * u.deg, 3 * u.deg
    )
    g23_area = calculate_area_of_rectangular_region(
        339 * u.deg, 351 * u.deg, -35 * u.deg, -30 * u.deg
    )

    gama_names = ["g09", "g12", "g15", "g23"]
    gama_areas = [g09_area, g12_area, g15_area, g23_area]
    main(ZMIN, ZMAX, MAG_LIM, gama_names, gama_areas)
