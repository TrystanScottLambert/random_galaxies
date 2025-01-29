"""
Putting catalog into form that fortran can run on.
"""

import pandas as pd
from geom_calcs import calculate_area_of_rectangular_region
import astropy.units as u
import numpy as np


if __name__ == '__main__':
    infile = 'cut_9.dat'
    df = pd.read_csv(infile, sep='\s+')
    area = calculate_area_of_rectangular_region(129*u.deg, 141*u.deg, -2*u.deg, 3*u.deg)
    ZMIN=0
    ZMAX=0.6
    MAG_LIM = 19.65
    SOLID_ANGLE = area.value # square degrees

    df = df[df['Rpetro'] < MAG_LIM]
    NCAT = int(len(df))
    with open('catalogue.txt', 'w', encoding='utf-8') as file:
        file.write('# Selection criteria:\n')
        file.write(f'# zmin=  {ZMIN:.4f} zmax=  {ZMAX:.4f}\n')
        file.write(f'# faint maglimit= {MAG_LIM:.4f}\n')
        file.write(f'# solid angle=  {SOLID_ANGLE:.4f}  square degrees\n')
        file.write(f'# Catalogue: ncat=       {NCAT}\n')
        file.write('mag     z      id\n')
        mags = np.array(df['Rpetro'])
        redshifts = np.array(df['Z'])
        ids = np.arange(len(mags)) +1
        for mag, redshift, id in zip(mags, redshifts, ids):
            file.write(f'{mag} {redshift} {id}\n')
