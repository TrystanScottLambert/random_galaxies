"""
user_routines
"""

from dataclasses import dataclass
import numpy as np

# 1) K-correction for GAMA r-band selected galaxies
def kcorr(z):
    """Polynomial fit to GAMA r-band k-correction as a function of redshift."""
    kcorrvals = [0.20848, 1.0226, 0.52366, 3.5902, 2.3843]
    dz = z - 0.2
    return sum(kcorrvals[i] * dz**i for i in range(len(kcorrvals)))

# 2) E-correction
def ecorr(z, u):
    """Luminosity evolution e-correction."""
    return (u - 1.75) * z

# 3) Derivative of e-correction w.r.t. u
def decorr_du(z, u):
    """Derivative of e-correction with respect to u."""
    return z + 0.0 * u  # Force usage of u

# 4) Number density evolution function
def P(z, a):
    """Number density evolution function."""
    return np.exp((a + 0.1) * z)

# 5) Galaxy datatype
@dataclass
class Galaxy:
    # Core properties
    z: float
    mag: float
    
    # Computed by rancat_jswml
    absmag: float = 0.0
    dm: float = 0.0
    zmax: float = 0.0
    
    # Internal quantities (used by rancat_jswml)
    dabsmag_du: float = 0.0
    pvmax: float = 0.0
    pvmaxeff: float = 0.0
    weight: float = 0.0
    pv: float = 0.0
    v: float = 0.0
    
    # User-defined extra properties (example)
    id: int = 0
    # Xmag: float = field(default=0.0)  # Uncomment if you have additional band magnitudes

# 6) User-defined transformations (placeholder)
def user_transformations(clone: Galaxy, original: Galaxy, u: float):
    """
    Transform redshift/distance-dependent user-defined properties.
    This is a placeholder function. You can modify it to update any custom property.
    
    For example, to update X-band apparent magnitudes:
    clone.Xmag = (original.Xmag - original.dm - kcorrX(original.z) - ecorr(original.z, u)
                  + clone.dm + kcorrX(clone.z) + ecorr(clone.z, u))
    """
    # Dummy line to suppress linter complaints
    clone.z = clone.z + 0.0 * u * original.z
