"""
Python module storing all the dataclasses that mimic those in the original fortran code.
"""

from dataclasses import dataclass

@dataclass
class ModelParameters:
    a: float
    u: float
    mu: float
    p_post: float

@dataclass
class SurveySpec:
    zmin: float
    zmax: float
    magfaint: float
    solid_angle: float
    survey_fractional_area: float

@dataclass
class PriorParameters:
    u: float
    a: float
    fourPiJ3: float
    spline: bool = False

@dataclass
class BinnedZ:
    z: float = 0
    n: float = 0
    n_ran: float = 0
    p_vol: float = 0
    vol: float = 0
    delta: float = 0

@dataclass
class Galaxy:
    