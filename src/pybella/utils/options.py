from enum import Enum  # ! Version > Python 3.4


class LimiterType(Enum):
    NONE = 0
    # MINMOD = 1
    # VANLEER = 2
    # VANLEERSmooth = 3
    # SUPERBEE = 4
    # MONOTONIZED_CENTRAL = 5
    # SWEBY_MUNZ = 6
    # RUPE = 7
    # NO_SLOPE = 8
    # NUMBER_OF_LIMITER = 9


class BdryType(Enum):
    """
    An enumeration class that defines the accepted boundary condition types.
    """

    WALL = "symmetric"
    PERIODIC = "wrap"
    RAYLEIGH = "radiation"
    #: lat-lon coordinate-singularity boundary (the pole of a spherical
    #: map). The pole is neither a wall nor periodic: a meridional flow
    #: crossing the pole at longitude lambda re-emerges at lambda + pi.
    POLE = "pole"
