from enum import Enum # ! Version > Python 3.4

class BdryType(Enum):
    """
    An enumeration class that defines the accepted boundary condition types.
    """
    
    WALL = 'symmetric'
    PERIODIC = 'wrap'
    RAYLEIGH = 'radiation'