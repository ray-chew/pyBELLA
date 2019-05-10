from enum import Enum # ! Version > Python 3.4

class BdryType(Enum):
    TUNIX = 0
    WALL = 'symmetric'
    INFLOW = 2
    OUTFLOW = 3
    PERIODIC = 'wrap'
    NEUMANN = 5
    DIRICHLET = 6
    OPEN = 7
    SLANTED_WALL = 8
    PRESCRIBED_FLUX = 9