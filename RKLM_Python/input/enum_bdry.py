from enum import Enum

class BdryType(Enum):
    TUNIX = 0
    WALL = 1
    INFLOW = 2
    OUTFLOW = 3
    PERIODIC = 4
    NEUMANN = 5
    DIRICHLET = 6
    OPEN = 7
    SLANTED_WALL = 8
    PRESCRIBED_FLUX = 9