from enum import Enum # ! Version > Python 3.4

# class RecoveryOrder(Enum):
#     FIRST = 0
#     SECOND = 1

# class TimeIntegrator(Enum):
#     OP_SPLIT = 0
#     OP_SPLIT_MD_UPDATE = 1
#     HUEN = 2
#     EXPL_MIDPT = 3
#     RK3_SKAMA = 4
#     RK3_TEST = 5
#     SI_MIDPT = 6
#     STRANG = 7

# class HillShapes(Enum):
#     SCHULTOW = 0
#     AGNESI = 1

# class MolecularTransport(Enum):
#     FULL_MOLECULAR_TRANSPORT = 0
#     STRAKA_DIFFUSION_MODEL = 1
#     NO_MOLECULAR_TRANSPORT = 2

# class BottomBC(Enum):
#     ZERO_ORDER_EXTRAPOL = 0
#     BOTTOM_BC_DEFAULT = 1

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
    
    WALL = 'symmetric'
    PERIODIC = 'wrap'
    RAYLEIGH = 'radiation'