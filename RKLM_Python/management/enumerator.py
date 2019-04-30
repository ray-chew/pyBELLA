from enum import Enum # ! Version > Python 3.4

class TimeIntegrator(Enum):
    OP_SPLIT = 0
    OP_SPLIT_MD_UPDATE = 1
    HUEN = 2
    EXPL_MIDPT = 3
    RK3_SKAMA = 4
    RK3_TEST = 5
    SI_MIDPT = 6
    STRANG = 7

class MolecularTransport(Enum):
    FULL_MOLECULAR_TRANSPORT = 0
    STRAKA_DIFFUSION_MODEL = 1
    NO_MOLECULAR_TRANSPORT = 2
