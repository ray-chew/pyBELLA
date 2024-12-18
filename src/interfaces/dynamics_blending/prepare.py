from . import schemes

def initialise(sst):
    bld = schemes.Blend(sst.ud)

    sst.interface_params.bld = bld


def init_da_window(sst, tout_old, outer_step):
    # In ensemble case, do blending for each DA window
    if sst.N > 1:
        blend = sst.interface_params.bld if tout_old in sst.da_params.dap.da_times else None
    else:
        blend = sst.interface_params.bld

    # initial blending?
    if sst.ud.initial_blending == True and (outer_step == 0 or outer_step == 1):
        blend = sst.interface_params.bld

    return blend