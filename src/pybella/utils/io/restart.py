"""Simulation restart and output-filename generation."""

import h5py


class read_input(object):
    def __init__(self, fn, path):
        self.fn = fn
        self.path = path

    def get_data(self, Sol, npf, time_tag, half=False):
        file = h5py.File(self.path + "/" + self.fn, "r")

        if half:
            half_tag = "_half"
        else:
            half_tag = ""

        Sol.rho[...] = file["rho" + half_tag]["rho" + half_tag + "_" + time_tag][:]
        Sol.rhou[...] = file["rhou" + half_tag]["rhou" + half_tag + "_" + time_tag][:]
        Sol.rhov[...] = file["rhov" + half_tag]["rhov" + half_tag + "_" + time_tag][:]
        Sol.rhow[...] = file["rhow" + half_tag]["rhow" + half_tag + "_" + time_tag][:]
        Sol.rhoY[...] = file["rhoY" + half_tag]["rhoY" + half_tag + "_" + time_tag][:]

        npf.p2_nodes[...] = file["p2_nodes" + half_tag][
            "p2_nodes" + half_tag + "_" + time_tag
        ][:]

        file.close()


def sim_restart(path, name, elem, node, ud, Sol, npf, restart_touts):
    """
    Function to restart simulation from a saved file. Dataset has to be structured in the same way as the output format of this file.

    Parameters
    ----------
    path : str
        path to the hdf5 file for simulation restart.

    """
    file = h5py.File(str(path), "r")

    Sol_data = ["rho", "rhou", "rhov", "rhow", "rhoX", "rhoY"]
    npf_data = ["p2_nodes"]

    for data in Sol_data:
        value = file[data][data + name][:]

        if hasattr(Sol, data):
            shp = getattr(Sol, data).shape
            setattr(Sol, data, value)
            assert getattr(Sol, data).shape == shp
        else:
            assert 0, "Sol attribute mismatch"

    for data in npf_data:
        value = file[data][data + name][:]
        if hasattr(npf, data):
            shp = getattr(npf, data).shape
            setattr(npf, data, value)
            assert getattr(npf, data).shape == shp
        else:
            assert 0, "npf attribute mismatch"

    t = restart_touts

    ud.output_suffix = "_%i_%i_%i_%.1f_%s" % (
        ud.inx - 1,
        ud.iny - 1,
        ud.inz - 1,
        t[-1],
        ud.aux,
    )

    file.close()
    return Sol, npf, t


def fn_gen(ud, dap, N):
    suffix = ""
    suffix += "_%i" % (ud.inx - 1)
    suffix += "_%i" % (ud.iny - 1)
    if ud.iny == 2:
        suffix += "_%i" % (ud.inz - 1)
    suffix += "_%.6f" % ud.tout[-1]
    suffix = "_ensemble=%i%s" % (N, suffix)

    if len(dap.da_times) > 0 and N > 1:
        suffix += "_wda"
        if dap.da_type == "rloc":
            suffix += "wloc"
        if len(dap.obs_attributes) < 5:
            for attr in dap.obs_attributes:
                suffix += "_%s" % attr
        else:
            suffix += "_all"

    if ud.aux is not None:
        suffix += "_" + ud.aux

    if ud.initial_blending:
        bw = int(ud.blending_weight * 16)
        suffix += "_ib-%i" % bw

    if ud.continuous_blending == True:
        suffix += "_cont_blend"
        suffix += "_fs=%i_ts=%i" % (ud.no_of_pi_initial, ud.no_of_pi_transition)

    return suffix
