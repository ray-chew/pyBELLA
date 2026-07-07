"""Elliptic pole-ring collapse (Stage F, F4).

At a lat-lon pole every longitude node of a given radius is ONE physical
point, so leaving them as independent pressure unknowns makes the discrete
Helmholtz system rank-deficient / inconsistent there. The fix is a
master-node Galerkin embedding: one MASTER unknown per (radius, hemisphere)
pole ring, with

    R  = scatter : master -> every longitude node of the ring
    R^T = gather : ring-sum the equations back into the master

The wrapped operator ``gather o A o scatter`` is the Galerkin projection
``R^T A R``; the reduced system is embedded in the full-size solve vector
with the non-master ring entries pinned to zero (they stay zero because the
rhs and every operator output has them zeroed), so the BiCGSTAB plumbing
and the reshape are untouched. The pole rows are assembled ONE-SIDED
(``lap3D`` already treats the non-periodic phi axis as a wall: it zeroes the
beyond-pole coefficient slabs and does no ghost exchange), which sidesteps
the fold's sign trap on the cross-metric terms — see
dev_notes/sphere_poles_plan.md (F4).

The solve box is ``node.isc`` = the full node array minus one layer per
side, so box index ``b`` maps to full-node index ``b + 1``: the pole phi
rows (full-node phi = igz and icz-1-igz, i.e. |phi| = pi/2) land at box-phi
``1`` and ``Lphi - 2``, and the longitude +pi seam duplicate lands at box
``Llam - 2`` (excluded from the ring-sum so each physical node counts once).
"""

import numpy as np

from ...utils import options as opts


def pole_axis_present(ud):
    bdry = getattr(ud, "bdry_type", None)
    return bdry is not None and any(bt == opts.BdryType.POLE for bt in bdry)


class PoleCollapse:
    """Precomputed scatter/gather index maps on the flat 3D solve vector."""

    def __init__(self, node):
        self.shape = tuple(int(s) for s in node.isc)
        Ll, Lr, Lp = self.shape
        self.n = Ll * Lr * Lp
        pole_phi = (1, Lp - 2)  # south (|phi|=pi/2) and north pole node rows
        r_int = range(1, Lr - 1)  # interior radius unknowns
        lam_all = range(1, Ll - 1)  # every interior longitude (incl. seam)
        lam_unique = range(1, Ll - 2)  # exclude the +pi seam duplicate
        master_lam = 1

        def flat(il, ir, ip):
            return (il * Lr + ir) * Lp + ip

        ring_all, uniq_mem, uniq_master = [], [], []
        scatter_src = np.arange(self.n)
        for ip in pole_phi:
            for ir in r_int:
                m = flat(master_lam, ir, ip)
                for il in lam_all:
                    idx = flat(il, ir, ip)
                    ring_all.append(idx)
                    scatter_src[idx] = m  # broadcast master over the ring
                for il in lam_unique:
                    uniq_mem.append(flat(il, ir, ip))
                    uniq_master.append(m)
        self.ring_all = np.asarray(ring_all, dtype=np.intp)
        self.uniq_mem = np.asarray(uniq_mem, dtype=np.intp)
        self.uniq_master = np.asarray(uniq_master, dtype=np.intp)
        self.scatter_src = scatter_src

    def scatter(self, x):
        """Broadcast each master value over its whole longitude ring."""
        return x[self.scatter_src]

    def gather(self, y):
        """Ring-sum the unique nodes into the master; zero the rest."""
        out = y.copy()
        out[self.ring_all] = 0.0
        np.add.at(out, self.uniq_master, y[self.uniq_mem])
        return out


def get(node):
    """Cached collapse for a node grid (depends only on node.isc)."""
    c = getattr(node, "_pole_collapse", None)
    if c is None:
        c = PoleCollapse(node)
        node._pole_collapse = c
    return c
