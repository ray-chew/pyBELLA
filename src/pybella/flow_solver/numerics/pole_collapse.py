"""Elliptic pole-ring collapse.

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
the fold's sign trap on the cross-metric terms.

Where the bare integers below come from. The solve runs not on the full
node array but on ``node.isc``, the same (lambda, r, phi) array with one
layer stripped from each side; call its shape (Llam, Lr, Lphi). An index
into it is one less than the matching full-node index.

``node.isc`` still keeps ONE ghost layer per side, so along each axis
index 0 and index L-1 are ghosts and the interior runs 1 .. L-2 — which
is where the -2 in the table comes from, not a second stripping. The
three positions the collapse needs:

                                  full node         node.isc
    south pole row                phi = igz         phi = 1
    north pole row                phi = icz-1-igz   phi = Lphi - 2
    lambda = +pi seam duplicate   lam = icx-1-igx   lam = Llam - 2

The seam column is a second copy of a longitude the box already holds, so
``scatter`` writes to it (it must carry the master's value like the rest of
the ring) but ``gather`` skips it, or that wedge would be counted twice.
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
