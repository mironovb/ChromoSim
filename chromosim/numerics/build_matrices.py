# chromosim/numerics/build_matrices.py
"""
Finite-volume derivative matrices in z for 1D packed-column models.

This is a direct Python/NumPy/SciPy port of the MATLAB build_matrices.m:

- First derivative (D1): backward/upwind stencil, inlet row is zeroed
  because the inlet boundary condition is handled explicitly in the ODE.
- Second derivative (D2): central stencil with:
    * inlet ghost-cell treatment: D2[0,0] = -2/dz^2 (Dirichlet-type)
    * Danckwerts-style outlet:    D2[-1,-1] =  1/dz^2

Returns the same 'par' object/dict with .D1 and .D2 as CSR sparse matrices.
"""

from __future__ import annotations
import numpy as np
from scipy.sparse import diags

def _get(par, key):
    """Support both dicts and objects with attributes."""
    try:
        return getattr(par, key)
    except AttributeError:
        return par[key]

def _set(par, key, value):
    """Support both dicts and objects with attributes."""
    try:
        setattr(par, key, value)
    except AttributeError:
        par[key] = value

def build_matrices(par):
    """
    Build derivative matrices and attach them to `par`.

    Parameters
    ----------
    par : dict or object with attributes
        Required fields:
          - N  : int   number of axial cells
          - dz : float cell length (m)

    Side effects
    ------------
    Sets:
      par.D1 : (N x N) scipy.sparse.csr_matrix
      par.D2 : (N x N) scipy.sparse.csr_matrix

    Returns
    -------
    par : same object/dict, for chaining.
    """
    N  = int(_get(par, "N"))
    dz = float(_get(par, "dz"))

    e = np.ones(N)

    # 1) First derivative: backward/upwind
    D1 = diags([-e, e], offsets=[-1, 0], shape=(N, N), format="lil") / dz
    D1[0, :] = 0.0  # inlet row is handled via the boundary condition
    D1 = D1.tocsr()

    # 2) Second derivative: central
    D2 = diags([e, -2*e, e], offsets=[-1, 0, 1], shape=(N, N), format="lil") / dz**2
    D2[0, 0]  = -2.0 / dz**2   # inlet ghost cell (Dirichlet)
    D2[-1, -1] =  1.0 / dz**2  # Danckwerts-type outlet
    D2 = D2.tocsr()

    _set(par, "D1", D1)
    _set(par, "D2", D2)
    return par

if __name__ == "__main__":
    # Minimal smoke test
    from types import SimpleNamespace
    par = SimpleNamespace(N=5, dz=0.01)
    par = build_matrices(par)
    print("D1:\n", par.D1.toarray())
    print("D2:\n", par.D2.toarray())
