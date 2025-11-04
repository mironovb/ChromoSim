from __future__ import annotations
from collections.abc import Mapping
import numpy as np
from scipy.sparse import diags
from .sparsity import sparsity_pattern_multi

def _get(par, key):
    """Get a parameter from either a Mapping (dict) or an object with attrs."""
    if isinstance(par, Mapping):
        return par[key]
    return getattr(par, key)

def _set(par, key, value):
    """Set a parameter on either a Mapping (dict) or an object with attrs."""
    if isinstance(par, Mapping):
        par[key] = value
    else:
        setattr(par, key, value)

def build_matrices(par):
    """
    Assemble finite-volume derivative matrices in z and attach to `par`.
    Requires: par.N (int), par.dz (float).
    Adds:
      par.D1  -- first derivative (backward/upwind)  [N x N] sparse
      par.D2  -- second derivative (central)         [N x N] sparse
    """
    N  = int(_get(par, "N"))
    dz = float(_get(par, "dz"))

    e = np.ones(N)

    # 1) Backward/upwind first derivative
    D1 = diags([-e, e], offsets=[-1, 0], shape=(N, N), format="lil") / dz
    D1[0, :] = 0  # inlet row handled via boundary condition
    D1 = D1.tocsc()

    # 2) Central second derivative
    D2 = diags([e, -2.0 * e, e], offsets=[-1, 0, 1], shape=(N, N), format="lil") / (dz**2)
    D2[0, 0]   = -2.0 / (dz**2)  # ghost for inlet Dirichlet
    D2[-1, -1] =  1.0 / (dz**2)  # Danckwerts outlet
    D2 = D2.tocsc()

    _set(par, "D1", D1)
    _set(par, "D2", D2)
    return par

def attach_jacobian_sparsity(par):
    """
    Attach Jacobian sparsity pattern for multi-species model.
    Requires: par.N, par.Nsp.
    Adds: par.Jpattern
    """
    N   = int(_get(par, "N"))
    Nsp = int(_get(par, "Nsp"))
    Jpat = sparsity_pattern_multi(N, Nsp)
    _set(par, "Jpattern", Jpat)
    return par
